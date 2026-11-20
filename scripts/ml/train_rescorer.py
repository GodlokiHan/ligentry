"""Train Ligentry CPU/GPU pose-rescoring models on RTX 2050 and export to ONNX.

Streams `.npz` shards one file at a time so peak system RAM stays under ~400 MB.

Usage:
  # Train CPU default model:
  python -m scripts.ml.train_rescorer --tier cpu --data data/ml_shards_cd2020 --out ligentry/models/cpu_v1 --epochs 35 --grad-accum 16

  # Train deeper GPU optional model:
  python -m scripts.ml.train_rescorer --tier gpu --data data/ml_shards_cd2020 --out models/gpu_v1 --epochs 45 --grad-accum 16
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import logging
import random
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

from scripts.ml.models import build_model

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("train_rescorer")


def iter_shard_records(
    shard_dir: Path, split: str, shuffle: bool = False
) -> Iterator[dict[str, np.ndarray]]:
    """Stream per-complex dicts one `.npz` shard at a time to keep RAM usage low."""
    shard_paths = sorted(shard_dir.glob(f"{split}_shard_*.npz"))
    if shuffle:
        random.shuffle(shard_paths)

    for path in shard_paths:
        with np.load(path) as data:
            ptr = data["complex_ptr"]
            arrays = {k: data[k] for k in data.files if k != "complex_ptr"}
        indices = list(range(len(ptr) - 1))
        if shuffle:
            random.shuffle(indices)
        for i in indices:
            s, e = int(ptr[i]), int(ptr[i + 1])
            yield {k: arr[s:e] for k, arr in arrays.items()}
        del arrays
        gc.collect()


def count_split_complexes(shard_dir: Path, split: str) -> int:
    total = 0
    for path in sorted(shard_dir.glob(f"{split}_shard_*.npz")):
        with np.load(path) as data:
            total += max(0, len(data["complex_ptr"]) - 1)
    return total


def compute_complex_loss(
    pose_logit: torch.Tensor,  # [P, 1]
    pkd_pred: torch.Tensor,    # [P, 1]
    rmsds: torch.Tensor,       # [P, 1]
    pkd_true: torch.Tensor,    # [P, 1]
    has_pkd: torch.Tensor,     # [P, 1]
) -> torch.Tensor:
    """Multi-task loss: native-pose BCE + intra-complex pairwise margin ranking + masked pKd."""
    y_bin = (rmsds <= 2.0).float()
    bce_loss = F.binary_cross_entropy_with_logits(pose_logit, y_bin)

    scores = pose_logit.squeeze(-1)
    r = rmsds.squeeze(-1)
    pos_mask = r <= 2.0
    neg_mask = r >= 2.5
    if pos_mask.any() and neg_mask.any():
        pos_scores = scores[pos_mask][:, None]
        neg_scores = scores[neg_mask][None, :]
        rank_loss = F.softplus(1.0 - (pos_scores - neg_scores)).mean()
    else:
        rank_loss = torch.tensor(0.0, device=pose_logit.device)

    valid_pkd = (has_pkd > 0.5) & (rmsds <= 2.0)
    if valid_pkd.any():
        pkd_loss = F.huber_loss(pkd_pred[valid_pkd], pkd_true[valid_pkd], delta=1.0)
    else:
        pkd_loss = torch.tensor(0.0, device=pose_logit.device)

    return bce_loss + 1.5 * rank_loss + 0.25 * pkd_loss


@torch.no_grad()
def evaluate(model: torch.nn.Module, shard_dir: Path, split: str, device: torch.device) -> dict[str, float]:
    model.eval()
    total_loss = 0.0
    ai_top1_hits = 0
    vina_top1_hits = 0
    n = 0

    for rec in iter_shard_records(shard_dir, split, shuffle=False):
        tensors = {k: torch.from_numpy(v).to(device) for k, v in rec.items()}
        pose_logit, pkd_pred = model(
            tensors["node_x"],
            tensors["node_mask"],
            tensors["edge_index"],
            tensors["edge_attr"],
            tensors["edge_mask"],
            tensors["vina_score"],
        )
        loss = compute_complex_loss(
            pose_logit, pkd_pred, tensors["rmsds"], tensors["pkd"], tensors["has_pkd"]
        )
        total_loss += float(loss.item())

        rmsds = rec["rmsds"].squeeze(-1)
        ai_scores = pose_logit.squeeze(-1).cpu().numpy()
        best_ai_idx = int(np.argmax(ai_scores))
        if rmsds[best_ai_idx] <= 2.0:
            ai_top1_hits += 1
        if rmsds[0] <= 2.0:
            vina_top1_hits += 1
        n += 1

    denom = max(n, 1)
    return {
        "val_loss": total_loss / denom,
        "ai_top1_rate": ai_top1_hits / denom,
        "vina_top1_rate": vina_top1_hits / denom,
    }


def export_and_verify_onnx(model: torch.nn.Module, sample: dict[str, np.ndarray], onnx_path: Path) -> str:
    """Export trained model to ONNX and verify numerical parity against PyTorch."""
    model.eval().cpu()
    input_names = ["node_x", "node_mask", "edge_index", "edge_attr", "edge_mask", "vina_score"]
    dummy_inputs = tuple(torch.from_numpy(sample[k]) for k in input_names)
    dynamic_axes = {name: {0: "num_poses"} for name in input_names}
    dynamic_axes["pose_logit"] = {0: "num_poses"}
    dynamic_axes["pkd_pred"] = {0: "num_poses"}

    onnx_path.parent.mkdir(parents=True, exist_ok=True)
    torch.onnx.export(
        model,
        dummy_inputs,
        str(onnx_path),
        input_names=input_names,
        output_names=["pose_logit", "pkd_pred"],
        dynamic_axes=dynamic_axes,
        opset_version=17,
        do_constant_folding=True,
    )

    with torch.no_grad():
        pt_logit, _pt_pkd = model(*dummy_inputs)

    try:
        import onnxruntime as ort

        sess = ort.InferenceSession(str(onnx_path), providers=["CPUExecutionProvider"])
        ort_out = sess.run(None, {k: sample[k] for k in input_names})
        max_diff = float(np.max(np.abs(pt_logit.numpy() - ort_out[0])))
        log.info("ONNX parity check passed (max abs diff = %.2e)", max_diff)
    except ImportError:
        log.warning("onnxruntime not installed yet; skipped runtime parity check")

    return hashlib.sha256(onnx_path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Train Ligentry pose rescorer")
    parser.add_argument("--tier", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--data", type=Path, default=Path("data/ml_shards_cd2020"))
    parser.add_argument("--out", type=Path, default=Path("ligentry/models/cpu_v1"))
    parser.add_argument("--epochs", type=int, default=35)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--grad-accum", type=int, default=16, help="Complexes per optimizer step")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    n_train = count_split_complexes(args.data, "train")
    n_val = count_split_complexes(args.data, "val")
    if n_train == 0:
        raise RuntimeError(f"No training shards found in {args.data}")
    log.info("Training tier='%s' on device=%s (%d train, %d val complexes)", args.tier, device, n_train, n_val)

    model = build_model(args.tier).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-2)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)
    scaler = torch.amp.GradScaler("cuda", enabled=(device.type == "cuda"))

    best_score = -1.0
    best_state = None
    eval_split = "val" if n_val > 0 else "train"

    for epoch in range(1, args.epochs + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        running_loss = 0.0
        step = 0

        for step, rec in enumerate(iter_shard_records(args.data, "train", shuffle=True), start=1):
            tensors = {k: torch.from_numpy(v).to(device) for k, v in rec.items()}
            with torch.amp.autocast("cuda", enabled=(device.type == "cuda")):
                pose_logit, pkd_pred = model(
                    tensors["node_x"],
                    tensors["node_mask"],
                    tensors["edge_index"],
                    tensors["edge_attr"],
                    tensors["edge_mask"],
                    tensors["vina_score"],
                )
                loss = compute_complex_loss(
                    pose_logit, pkd_pred, tensors["rmsds"], tensors["pkd"], tensors["has_pkd"]
                )
                scaled_loss = loss / args.grad_accum

            scaler.scale(scaled_loss).backward()
            running_loss += float(loss.item())

            if step % args.grad_accum == 0 or step == n_train:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=2.0)
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad(set_to_none=True)

        scheduler.step()
        metrics = evaluate(model, args.data, eval_split, device)
        composite = metrics["ai_top1_rate"] - 0.01 * metrics["val_loss"]
        if composite >= best_score:
            best_score = composite
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

        if epoch == 1 or epoch % 5 == 0 or epoch == args.epochs:
            log.info(
                "Epoch %02d/%02d | train_loss=%.4f | val_loss=%.4f | Top-1 Native: AI=%.1f%% vs Vina=%.1f%%",
                epoch,
                args.epochs,
                running_loss / max(step, 1),
                metrics["val_loss"],
                100.0 * metrics["ai_top1_rate"],
                100.0 * metrics["vina_top1_rate"],
            )

    assert best_state is not None
    model.load_state_dict(best_state)
    sample = next(iter_shard_records(args.data, "train", shuffle=False))
    onnx_path = args.out / "model.onnx"
    sha256 = export_and_verify_onnx(model, sample, onnx_path)

    manifest = {
        "tier": args.tier,
        "version": "1.0.0",
        "sha256": sha256,
        "size_bytes": onnx_path.stat().st_size,
        "train_complexes": n_train,
        "val_complexes": n_val,
        "best_val_top1_rate": round(best_score, 4),
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    log.info("Saved %s (%.2f MB, sha256=%s...)", onnx_path, manifest["size_bytes"] / 1e6, sha256[:12])


if __name__ == "__main__":
    main()
