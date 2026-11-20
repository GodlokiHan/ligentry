"""Neural pose rescoring plugin (ONNX Runtime: CPU default + optional GPU pack)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

from ligentry.core.model_manager import get_model_status
from ligentry.core.plugin import Plugin, register
from ligentry.plugins._rescore_features import featurize_from_files

log = logging.getLogger("ligentry.rescorer")

_SESSION_CACHE: dict[tuple[str, str], Any] = {}


def _get_session(model_path: Path, tier: str):
    key = (str(model_path), tier)
    if key not in _SESSION_CACHE:
        import onnxruntime as ort

        providers = (
            ["CUDAExecutionProvider", "CPUExecutionProvider"]
            if tier == "gpu"
            else ["CPUExecutionProvider"]
        )
        opts = ort.SessionOptions()
        opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        _SESSION_CACHE[key] = ort.InferenceSession(
            str(model_path), sess_options=opts, providers=providers
        )
    return _SESSION_CACHE[key]


def _split_pdbqt_models(pdbqt_text: str) -> tuple[list[str], list[float]]:
    """Split a multi-model Vina PDBQT string into per-MODEL blocks and Vina energies."""
    blocks: list[str] = []
    scores: list[float] = []
    cur_lines: list[str] = []
    in_model = False

    for line in pdbqt_text.splitlines():
        if line.startswith("MODEL"):
            in_model = True
            cur_lines = [line]
        elif line.startswith("ENDMDL"):
            cur_lines.append(line)
            blocks.append("\n".join(cur_lines))
            cur_lines = []
            in_model = False
        else:
            if line.startswith("REMARK VINA RESULT"):
                try:
                    scores.append(float(line.split()[3]))
                except (IndexError, ValueError):
                    scores.append(0.0)
            if in_model:
                cur_lines.append(line)

    if cur_lines:
        blocks.append("\n".join(cur_lines))
    return blocks, scores


def _rewrite_ranked_pdbqt(
    blocks: list[str],
    order: list[int],
    ai_probs: np.ndarray,
    pkd_preds: np.ndarray,
) -> str:
    """Re-number MODEL 1..K in AI-ranked order and annotate each pose with REMARK lines."""
    out_models: list[str] = []
    for new_rank, orig_idx in enumerate(order, start=1):
        lines = blocks[orig_idx].splitlines()
        new_lines: list[str] = []
        inserted_remark = False
        ai_remark = (
            f"REMARK LIGENTRY AI SCORE: {ai_probs[orig_idx]:.4f}  "
            f"PRED_PKD: {pkd_preds[orig_idx]:.2f}"
        )
        for line in lines:
            if line.startswith("MODEL"):
                new_lines.append(f"MODEL {new_rank}")
                new_lines.append(ai_remark)
                inserted_remark = True
            else:
                new_lines.append(line)
        if not inserted_remark:
            new_lines.insert(0, ai_remark)
        out_models.append("\n".join(new_lines))
    return "\n".join(out_models) + "\n"


@register
class PoseRescorer(Plugin):
    name = "pose_rescorer"
    description = "Neural pose rescorer (ONNX Runtime, CPU default + optional GPU pack)"

    def run(self, inputs: dict[str, Any], workdir: Path) -> dict[str, Any]:
        receptor_path = Path(inputs["receptor_pdbqt"])
        poses_path = Path(inputs["poses_pdbqt"])
        prefer_gpu = str(inputs.get("device", "auto")).lower() != "cpu"

        status = get_model_status(prefer_gpu=prefer_gpu)
        if status.model_path is None or not status.model_path.is_file():
            raise RuntimeError(
                f"No ONNX rescorer weights found at {status.cpu_pack_dir / 'model.onnx'}"
            )

        raw_text = poses_path.read_text()
        blocks, vina_scores = _split_pdbqt_models(raw_text)
        if not blocks:
            raise ValueError(f"No MODEL blocks found in {poses_path}")

        feats = featurize_from_files(receptor_path, poses_path, vina_scores=vina_scores)
        sess = _get_session(status.model_path, status.active_tier)

        feed = {
            "node_x": feats["node_x"],
            "node_mask": feats["node_mask"],
            "edge_index": feats["edge_index"],
            "edge_attr": feats["edge_attr"],
            "edge_mask": feats["edge_mask"],
            "vina_score": feats["vina_score"],
        }
        pose_logit, pkd_pred = sess.run(None, feed)
        logits = pose_logit.reshape(-1)
        probs = 1.0 / (1.0 + np.exp(-np.clip(logits, -30.0, 30.0)))
        pkds = pkd_pred.reshape(-1)

        order = [int(i) for i in np.argsort(-probs)]
        best_idx = order[0]

        out_path = Path(inputs.get("out_pdbqt", poses_path))
        out_path.write_text(_rewrite_ranked_pdbqt(blocks, order, probs, pkds))

        return {
            "ranked_indices": order,
            "ai_scores": [round(float(probs[i]), 4) for i in order],
            "vina_scores_ranked": [
                round(float(vina_scores[i]), 3) if i < len(vina_scores) else None for i in order
            ],
            "top_pose_index": best_idx,
            "top_ai_score": round(float(probs[best_idx]), 4),
            "top_vina_score": round(float(vina_scores[best_idx]), 3) if vina_scores else None,
            "top_pred_pkd": round(float(pkds[best_idx]), 2),
            "poses_pdbqt": str(out_path),
            "backend": status.active_tier,
        }
