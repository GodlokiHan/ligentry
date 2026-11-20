"""Build leakage-free training/validation shards for Ligentry pose rescoring.

Usage examples:
  # 1. Generate 200 Ligentry-native redocking complexes (fast, exact pipeline match):
  python scripts/ml/build_dataset.py --mode local-redock --out data/ml_shards --max-complexes 200

  # 2. Convert extracted CrossDocked2020 v1.3 subset:
  python scripts/ml/build_dataset.py --mode crossdocked \
      --crossdocked-root data/crossdocked \
      --types-file data/crossdocked/types/it2_tt_v1.3_10p20n_train0.types \
      --out data/ml_shards --max-complexes 8000
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import logging
import random
import re
import tempfile
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

import ligentry.plugins  # noqa: F401
from ligentry.core.plugin import get_plugin
from ligentry.data.pdb import download_pdb, extract_pocket
from ligentry.plugins._rescore_features import (
    _encode_atom_type,
    featurize_pose_batch,
    parse_pdb_receptor,
    parse_pdbqt_atoms,
)
from ligentry.workflows.benchmark import BENCHMARK_SET
from ligentry.workflows.redock import (
    _crystal_ligand_atoms,
    _fetch_ccd_smiles,
    _mol_from_atoms,
    _pose_models,
    _template_from_smiles,
)

RDLogger.DisableLog("rdApp.warning")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("build_dataset")

# Never allow benchmark complexes (or their exact PDB codes) into train/val
EXCLUDED_PDBS = frozenset(entry[0].upper() for entry in BENCHMARK_SET)

# Curated diverse high-resolution (<2.0 A) PDB complexes for fast local redocking
CURATED_TRAIN_PDBS = [
    "1A4G", "1A4Q", "1B9V", "1C1E", "1C5C", "1D3H", "1E66", "1F0R", "1G9V", "1H1P",
    "1HWI", "1I7Z", "1J4H", "1K3U", "1L2S", "1M2Z", "1N2V", "1O3F", "1P62", "1Q41",
    "1R58", "1S19", "1T46", "1U4D", "1V4S", "1W1P", "1X8X", "1YVF", "1Z95", "2AOG",
    "2BSM", "2C3I", "2D3U", "2E9U", "2FVD", "2G94", "2H4N", "2I0E", "2J7H", "2K3W",
    "2NRY", "2O0U", "2P16", "2QBS", "2R2M", "2V7A", "2WCG", "2X00", "2Y5H", "2Z5X",
    "3B27", "3C4F", "3D4Q", "3E93", "3F9M", "3G0E", "3H9F", "3I25", "3K5V", "3L3M",
    "3M2W", "3N7O", "3O9I", "3P2V", "3Q4Z", "3R88", "3S8O", "3T0W", "3U9Q", "3V2W",
    "3W32", "3ZSO", "4A7I", "4B7P", "4C1U", "4DJV", "4E5W", "4F9W", "4G8Y", "4HNF",
    "4I74", "4J21", "4K77", "4L6R", "4M0Z", "4N6H", "4O2B", "4P6W", "4Q0K", "4R06",
    "4S0V", "4TWP", "4U5J", "4W52", "4X6P", "4Y6D", "4Z0K", "5A5S", "5C2H", "5DWR",
]

MAX_LIG_ATOMS_CHECK = 80
_GNINA_REC_RE = re.compile(r"_\d+\.gninatypes$")
_GNINA_LIG_RE = re.compile(r"_(\d+)\.gninatypes$")


def _assign_split(group_key: str) -> str:
    """Deterministic 85/15 train/val split keyed by receptor/cluster ID."""
    h = int(hashlib.sha256(group_key.encode()).hexdigest(), 16) % 100
    return "val" if h < 15 else "train"


def _rdkit_mol_to_feats(mol: Chem.Mol) -> tuple[np.ndarray, np.ndarray]:
    """Convert an RDKit 3D conformer (heavy atoms) to (coords [N,3], feats [N,16])."""
    mol_no_h = Chem.RemoveHs(mol)
    conf = mol_no_h.GetConformer()
    coords = []
    feats = []
    for atom in mol_no_h.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        el = atom.GetSymbol()
        arom = 1.0 if atom.GetIsAromatic() else 0.0
        h_count = atom.GetTotalNumHs()
        donor = 1.0 if el in ("N", "O", "S") and h_count > 0 else 0.0
        acc = 1.0 if el in ("N", "O", "S", "F") else 0.0
        ad_guess = (
            "A" if (el == "C" and arom)
            else ("OA" if el == "O" else ("NA" if el == "N" else el))
        )
        f = _encode_atom_type(ad_guess, 0.0, is_ligand=True)
        f[11], f[12], f[13] = donor, acc, arom
        coords.append([pos.x, pos.y, pos.z])
        feats.append(f)
    return np.asarray(coords, dtype=np.float32), np.asarray(feats, dtype=np.float32)


def _resolve_crossdocked_rec(root: Path, rec_rel: str) -> Path | None:
    """Map `cluster/5f74_A_rec_0.gninatypes` -> `cluster/5f74_A_rec.pdb` (or `.pdb.gz`)."""
    stem = _GNINA_REC_RE.sub("", rec_rel)
    for ext in (".pdb", ".pdb.gz", ".pdbqt", ".pdbqt.gz", ""):
        cand = root / f"{stem}{ext}"
        if cand.is_file():
            return cand
    return None


def _resolve_crossdocked_lig(root: Path, lig_rel: str) -> tuple[Path | None, int]:
    """Map `cluster/..._docked_5.gninatypes` -> (`cluster/..._docked.sdf.gz`, 5)."""
    if ":" in lig_rel:
        base_rel, idx_s = lig_rel.rsplit(":", 1)
        pose_idx = int(idx_s)
    else:
        m = _GNINA_LIG_RE.search(lig_rel)
        if m:
            pose_idx = int(m.group(1))
            base_rel = _GNINA_LIG_RE.sub("", lig_rel)
        else:
            base_rel = lig_rel
            pose_idx = 0

    for ext in (".sdf.gz", ".sdf", ".mol2.gz", ".mol2", ""):
        cand = root / f"{base_rel}{ext}"
        if cand.is_file():
            return cand, pose_idx
    return None, pose_idx


def build_local_redock_samples(pdb_ids: list[str], seeds: tuple[int, ...] = (42, 7)) -> list[dict]:
    """Run Ligentry's own docking pipeline on non-benchmark PDBs to generate exact training poses."""
    samples: list[dict] = []
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.upper()
        if pdb_id in EXCLUDED_PDBS:
            continue
        with tempfile.TemporaryDirectory(prefix=f"lig_ml_{pdb_id}_") as tmp:
            workdir = Path(tmp)
            try:
                pdb_file = download_pdb(pdb_id, workdir)
                pocket = extract_pocket(pdb_file, workdir)
                resname = pocket.ligand_resname
                crystal_atoms = _crystal_ligand_atoms(pdb_file, resname)
                if not (6 <= len(crystal_atoms) <= 75):
                    continue
                template = _template_from_smiles(_fetch_ccd_smiles(resname))
                if template.GetNumHeavyAtoms() != len(crystal_atoms):
                    continue
                crystal = _mol_from_atoms(crystal_atoms, template)
                Chem.AssignStereochemistryFrom3D(crystal)
                dock_smiles = Chem.MolToSmiles(crystal)

                rec_out = get_plugin("receptor_prep").run(
                    {"receptor_pdb": str(pocket.receptor_pdb)}, workdir
                )
                rec_coords, rec_feats = parse_pdbqt_atoms(
                    Path(rec_out["receptor_pdbqt"]).read_text(), is_ligand=False
                )[0]

                for seed in seeds:
                    lig_out = get_plugin("ligand_prep").run(
                        {"smiles": dock_smiles, "seed": seed}, workdir
                    )
                    poses_path = workdir / f"poses_{seed}.pdbqt"
                    get_plugin("vina").run(
                        {
                            "receptor_pdbqt": rec_out["receptor_pdbqt"],
                            "ligand_pdbqt": lig_out["pdbqt"],
                            "center": pocket.center,
                            "size": pocket.size,
                            "exhaustiveness": 12,
                            "n_poses": 9,
                            "seed": seed,
                            "out_pdbqt": str(poses_path),
                        },
                        workdir,
                    )
                    models, scores = _pose_models(str(poses_path))
                    pose_parsed = parse_pdbqt_atoms(poses_path.read_text(), is_ligand=True)
                    if len(models) < 2 or len(models) != len(pose_parsed):
                        continue

                    rmsds = []
                    for atoms in models:
                        pose_mol = _mol_from_atoms(atoms, template)
                        rmsds.append(float(rdMolAlign.CalcRMS(pose_mol, crystal)))

                    if min(rmsds) > 2.0 or max(rmsds) < 2.5:
                        continue

                    feats = featurize_pose_batch(
                        rec_coords,
                        rec_feats,
                        [p[0] for p in pose_parsed],
                        [p[1] for p in pose_parsed],
                        scores,
                    )
                    feats["rmsds"] = np.asarray(rmsds, dtype=np.float32)[:, None]
                    feats["pkd"] = np.zeros((len(rmsds), 1), dtype=np.float32)
                    feats["has_pkd"] = np.zeros((len(rmsds), 1), dtype=np.float32)
                    feats["split"] = _assign_split(pdb_id)
                    feats["complex_id"] = f"{pdb_id}_{resname}_s{seed}"
                    samples.append(feats)
                    log.info(
                        "Collected %s (%d poses, min RMSD=%.2f A)",
                        feats["complex_id"], len(rmsds), min(rmsds),
                    )
            except Exception as exc:
                log.debug("Skipped %s: %s", pdb_id, exc)
    return samples


def build_crossdocked_samples(
    root: Path, types_file: Path, max_complexes: int
) -> list[dict]:
    """Group CrossDocked2020 .types entries by (receptor_pdb, ligand_sdf) and featurize poses."""
    # (rec_pdb_str, lig_sdf_str) -> list of (pose_idx, label, pkd, rmsd, vina_score)
    grouped: dict[tuple[str, str], list[tuple[int, int, float, float, float]]] = {}
    rec_cache: dict[str, Path | None] = {}
    lig_cache: dict[str, tuple[Path | None, int]] = {}

    log.info("Indexing CrossDocked2020 types file %s ...", types_file)
    with open(types_file) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            label = int(parts[0])
            pkd = float(parts[1])
            rmsd = float(parts[2])
            rec_rel = parts[3]
            lig_rel = parts[4]
            vina_score = 0.0
            if len(parts) >= 6 and parts[5].startswith("#"):
                try:
                    vina_score = float(parts[5][1:])
                except ValueError:
                    vina_score = 0.0

            pdb_code = Path(rec_rel).name[:4].upper()
            if pdb_code in EXCLUDED_PDBS:
                continue

            if rec_rel not in rec_cache:
                rec_cache[rec_rel] = _resolve_crossdocked_rec(root, rec_rel)
            rec_path = rec_cache[rec_rel]
            if rec_path is None:
                continue

            if lig_rel not in lig_cache:
                lig_cache[lig_rel] = _resolve_crossdocked_lig(root, lig_rel)
            lig_path, pose_idx = lig_cache[lig_rel]
            if lig_path is None:
                continue

            grouped.setdefault((str(rec_path), str(lig_path)), []).append(
                (pose_idx, label, pkd, rmsd, vina_score)
            )

    log.info("Matched %d receptor-ligand pose groups on disk.", len(grouped))
    keys = list(grouped.keys())
    random.Random(42).shuffle(keys)

    # Cache parsed receptor coordinates across ligands docked into the same receptor
    parsed_receptors: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    samples: list[dict] = []

    for rec_str, lig_str in keys:
        if len(samples) >= max_complexes:
            break
        entries = sorted(grouped[(rec_str, lig_str)], key=lambda x: x[0])
        if len(entries) < 3:
            continue
        rmsds = [e[3] for e in entries]
        if min(rmsds) > 2.0 or max(rmsds) < 2.5:
            continue

        rec_path = Path(rec_str)
        lig_path = Path(lig_str)
        try:
            if rec_str not in parsed_receptors:
                rec_text = (
                    gzip.open(rec_path, "rt").read()
                    if rec_path.suffix == ".gz"
                    else rec_path.read_text()
                )
                parsed_receptors[rec_str] = (
                    parse_pdbqt_atoms(rec_text, is_ligand=False)[0]
                    if ".pdbqt" in rec_path.name
                    else parse_pdb_receptor(rec_text)
                )
            rec_coords, rec_feats = parsed_receptors[rec_str]

            sdf_bytes = (
                gzip.open(lig_path, "rb").read()
                if lig_path.suffix == ".gz"
                else lig_path.read_bytes()
            )
            suppl = Chem.SDMolSupplier()
            suppl.SetData(sdf_bytes.decode("utf-8", errors="ignore"), removeHs=False)
            mols = [m for m in suppl if m is not None]

            pose_coords, pose_feats, valid_rmsds, valid_pkds, valid_vinas = [], [], [], [], []
            for p_idx, _lbl, pkd, rmsd, v_score in entries:
                if p_idx >= len(mols):
                    continue
                c, f = _rdkit_mol_to_feats(mols[p_idx])
                if len(c) > MAX_LIG_ATOMS_CHECK:
                    continue
                pose_coords.append(c)
                pose_feats.append(f)
                valid_rmsds.append(rmsd)
                valid_pkds.append(abs(pkd))
                valid_vinas.append(v_score)

            if len(pose_coords) < 3 or min(valid_rmsds) > 2.0 or max(valid_rmsds) < 2.5:
                continue

            feats = featurize_pose_batch(
                rec_coords, rec_feats, pose_coords, pose_feats, valid_vinas
            )
            feats["rmsds"] = np.asarray(valid_rmsds, dtype=np.float32)[:, None]
            pkd_arr = np.asarray(valid_pkds, dtype=np.float32)[:, None]
            feats["pkd"] = pkd_arr
            feats["has_pkd"] = (pkd_arr > 0.0).astype(np.float32)
            pocket_cluster = rec_path.parent.name
            feats["split"] = _assign_split(pocket_cluster)
            feats["complex_id"] = f"{pocket_cluster}_{lig_path.stem}"
            samples.append(feats)
            if len(samples) % 250 == 0:
                log.info("Featurized %d / %d complexes ...", len(samples), max_complexes)
        except Exception:
            continue
    return samples


def save_shards(samples: list[dict], out_dir: Path, shard_size: int = 128) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for split in ("train", "val"):
        split_samples = [s for s in samples if s["split"] == split]
        log.info("Saving %d complexes to %s split...", len(split_samples), split)
        for start in range(0, len(split_samples), shard_size):
            chunk = split_samples[start : start + shard_size]
            complex_ptr = [0]
            for s in chunk:
                complex_ptr.append(complex_ptr[-1] + s["node_x"].shape[0])

            shard_path = out_dir / f"{split}_shard_{start // shard_size:04d}.npz"
            np.savez_compressed(
                shard_path,
                node_x=np.concatenate([s["node_x"] for s in chunk], axis=0),
                node_mask=np.concatenate([s["node_mask"] for s in chunk], axis=0),
                edge_index=np.concatenate([s["edge_index"] for s in chunk], axis=0),
                edge_attr=np.concatenate([s["edge_attr"] for s in chunk], axis=0),
                edge_mask=np.concatenate([s["edge_mask"] for s in chunk], axis=0),
                vina_score=np.concatenate([s["vina_score"] for s in chunk], axis=0),
                rmsds=np.concatenate([s["rmsds"] for s in chunk], axis=0),
                pkd=np.concatenate([s["pkd"] for s in chunk], axis=0),
                has_pkd=np.concatenate([s["has_pkd"] for s in chunk], axis=0),
                complex_ptr=np.asarray(complex_ptr, dtype=np.int64),
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Ligentry pose-rescorer dataset shards")
    parser.add_argument("--mode", choices=("local-redock", "crossdocked"), required=True)
    parser.add_argument("--out", type=Path, default=Path("data/ml_shards"))
    parser.add_argument("--max-complexes", type=int, default=100)
    parser.add_argument("--crossdocked-root", type=Path, default=None)
    parser.add_argument("--types-file", type=Path, default=None)
    args = parser.parse_args()

    if args.mode == "local-redock":
        samples = build_local_redock_samples(CURATED_TRAIN_PDBS[: args.max_complexes])
    else:
        if not args.crossdocked_root or not args.types_file:
            parser.error("--crossdocked-root and --types-file are required for crossdocked mode")
        samples = build_crossdocked_samples(args.crossdocked_root, args.types_file, args.max_complexes)

    save_shards(samples, args.out)
    log.info("Done. Total complexes saved: %d", len(samples))


if __name__ == "__main__":
    main()
