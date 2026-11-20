"""Shared feature extraction for neural pose rescoring (Qt-free, NumPy/SciPy only).

Extracts a bipartite protein-ligand radial-basis interaction graph from:
  - Runtime: receptor PDBQT + multi-model Vina output PDBQT
  - Training: receptor PDB/PDBQT + multi-conformer docked SDF/PDBQT

Keeping feature extraction in a pure NumPy/SciPy module guarantees zero
train/serve skew and fast CPU inference without requiring PyTorch at runtime.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.spatial.distance import cdist

# Canonical element vocabulary (11 classes: 0..9 + 10=Other)
ELEMENTS = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "METAL")
ELEMENT_TO_IDX = {el: i for i, el in enumerate(ELEMENTS)}
METALS = frozenset({"ZN", "MG", "FE", "CA", "MN", "CU", "CO", "NI", "NA", "K"})

# AutoDock atom types -> (element, is_donor, is_acceptor, is_aromatic)
_AD_PROPS: dict[str, tuple[str, float, float, float]] = {
    "C":  ("C",  0.0, 0.0, 0.0),
    "A":  ("C",  0.0, 0.0, 1.0),
    "N":  ("N",  0.0, 0.0, 0.0),
    "NA": ("N",  0.0, 1.0, 1.0),
    "NS": ("N",  1.0, 0.0, 1.0),
    "OA": ("O",  0.0, 1.0, 0.0),
    "OS": ("O",  0.0, 1.0, 0.0),
    "S":  ("S",  0.0, 0.0, 0.0),
    "SA": ("S",  0.0, 1.0, 0.0),
    "F":  ("F",  0.0, 1.0, 0.0),
    "CL": ("Cl", 0.0, 0.0, 0.0),
    "BR": ("Br", 0.0, 0.0, 0.0),
    "I":  ("I",  0.0, 0.0, 0.0),
    "P":  ("P",  0.0, 0.0, 0.0),
}

# Fixed graph tensor budgets (covers >99.5% of drug-like ligands & 6 A pockets)
MAX_LIG_ATOMS = 80
MAX_REC_ATOMS = 320
MAX_EDGES = 1600
POCKET_CUTOFF = 6.0   # Angstrom: crop receptor atoms within 6 A of any pose
INTER_CUTOFF = 5.5    # Angstrom: protein-ligand interaction edges
INTRA_CUTOFF = 4.0    # Angstrom: intra-ligand spatial/covalent edges
NUM_RBF = 16          # Gaussian radial basis functions over [1.5, 5.5] A

# Node feature dim: 11 element one-hot + [donor, acceptor, aromatic, partial_q, is_ligand] = 16
NODE_FEAT_DIM = 16
# Edge feature dim: 16 RBF + [is_inter, hbond_compat, hydrophobic_pair, charge_prod] = 20
EDGE_FEAT_DIM = 20


def _encode_atom_type(ad_type: str, charge: float, is_ligand: bool) -> np.ndarray:
    """Build 16-d node feature vector from AutoDock atom type + partial charge."""
    t = ad_type.strip().upper()
    if t in _AD_PROPS:
        el, donor, acc, arom = _AD_PROPS[t]
    elif t in METALS:
        el, donor, acc, arom = "METAL", 0.0, 0.0, 0.0
    else:
        el, donor, acc, arom = t.capitalize(), 0.0, 0.0, 0.0

    feat = np.zeros(NODE_FEAT_DIM, dtype=np.float32)
    idx = ELEMENT_TO_IDX.get(el, len(ELEMENTS))
    feat[idx] = 1.0
    feat[11] = donor
    feat[12] = acc
    feat[13] = arom
    feat[14] = float(np.clip(charge, -2.0, 2.0))
    feat[15] = 1.0 if is_ligand else 0.0
    return feat


def parse_pdbqt_atoms(text: str, is_ligand: bool) -> list[tuple[np.ndarray, np.ndarray]]:
    """Parse single- or multi-model PDBQT text into [(coords [N,3], feats [N,16]), ...]."""
    models: list[tuple[np.ndarray, np.ndarray]] = []
    coords: list[list[float]] = []
    feats: list[np.ndarray] = []

    def _flush() -> None:
        if coords:
            models.append((np.asarray(coords, dtype=np.float32), np.asarray(feats, dtype=np.float32)))
            coords.clear()
            feats.clear()

    for line in text.splitlines():
        if line.startswith("MODEL"):
            _flush()
        elif line.startswith("ENDMDL"):
            _flush()
        elif line.startswith(("ATOM", "HETATM")):
            parts = line.split()
            ad_type = parts[-1]
            if ad_type.upper() in ("H", "HD", "HS"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            try:
                q = float(line[70:76])
            except ValueError:
                q = 0.0
            coords.append([x, y, z])
            feats.append(_encode_atom_type(ad_type, q, is_ligand=is_ligand))

    _flush()
    return models


def parse_pdb_receptor(pdb_text: str) -> tuple[np.ndarray, np.ndarray]:
    """Fallback parser for raw PDB files (e.g., CrossDocked2020 receptors)."""
    coords: list[list[float]] = []
    feats: list[np.ndarray] = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[16] not in (" ", "A"):
            continue
        el = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not el:
            el = re.sub(r"[^A-Za-z]", "", line[12:16])[:1].upper()
        if el == "H":
            continue
        atom_name = line[12:16].strip().upper()
        # Heuristic pharmacophore flags for standard protein heavy atoms
        donor = 1.0 if el == "N" and atom_name != "N" else 0.0
        acc = 1.0 if el in ("O", "N") else 0.0
        arom = 1.0 if line[17:20].strip() in ("PHE", "TYR", "TRP", "HIS") and atom_name not in ("N", "CA", "C", "O") else 0.0
        ad_guess = "OA" if el == "O" else ("NA" if el == "N" else ("A" if arom else el))
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        feats.append(_encode_atom_type(ad_guess, 0.0, is_ligand=False))
    return np.asarray(coords, dtype=np.float32), np.asarray(feats, dtype=np.float32)


def _rbf(distances: np.ndarray) -> np.ndarray:
    """Expand distances [E] into 16 Gaussian RBF channels over [1.5, 5.5] A."""
    centers = np.linspace(1.5, 5.5, NUM_RBF, dtype=np.float32)
    gamma = 1.0 / ((centers[1] - centers[0]) ** 2)
    diff = distances[:, None] - centers[None, :]
    return np.exp(-gamma * (diff ** 2))


def featurize_pose_batch(
    rec_coords: np.ndarray,
    rec_feats: np.ndarray,
    pose_coords_list: Sequence[np.ndarray],
    pose_feats_list: Sequence[np.ndarray],
    vina_scores: Sequence[float] | None = None,
) -> dict[str, np.ndarray]:
    """Build fixed-shape padded graph tensors for P poses of a single complex.

    Returns dict of NumPy arrays ready for ONNX Runtime or PyTorch batching:
      node_x:     [P, MAX_NODES, 16] float32
      node_mask:  [P, MAX_NODES]     float32
      edge_index: [P, 2, MAX_EDGES]  int64
      edge_attr:  [P, MAX_EDGES, 20] float32
      edge_mask:  [P, MAX_EDGES]     float32
      vina_score: [P, 1]             float32
    """
    n_poses = len(pose_coords_list)
    max_nodes = MAX_LIG_ATOMS + MAX_REC_ATOMS

    # Crop receptor once using union of all candidate poses for speed & consistency
    all_lig_xyz = np.concatenate(pose_coords_list, axis=0)
    min_dist_to_lig = cdist(rec_coords, all_lig_xyz).min(axis=1)
    pocket_idx = np.where(min_dist_to_lig <= POCKET_CUTOFF)[0]
    if len(pocket_idx) > MAX_REC_ATOMS:
        # Keep closest MAX_REC_ATOMS heavy atoms
        pocket_idx = pocket_idx[np.argsort(min_dist_to_lig[pocket_idx])[:MAX_REC_ATOMS]]

    p_rec_xyz = rec_coords[pocket_idx]
    p_rec_feat = rec_feats[pocket_idx]
    n_rec = len(p_rec_xyz)

    node_x = np.zeros((n_poses, max_nodes, NODE_FEAT_DIM), dtype=np.float32)
    node_mask = np.zeros((n_poses, max_nodes), dtype=np.float32)
    edge_index = np.zeros((n_poses, 2, MAX_EDGES), dtype=np.int64)
    edge_attr = np.zeros((n_poses, MAX_EDGES, EDGE_FEAT_DIM), dtype=np.float32)
    edge_mask = np.zeros((n_poses, MAX_EDGES), dtype=np.float32)
    v_scores = np.zeros((n_poses, 1), dtype=np.float32)

    for p_idx, (l_xyz, l_feat) in enumerate(zip(pose_coords_list, pose_feats_list)):
        n_lig = min(len(l_xyz), MAX_LIG_ATOMS)
        l_xyz = l_xyz[:n_lig]
        l_feat = l_feat[:n_lig]
        n_tot = n_lig + n_rec

        xyz = np.concatenate([l_xyz, p_rec_xyz], axis=0)
        feat = np.concatenate([l_feat, p_rec_feat], axis=0)
        node_x[p_idx, :n_tot] = feat
        node_mask[p_idx, :n_tot] = 1.0
        if vina_scores is not None and p_idx < len(vina_scores):
            v_scores[p_idx, 0] = float(vina_scores[p_idx])

        # Pairwise distances
        dmat = cdist(xyz, xyz)
        is_lig = np.arange(n_tot) < n_lig
        inter_pair = is_lig[:, None] != is_lig[None, :]
        intra_lig_pair = is_lig[:, None] & is_lig[None, :]

        valid = (dmat > 1e-4) & (
            (inter_pair & (dmat <= INTER_CUTOFF)) | (intra_lig_pair & (dmat <= INTRA_CUTOFF))
        )
        src, dst = np.where(valid)
        if len(src) > MAX_EDGES:
            order = np.argsort(dmat[src, dst])[:MAX_EDGES]
            src, dst = src[order], dst[order]

        n_e = len(src)
        if n_e == 0:
            continue

        dists = dmat[src, dst].astype(np.float32)
        rbf = _rbf(dists)

        # Interaction chemistry descriptors
        is_inter = inter_pair[src, dst].astype(np.float32)
        hbond = (
            (feat[src, 11] * feat[dst, 12] > 0.5) | (feat[src, 12] * feat[dst, 11] > 0.5)
        ).astype(np.float32) * is_inter
        hydrophobic = ((feat[src, 0] > 0.5) & (feat[dst, 0] > 0.5)).astype(np.float32) * is_inter
        q_prod = (feat[src, 14] * feat[dst, 14]).astype(np.float32)

        e_feat = np.concatenate(
            [rbf, is_inter[:, None], hbond[:, None], hydrophobic[:, None], q_prod[:, None]],
            axis=1,
        )
        edge_index[p_idx, 0, :n_e] = src
        edge_index[p_idx, 1, :n_e] = dst
        edge_attr[p_idx, :n_e] = e_feat
        edge_mask[p_idx, :n_e] = 1.0

    return {
        "node_x": node_x,
        "node_mask": node_mask,
        "edge_index": edge_index,
        "edge_attr": edge_attr,
        "edge_mask": edge_mask,
        "vina_score": v_scores,
    }


def featurize_from_files(
    receptor_path: str | Path,
    poses_pdbqt_path: str | Path,
    vina_scores: Sequence[float] | None = None,
) -> dict[str, np.ndarray]:
    """Convenience entry point for runtime inference on PDBQT files."""
    rec_text = Path(receptor_path).read_text()
    rec_models = parse_pdbqt_atoms(rec_text, is_ligand=False)
    if rec_models:
        rec_coords, rec_feats = rec_models[0]
    else:
        rec_coords, rec_feats = parse_pdb_receptor(rec_text)

    pose_models = parse_pdbqt_atoms(Path(poses_pdbqt_path).read_text(), is_ligand=True)
    if not pose_models:
        raise ValueError(f"No ligand poses found in {poses_pdbqt_path}")

    pose_coords = [m[0] for m in pose_models]
    pose_feats = [m[1] for m in pose_models]
    return featurize_pose_batch(rec_coords, rec_feats, pose_coords, pose_feats, vina_scores)
