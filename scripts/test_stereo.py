"""Offline checks for stereoisomer enumeration and per-isomer docking.

Run inside the biosim env:  python scripts/test_stereo.py
Uses real RDKit + Meeko; Vina and receptor prep are faked, so no network is needed.
"""
from __future__ import annotations

import csv
import tempfile
import zlib
from pathlib import Path
from types import SimpleNamespace

from rdkit import RDLogger

import biosim.workflows.docking as docking
from biosim.plugins.ligand_prep import enumerate_isomers

RDLogger.DisableLog("rdApp.warning")

failures: list[str] = []


def check(label: str, condition: bool, detail: str = "") -> None:
    print(("PASS  " if condition else "FAIL  ") + label + (f"   {detail}" if detail else ""))
    if not condition:
        failures.append(label)


# ---- enumerate_isomers ---------------------------------------------------------
isos, info = enumerate_isomers("CC(O)C(N)C(=O)O")  # 2 unspecified centres
check("2 centres -> 4 isomers", len(isos) == 4 and not info["sampled"], str(info))
check("isomers are unique", len(set(isos)) == len(isos))

spec = "C[C@H](O)C(=O)O"
isos, info = enumerate_isomers(spec)
check("fully specified -> unchanged, 1 isomer", isos == [spec] and info["n_unspecified"] == 0)

isos, info = enumerate_isomers("CC(O)C(N)C(O)C(N)C(=O)O")  # 4 centres -> 16 possible
check("4 centres capped at 8, sampled", len(isos) == 8 and info["sampled"], str(info))
again, _ = enumerate_isomers("CC(O)C(N)C(O)C(N)C(=O)O")
check("sampling reproducible for same seed", isos == again)
more, info = enumerate_isomers("CC(O)C(N)C(O)C(N)C(=O)O", max_isomers=0)
check("max_isomers=0 -> all 16", len(more) == 16 and not info["sampled"])

isos, info = enumerate_isomers("CC(O)C(N)C(=O)O.[Na+].[Cl-]")
check("salts stripped before counting", len(isos) == 4 and all("." not in s for s in isos))

isos, info = enumerate_isomers("C[C@H](O)C(N)C")  # 1 specified + 1 unspecified
check(
    "specified centre kept, only the other varied",
    len(isos) == 2 and all("@" in s for s in isos) and len(set(isos)) == 2,
    str(isos),
)

isos, info = enumerate_isomers("CC=CC")
check("unspecified C=C -> E and Z", len(isos) == 2, str(isos))

try:
    enumerate_isomers("not a smiles")
    check("invalid SMILES raises ValueError", False)
except ValueError:
    check("invalid SMILES raises ValueError", True)

# ---- run_docking with a fake Vina ----------------------------------------------
vina_calls: list[str] = []


class FakeReceptor:
    def run(self, inputs, workdir):
        return {"receptor_pdbqt": "fake_receptor.pdbqt"}


class FakeVina:
    def run(self, inputs, workdir):
        score = -5.0 - (zlib.crc32(inputs["ligand_pdbqt"].encode()) % 1000) / 100.0
        out = Path(inputs["out_pdbqt"])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(f"SCORE {score}\n")
        vina_calls.append(out.name)
        return {"score": score, "poses_pdbqt": str(out)}


_real_get_plugin = docking.get_plugin


def fake_get_plugin(name):
    if name == "receptor_prep":
        return FakeReceptor()
    if name == "vina":
        return FakeVina()
    return _real_get_plugin(name)  # real ligand_prep (RDKit + Meeko)


with tempfile.TemporaryDirectory() as tmp:
    tmp_path = Path(tmp)
    (tmp_path / "rec.pdb").write_text("END\n")
    docking.get_plugin = fake_get_plugin
    docking.download_pdb = lambda pdb_id, out_dir: tmp_path / "rec.pdb"
    docking.extract_pocket = lambda *a, **k: SimpleNamespace(
        receptor_pdb=tmp_path / "rec.pdb", ligand_resname="LIG",
        center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0),
    )

    ligands = [
        ("two_centres", "CC(O)C(N)C(=O)O"),       # 4 isomers
        ("specified", "C[C@H](O)C(=O)O"),         # 1
        ("salt", "CC(O)C(N)C(=O)O.[Na+].[Cl-]"),  # 4
        ("bad", "not a smiles"),                  # fails
    ]

    vina_calls.clear()
    out = tmp_path / "run_on"
    rows = docking.run_docking("1ABC", ligands, out, enumerate_stereo=True)
    by_name = {r["name"]: r for r in rows}
    check("enumeration on: 4 + 1 + 4 Vina calls", len(vina_calls) == 9, str(len(vina_calls)))
    check("n_isomers recorded", [by_name[n]["n_isomers"] for n in
          ("two_centres", "specified", "salt")] == [4, 1, 4])
    check("bad SMILES recorded as failed", by_name["bad"]["status"].startswith("failed"))
    check("failed ligand ranked last", rows[-1]["name"] == "bad")

    for name in ("two_centres", "salt", "specified"):
        pose = out / "poses" / f"{name}.pdbqt"
        ok = pose.exists() and round(float(pose.read_text().split()[1]), 3) == by_name[name][
            "vina_score_kcal_mol"
        ]
        check(f"{name}: poses/<name>.pdbqt is the best isomer's pose", ok)
    check("no leftover __iso files", not list((out / "poses").glob("*__iso*")))

    with open(out / "results.csv") as fh:
        csv_rows = list(csv.DictReader(fh))
    check("results.csv has isomer columns",
          {"best_isomer_smiles", "n_isomers"} <= set(csv_rows[0]))

    vina_calls.clear()
    off = tmp_path / "run_off"
    rows_off = docking.run_docking("1ABC", ligands, off, enumerate_stereo=False)
    check("enumeration off: one Vina call per valid ligand", len(vina_calls) == 3,
          str(len(vina_calls)))
    check("enumeration off: pose written directly, n_isomers=1",
          (off / "poses" / "two_centres.pdbqt").exists()
          and all(r["n_isomers"] == 1 for r in rows_off))

print()
print("ALL PASSED" if not failures else f"{len(failures)} FAILED: {failures}")
raise SystemExit(1 if failures else 0)
