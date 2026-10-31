"""Checks for biosim.gui.runio (pure Python). Run: python scripts/test_runio.py"""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path

from biosim.gui import runio

failures: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(("PASS  " if ok else "FAIL  ") + label + (f"   {detail}" if detail else ""))
    if not ok:
        failures.append(label)


def raises(exc, fn, *a) -> bool:
    try:
        fn(*a)
    except exc:
        return True
    return False


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)

    # old-format run (before stereoisomer enumeration): no isomer columns
    old = root / "old_run"
    (old / "poses").mkdir(parents=True)
    (old / "1HSG_protein.pdb").write_text("END\n")
    (old / "poses" / "CHEMBL1.pdbqt").write_text("x\n")
    (old / "results.csv").write_text(
        "rank,name,smiles,vina_score_kcal_mol,status\n"
        "1,CHEMBL1,CCO,-7.5,ok\n"
        "2,CHEMBL2,CCN,-6.1,ok\n"
        "3,CHEMBL3,bad,,failed: RDKit could not parse SMILES\n"
    )
    pdb_id, rows = runio.load_run(old)
    check("pdb id taken from *_protein.pdb", pdb_id == "1HSG", pdb_id)
    check("3 rows loaded", len(rows) == 3)
    check("scores parsed, blank -> None",
          rows[0]["vina_score_kcal_mol"] == -7.5 and rows[2]["vina_score_kcal_mol"] is None)
    check("old CSV gets defaults", rows[0]["n_isomers"] == 1
          and rows[0]["best_isomer_smiles"] == "CCO" and rows[0]["notes"] == "")
    check("failed status kept with its comma-free text",
          rows[2]["status"].startswith("failed"))

    # notes + save round trip
    rows[1]["notes"] = 'promising, check "H-bond" later'
    runio.save_results(old, rows)
    pdb_id2, rows2 = runio.load_run(old)
    check("notes survive save/load (commas and quotes)",
          rows2[1]["notes"] == 'promising, check "H-bond" later')
    check("scores survive round trip", [r["vina_score_kcal_mol"] for r in rows2]
          == [-7.5, -6.1, None])
    check("no .tmp file left behind", not list(old.glob("*.tmp")))

    # removal then save renumbers ranks
    del rows2[0]
    runio.save_results(old, rows2)
    with open(old / "results.csv") as fh:
        saved = list(csv.DictReader(fh))
    check("ranks renumbered after removal", [r["rank"] for r in saved] == ["1", "2"])

    # export to another path
    out = runio.export_csv(root / "export.csv", rows2)
    check("export writes a CSV", out.exists() and "name" in out.read_text().splitlines()[0])

    # copy run
    dest = root / "copy_of_run"
    runio.copy_run(old, dest)
    check("copy includes receptor and poses",
          (dest / "1HSG_protein.pdb").exists() and (dest / "poses" / "CHEMBL1.pdbqt").exists())
    check("copy onto existing folder refused", raises(runio.RunFolderError,
          runio.copy_run, old, dest))
    check("copy into itself refused", raises(runio.RunFolderError,
          runio.copy_run, old, old / "inside"))

    # error handling
    empty = root / "empty"
    empty.mkdir()
    check("folder without results.csv rejected", raises(runio.RunFolderError,
          runio.load_run, empty))
    check("missing folder rejected", raises(runio.RunFolderError,
          runio.load_run, root / "nope"))
    (empty / "results.csv").write_text("foo,bar\n1,2\n")
    check("wrong CSV rejected", raises(runio.RunFolderError, runio.load_run, empty))
    check("save into missing folder refused", raises(runio.RunFolderError,
          runio.save_results, root / "nope", []))

    # run without receptor still loads
    (empty / "results.csv").write_text("name,smiles\nA,CC\n")
    pdb_id3, rows3 = runio.load_run(empty)
    check("no receptor -> empty pdb id, rows still load", pdb_id3 == "" and len(rows3) == 1)

print()
print("ALL PASSED" if not failures else f"{len(failures)} FAILED: {failures}")
raise SystemExit(1 if failures else 0)
