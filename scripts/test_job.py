"""Checks for ligentry.gui.job (the separate-process docking job). Run: python scripts/test_job.py"""
from __future__ import annotations

import contextlib
import io
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from ligentry.gui import job

failures: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(("PASS  " if ok else "FAIL  ") + label + (f"   {detail}" if detail else ""))
    if not ok:
        failures.append(label)


def run_main(spec_path: Path, **kw) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
        code = job.main([str(spec_path)], **kw)
    return code, buf.getvalue()


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    work = root / "run"
    work.mkdir()

    def write_spec(**over) -> Path:
        spec = {
            "pdb_id": "1HSG", "chembl_target": "", "smiles_text": "CCO ethanol\nCCN",
            "n_ligands": 5, "exhaustiveness": 12, "workdir": str(work),
            "enumerate_stereo": False, "max_isomers": 4,
        }
        spec.update(over)
        p = root / "job.json"
        p.write_text(json.dumps(spec))
        return p

    # SMILES job: options forwarded, STOP file wired to should_stop
    seen: dict = {}

    def fake_run(pdb_id, ligands, workdir, **kw):
        seen.update(pdb_id=pdb_id, ligands=ligands, workdir=workdir, **kw)
        seen["stop_before"] = kw["should_stop"]()
        (Path(workdir) / "STOP").write_text("x")
        seen["stop_after"] = kw["should_stop"]()
        return []

    code, out = run_main(write_spec(), run=fake_run)
    check("success returns 0", code == 0, out)
    check("ligands parsed (names and default names)",
          seen["ligands"] == [("ethanol", "CCO"), ("lig2", "CCN")], str(seen.get("ligands")))
    check("options forwarded", seen["exhaustiveness"] == 12
          and seen["enumerate_stereo"] is False and seen["max_isomers"] == 4)
    check("should_stop follows the STOP file",
          seen["stop_before"] is False and seen["stop_after"] is True)

    # ChEMBL job
    got: dict = {}

    def fake_fetch(target, limit):
        got.update(target=target, limit=limit)
        return [("CHEMBL1", "CCO")]

    (work / "STOP").unlink()
    code, out = run_main(write_spec(chembl_target="CHEMBL243", n_ligands=7),
                         run=lambda *a, **k: [], fetch=fake_fetch)
    check("ChEMBL target fetched with limit",
          code == 0 and got == {"target": "CHEMBL243", "limit": 7}, str(got))
    check("fetch announced on stdout", "Fetching actives for CHEMBL243" in out)

    # error paths
    code, out = run_main(write_spec(smiles_text="  \n# only a comment"), run=fake_run)
    check("no ligands -> exit 1 with ERROR line",
          code == 1 and "ERROR: No ligands to dock." in out, out)

    def boom(*a, **k):
        raise RuntimeError("boom")

    code, out = run_main(write_spec(), run=boom)
    check("exception -> exit 1 with ERROR line", code == 1 and "ERROR: boom" in out, out)

    bad = root / "bad.json"
    bad.write_text("not json")
    code, out = run_main(bad)
    check("bad job file -> exit 2", code == 2 and out.startswith("ERROR: bad job file"), out)
    code, out = run_main(root / "missing.json")
    check("missing job file -> exit 2", code == 2)

    # real subprocess, exactly how the GUI starts it: imports cleanly, fails gracefully
    p = write_spec(pdb_id="ZZZZ")
    proc = subprocess.run(
        [sys.executable, "-u", "-m", "ligentry.gui.job", str(p)],
        capture_output=True, text=True, timeout=120,
    )
    check("subprocess: bad PDB id -> exit 1 and an ERROR line (module imports outside Qt)",
          proc.returncode == 1 and "ERROR:" in proc.stdout, (proc.stdout + proc.stderr)[-300:])

print()
print("ALL PASSED" if not failures else f"{len(failures)} FAILED: {failures}")
raise SystemExit(1 if failures else 0)
