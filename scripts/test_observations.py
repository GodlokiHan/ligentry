"""Checks for biosim.pkpd.observations (no Qt, numpy + scipy only).

    python scripts/test_observations.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from biosim.pkpd.fit import fit_concentration
from biosim.pkpd.model import PKParameters, Regimen, simulate
from biosim.pkpd.observations import (
    ObservationError,
    load_observations,
    parse_observations,
)

_results: list[bool] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    _results.append(bool(ok))
    print(f"{'ok  ' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""))


def close(a, b) -> bool:
    return np.allclose(np.asarray(a, dtype=float), np.asarray(b, dtype=float), rtol=0, atol=1e-12)


def raises(text: str, fragment: str) -> tuple[bool, str]:
    try:
        parse_observations(text)
    except ObservationError as exc:
        return fragment.lower() in str(exc).lower(), str(exc)
    return False, "no error raised"


def main() -> int:
    # ---- formats -----------------------------------------------------------
    obs = parse_observations("time,conc\n0.5,1.8\n1,3.1\n2,2.4\n")
    check("header csv: times", close(obs.times, [0.5, 1, 2]))
    check("header csv: conc", close(obs.conc, [1.8, 3.1, 2.4]))
    check("len()", len(obs) == 3)

    obs = parse_observations("0.5,1.8\n1,3.1\n")
    check("no header: first column is time", close(obs.times, [0.5, 1]) and close(obs.conc, [1.8, 3.1]))

    obs = parse_observations("time;conc\n0,5;1,5\n1;2,25\n")
    check("semicolon + decimal comma", close(obs.times, [0.5, 1]) and close(obs.conc, [1.5, 2.25]))

    obs = parse_observations("Time (h)\tConc (mg/L)\n0.5\t1.8\n1\t3.1\n")
    check("tab separated with units in header", close(obs.times, [0.5, 1]) and close(obs.conc, [1.8, 3.1]))

    obs = parse_observations("0.5 1.8\n1   3.1\n")
    check("whitespace separated", close(obs.times, [0.5, 1]) and close(obs.conc, [1.8, 3.1]))

    obs = parse_observations("# study 12\n\ntime,conc\n\n1,2\n# note\n2,1\n")
    check("comments and blank lines skipped", close(obs.times, [1, 2]) and close(obs.conc, [2, 1]))

    obs = parse_observations("2,1\n0.5,1.8\n1,3.1\n")
    check("rows sorted by time", close(obs.times, [0.5, 1, 2]) and close(obs.conc, [1.8, 3.1, 1]))

    obs = parse_observations("subject,conc,time\nA,1.8,0.5\nA,3.1,1\n")
    check("named columns in any order, extra columns ignored",
          close(obs.times, [0.5, 1]) and close(obs.conc, [1.8, 3.1]))

    obs = parse_observations("hours,cp\n1,2\n2,1\n")
    check("header aliases (hours, cp)", close(obs.times, [1, 2]) and close(obs.conc, [2, 1]))

    obs = parse_observations("x,y\n1,2\n2,1\n")
    check("two unrecognised header names -> time, conc", close(obs.times, [1, 2]) and close(obs.conc, [2, 1]))

    obs = parse_observations("1,2\n1,2.2\n2,1\n")
    check("duplicate times kept", len(obs) == 3)

    obs = parse_observations("0,0\n1,3\n")
    check("zero time and zero concentration allowed", close(obs.conc, [0, 3]))

    # ---- errors ------------------------------------------------------------
    ok, msg = raises("", "empty")
    check("empty input", ok, msg)
    ok, msg = raises("# only a comment\n\n", "empty")
    check("comment-only input", ok, msg)
    ok, msg = raises("time,conc\n", "header")
    check("header only", ok, msg)
    ok, msg = raises("time,conc\n1,2\n-1,3\n", "line 3")
    check("negative time reports its line", ok, msg)
    ok, msg = raises("time,conc\n1,-2\n", "negative")
    check("negative concentration", ok, msg)
    ok, msg = raises("time,conc\n1,abc\n", "line 2")
    check("non-numeric concentration reports its line", ok, msg)
    ok, msg = raises("time,conc\nabc,1\n", "time")
    check("non-numeric time", ok, msg)
    ok, msg = raises("time,conc\n1,nan\n", "finite")
    check("nan rejected", ok, msg)
    ok, msg = raises("time,conc\n1,inf\n", "finite")
    check("inf rejected", ok, msg)
    ok, msg = raises("a,b,c\n1,2,3\n", "which columns")
    check("three unrecognised header names", ok, msg)
    ok, msg = raises("time,conc\n1\n", "columns")
    check("missing column in a row", ok, msg)

    # ---- files ---------------------------------------------------------------
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "data.csv"
        path.write_bytes("\ufefftime,conc\n1,2\n2,1\n".encode("utf-8"))  # Excel-style BOM
        obs = load_observations(path)
        check("file with BOM", close(obs.times, [1, 2]) and obs.source == str(path))
        try:
            load_observations(Path(tmp) / "missing.csv")
            check("missing file raises ObservationError", False)
        except ObservationError:
            check("missing file raises ObservationError", True)

        # ---- round trip: simulate -> write CSV -> load -> fit ------------------
        truth = PKParameters(ka=1.2, volume=40.0, clearance=6.0)
        regimen = Regimen(dose=100.0)
        times = np.array([0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24])
        conc = simulate(truth, regimen, times).conc
        rng = np.random.default_rng(7)
        noisy = conc * (1 + 0.03 * rng.standard_normal(conc.size))
        csv_path = Path(tmp) / "pk.csv"
        csv_path.write_text(
            "time,conc\n" + "\n".join(f"{t:g},{c:.6g}" for t, c in zip(times[::-1], noisy[::-1]))
        )
        obs = load_observations(csv_path)
        check("round trip: sorted after load", bool(np.all(np.diff(obs.times) > 0)))
        start = PKParameters(ka=0.6, volume=25.0, clearance=3.0)
        fit = fit_concentration(obs.times, obs.conc, regimen, start)
        check("round trip: fit succeeds", fit.success, fit.message)
        for name, true_value in (("ka", 1.2), ("volume", 40.0), ("clearance", 6.0)):
            est = fit.estimates[name]
            check(f"round trip: {name} within 15% ({est:.3g} vs {true_value})",
                  abs(est - true_value) / true_value < 0.15)

    passed = sum(_results)
    print(f"\n{passed}/{len(_results)} checks passed")
    return 0 if all(_results) else 1


if __name__ == "__main__":
    sys.exit(main())
