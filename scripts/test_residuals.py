"""Checks for ligentry.pkpd.residuals (no Qt, numpy + scipy only).

    python scripts/test_residuals.py
"""
from __future__ import annotations

import sys
from dataclasses import replace

import numpy as np

from ligentry.pkpd.model import simulate
from ligentry.pkpd.observations import Observations
from ligentry.pkpd.residuals import residuals_for
from ligentry.pkpd.workflow import ModelInputs, build_model

_results: list[bool] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    _results.append(bool(ok))
    print(f"{'ok  ' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""))


def main() -> int:
    inputs = ModelInputs(ka=1.2, volume=40.0, clearance=6.0)
    pk, regimen, _ = build_model(inputs)
    times = np.array([0.5, 1, 2, 4, 8, 12, 24], dtype=float)
    conc = simulate(pk, regimen, times).conc
    obs = Observations(times=times, conc=conc)

    res = residuals_for(inputs, obs)
    check("one residual per observation", len(res) == times.size)
    check("noise-free data from the same model: residuals are zero", np.allclose(res.residual, 0.0, atol=1e-12))
    check("noise-free data: rmse is zero", res.rmse < 1e-12)
    check("predicted matches a direct simulation", np.allclose(res.predicted, conc, rtol=1e-12))

    bumped = Observations(times=times, conc=conc + 0.1)
    res = residuals_for(inputs, bumped)
    check("residual is observed minus predicted", np.allclose(res.residual, 0.1))
    check("rmse of a constant offset equals the offset", abs(res.rmse - 0.1) < 1e-9)

    low = Observations(times=times, conc=conc * 0.5)
    res = residuals_for(inputs, low)
    check("data below the model gives negative residuals", bool(np.all(res.residual < 0)))

    wrong = replace(inputs, clearance=3.0)
    check("a model that does not match gives a larger rmse",
          residuals_for(wrong, obs).rmse > residuals_for(inputs, obs).rmse + 1e-3)

    repeated = replace(inputs, n_doses=3, interval=12.0)
    check("the dosing in the form is used (single-dose data vs 3 doses)",
          residuals_for(repeated, obs).rmse > 0.05)

    unsorted = Observations(times=times[::-1].copy(), conc=conc[::-1].copy())
    res = residuals_for(inputs, unsorted)
    check("unsorted observation times are handled", np.allclose(res.residual, 0.0, atol=1e-12))

    res = residuals_for(inputs, obs)
    check("inputs are copies, not views of the data", res.observed is not obs.conc and res.times is not obs.times)

    passed = sum(_results)
    print(f"\n{passed}/{len(_results)} checks passed")
    return 0 if all(_results) else 1


if __name__ == "__main__":
    sys.exit(main())
