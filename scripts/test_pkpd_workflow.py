"""Checks for ligentry.pkpd.workflow (no Qt, numpy + scipy only).

    python scripts/test_pkpd_workflow.py
"""
from __future__ import annotations

import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np

from ligentry.pkpd.model import simulate
from ligentry.pkpd.observations import Observations, load_observations
from ligentry.pkpd.workflow import (
    LINEAR,
    MICHAELIS_MENTEN,
    ModelInputs,
    apply_fit,
    build_model,
    default_time_span,
    fit_notes,
    fit_table_rows,
    free_choices,
    run_fit,
    simulate_inputs,
    simulation_to_csv,
)

_results: list[bool] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    _results.append(bool(ok))
    print(f"{'ok  ' if ok else 'FAIL'} {name}" + (f"  [{detail}]" if detail and not ok else ""))


def raises_value_error(fn, fragment: str = "") -> bool:
    try:
        fn()
    except ValueError as exc:
        return fragment.lower() in str(exc).lower()
    return False


def main() -> int:
    base = ModelInputs()  # 100 mg, ka 1, V 40, CL 5, single dose

    # ---- build_model -------------------------------------------------------------
    pk, regimen, effect = build_model(base)
    check("linear model built from defaults",
          pk.elimination == LINEAR and pk.clearance == 5.0 and regimen.n_doses == 1)
    check("no effect unless requested", effect is None)

    pk, _, effect = build_model(replace(base, use_effect=True, e0=2.0, emax=-1.0, ec50=3.0, hill=2.0))
    check("effect parameters passed through",
          effect is not None and effect.e0 == 2.0 and effect.emax == -1.0 and effect.hill == 2.0)

    pk, _, _ = build_model(replace(base, elimination=MICHAELIS_MENTEN, vmax=10.0, km=1.5))
    check("Michaelis-Menten model built", pk.elimination == MICHAELIS_MENTEN and pk.vmax == 10.0)

    check("bad dose rejected", raises_value_error(lambda: build_model(replace(base, dose=0)), "dose"))
    check("bad ka rejected", raises_value_error(lambda: build_model(replace(base, ka=-1)), "ka"))
    check("multiple doses need an interval",
          raises_value_error(lambda: build_model(replace(base, n_doses=3, interval=0)), "interval"))
    check("unknown elimination rejected",
          raises_value_error(lambda: build_model(replace(base, elimination="x")), "unknown"))

    # ---- free_choices ------------------------------------------------------------
    check("fit choices, linear", free_choices(LINEAR) == ("ka", "volume", "clearance"))
    check("fit choices, Michaelis-Menten",
          free_choices(MICHAELIS_MENTEN) == ("ka", "volume", "vmax", "km"))
    check("fit choices, unknown rejected", raises_value_error(lambda: free_choices("x"), "unknown"))

    # ---- time span ---------------------------------------------------------------
    pk, regimen, _ = build_model(base)
    span = default_time_span(pk, regimen)
    check("default span covers 5 half-lives", span >= 5 * np.log(2) / pk.ke)
    check("default span is a tidy number", span == 36.0, f"got {span}")
    pk_slow = build_model(replace(base, clearance=0.5))[0]
    check("slower elimination gives a longer span",
          default_time_span(pk_slow, regimen) > span)
    pk_mm, reg_mm, _ = build_model(replace(base, elimination=MICHAELIS_MENTEN))
    check("Michaelis-Menten span is positive and finite",
          0 < default_time_span(pk_mm, reg_mm) <= 100_000)
    _, reg_multi, _ = build_model(replace(base, n_doses=4, interval=12.0))
    check("span extends past the last dose", default_time_span(pk, reg_multi) > 36.0)

    # ---- simulate_inputs -----------------------------------------------------------
    res = simulate_inputs(replace(base, t_end=200.0))
    check("AUC equals F*dose/CL over a long window",
          abs(res.metrics["auc"] - 100.0 / 5.0) / 20.0 < 0.005, f"auc={res.metrics['auc']}")
    check("t_end is honoured", res.t_end == 200.0 and abs(res.sim.t[-1] - 200.0) < 1e-12)
    check("no effect curve without effect", res.sim.effect is None)

    res = simulate_inputs(replace(base, n_doses=3, interval=12.0, t_end=60.0))
    dose_times = [0.0, 12.0, 24.0]
    check("dose times are on the time grid", all(np.any(np.isclose(res.sim.t, d)) for d in dose_times))
    check("time grid is sorted and unique", bool(np.all(np.diff(res.sim.t) > 0)))

    res = simulate_inputs(replace(base, use_effect=True, t_end=48.0))
    check("effect curve present with effect", res.sim.effect is not None and res.sim.effect.size == res.sim.t.size)
    check("effect starts at E0", abs(res.sim.effect[0] - 0.0) < 1e-12)

    res = simulate_inputs(base)
    check("automatic window is used when t_end is None", res.t_end == 36.0)
    res = simulate_inputs(base, min_t_end=100.0)
    check("automatic window stretches to cover data", res.t_end >= 100.0)
    res = simulate_inputs(base, min_t_end=10.0)
    check("data inside the window does not shrink it", res.t_end == 36.0)

    check("end time of zero rejected",
          raises_value_error(lambda: simulate_inputs(replace(base, t_end=0.0)), "end time"))
    check("huge end time rejected",
          raises_value_error(lambda: simulate_inputs(replace(base, t_end=1e9)), "end time"))

    # ---- fitting -------------------------------------------------------------------
    truth = replace(base, ka=1.2, volume=40.0, clearance=6.0)
    pk_t, reg_t, _ = build_model(truth)
    times = np.array([0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24], dtype=float)
    conc = simulate(pk_t, reg_t, times).conc
    obs = Observations(times=times, conc=conc, source="synthetic")
    start = replace(base, ka=0.6, volume=25.0, clearance=3.0)
    result = run_fit(start, obs, ("ka", "volume", "clearance"))
    check("fit succeeds on noise-free data", result.success, result.message)
    check("fit recovers ka, V, CL to 1%",
          all(abs(result.estimates[n] - v) / v < 0.01
              for n, v in (("ka", 1.2), ("volume", 40.0), ("clearance", 6.0))))

    applied = apply_fit(start, result)
    check("apply_fit copies the estimates",
          abs(applied.ka - result.estimates["ka"]) < 1e-12
          and abs(applied.clearance - result.estimates["clearance"]) < 1e-12)
    check("apply_fit leaves other values alone",
          applied.dose == start.dose and applied.interval == start.interval)

    rows = fit_table_rows(result)
    check("table has one row per fitted parameter", [r[0] for r in rows] == ["ka", "V", "CL"])
    check("table units", [r[2] for r in rows] == ["1/h", "L", "L/h"])
    check("table estimates are formatted text", all(isinstance(r[1], str) and r[1] for r in rows))
    check("table shows CV as text", all(isinstance(r[3], str) for r in rows))

    check("fit with a subset of parameters", run_fit(start, obs, ("clearance",)).free == ("clearance",))
    check("fit rejects an inapplicable parameter",
          raises_value_error(lambda: run_fit(start, obs, ("vmax",)), "no starting value"))

    # too few observations for an uncertainty estimate -> notes say so
    small = Observations(times=times[:3], conc=conc[:3])
    small_fit = run_fit(start, small, ("ka", "volume", "clearance"))
    check("no CV with as many points as parameters", small_fit.cv_percent is None)
    check("notes mention missing uncertainty",
          any("too few observations" in n.lower() for n in fit_notes(small_fit)))
    check("table shows n/a without CV", all(r[3] == "n/a" for r in fit_table_rows(small_fit)))
    check("apparent V and CL are flagged",
          any("apparent" in n.lower() for n in fit_notes(result)))
    check("notes are plain strings", all(isinstance(n, str) for n in fit_notes(result)))

    flip = replace(result, params=replace(result.params, ka=0.05))
    check("flip-flop warning when ka is below ke",
          any("flip-flop" in n.lower() for n in fit_notes(flip)))
    check("no flip-flop warning for a sensible fit",
          not any("flip-flop" in n.lower() for n in fit_notes(result)))

    # ---- CSV export round trip -------------------------------------------------------
    with tempfile.TemporaryDirectory() as tmp:
        res = simulate_inputs(replace(base, use_effect=True, t_end=24.0))
        path = simulation_to_csv(Path(tmp) / "sim.csv", res)
        header = path.read_text().splitlines()[0]
        check("CSV header has the effect column", header == "time_h,conc_mg_per_L,gut_mg,central_mg,effect")
        loaded = load_observations(path)
        check("exported CSV loads as observations", len(loaded) == res.sim.t.size)
        check("exported concentrations survive the round trip",
              np.allclose(loaded.conc, res.sim.conc, rtol=1e-6, atol=1e-9))
        res2 = simulate_inputs(replace(base, t_end=24.0))
        path2 = simulation_to_csv(Path(tmp) / "sim2.csv", res2)
        check("CSV without effect has no effect column",
              path2.read_text().splitlines()[0] == "time_h,conc_mg_per_L,gut_mg,central_mg")

    passed = sum(_results)
    print(f"\n{passed}/{len(_results)} checks passed")
    return 0 if all(_results) else 1


if __name__ == "__main__":
    sys.exit(main())
