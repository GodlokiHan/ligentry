"""Checks for effect-data loading, effect fitting and its glue (biosim.pkpd.effect_*).

    python scripts/test_effect.py

Needs only numpy and scipy. No network, no RDKit, no Qt.
"""
from __future__ import annotations

import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np

from biosim.pkpd.effect_data import (
    load_effect_observations,
    parse_effect_observations,
)
from biosim.pkpd.effect_fit import fit_effect, guess_effect_start
from biosim.pkpd.effect_workflow import (
    apply_effect_fit,
    concentrations_at,
    effect_fit_notes,
    effect_fit_table_rows,
    predict_effects,
    run_effect_fit,
)
from biosim.pkpd.model import EffectParameters, effect_from_conc
from biosim.pkpd.observations import ObservationError
from biosim.pkpd.workflow import ModelInputs

PASSED = 0
FAILED = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASSED, FAILED
    if condition:
        PASSED += 1
        print(f"ok   - {name}")
    else:
        FAILED += 1
        print(f"FAIL - {name} {detail}")


def raises(fn, exc=ValueError) -> bool:
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


def rel_errors(result, truth: EffectParameters) -> dict[str, float]:
    return {
        n: abs(v - getattr(truth, n)) / max(abs(getattr(truth, n)), 1e-12)
        for n, v in result.estimates.items()
    }


# ---------------------------------------------------------------- loader
obs = parse_effect_observations("time,effect\n0,2.0\n1,15.5\n2,30\n")
check("loads a simple file", len(obs) == 3 and obs.effect[1] == 15.5 and obs.times[2] == 2.0)

obs = parse_effect_observations("response;hours\n1,5;3\n0,2;1\n")
check("finds columns by name in any order, decimal comma", obs.times.tolist() == [1.0, 3.0] and obs.effect.tolist() == [0.2, 1.5])
check("rows are sorted by time", np.all(np.diff(obs.times) >= 0))

obs = parse_effect_observations("time,conc,effect\n0,0,1\n1,5,9\n")
check("picks the effect column when a conc column is present", obs.effect.tolist() == [1.0, 9.0])

obs = parse_effect_observations("# comment\n\n0 1.5\n2 3.5\n")
check("headerless whitespace file reads time then effect", obs.effect.tolist() == [1.5, 3.5])

obs = parse_effect_observations("time,effect\n0,-4.5\n1,-20\n")
check("negative effects are allowed", obs.effect.tolist() == [-4.5, -20.0])

obs = parse_effect_observations("a,b\n0,1\n1,2\n")
check("two unrecognised header names are read as time, effect", obs.effect.tolist() == [1.0, 2.0])

obs = parse_effect_observations("\ufefftime\teffect\n0\t1,5\n")
check("tab separated with decimal comma", obs.effect.tolist() == [1.5])

check("rejects empty input", raises(lambda: parse_effect_observations(""), ObservationError))
check("rejects a header-only file", raises(lambda: parse_effect_observations("time,effect\n"), ObservationError))
check("rejects a headerless file with three columns", raises(lambda: parse_effect_observations("0,1,2\n1,2,3\n"), ObservationError))
check("rejects negative times", raises(lambda: parse_effect_observations("time,effect\n-1,3\n"), ObservationError))
check("rejects non-numeric effect", raises(lambda: parse_effect_observations("time,effect\n1,abc\n"), ObservationError))
check("rejects non-finite values", raises(lambda: parse_effect_observations("time,effect\n1,inf\n"), ObservationError))
check("rejects an ambiguous 3-column header", raises(lambda: parse_effect_observations("x,y,z\n1,2,3\n"), ObservationError))
try:
    parse_effect_observations("time,effect\n0,1\n1,oops\n")
except ObservationError as exc:
    check("errors name the line number", "Line 3" in str(exc))

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "pd.csv"
    path.write_text("time,effect\n0,1\n1,2\n")
    check("load_effect_observations reads a file", len(load_effect_observations(path)) == 2)
    check("load_effect_observations reports a missing file", raises(lambda: load_effect_observations(Path(tmp) / "nope.csv"), ObservationError))

# ---------------------------------------------------------------- start guess
conc = np.array([0.0, 0.2, 0.5, 1, 2, 4, 8, 16, 32.0])
truth = EffectParameters(e0=5.0, emax=40.0, ec50=3.0, hill=1.0)
y = effect_from_conc(conc, truth)
g = guess_effect_start(conc, y)
check("start guess puts the baseline near the lowest-concentration effects", abs(g.e0 - 5.0) < 3.0)
check("start guess has the right sign and order of magnitude for Emax", 0.3 * 40 < g.emax < 3 * 40)
check("start guess for EC50 is inside the concentration range", 0.2 < g.ec50 < 32)
g_inh = guess_effect_start(conc, effect_from_conc(conc, EffectParameters(e0=100.0, emax=-80.0, ec50=3.0)))
check("start guess recognises an inhibitory effect", g_inh.emax < 0)

# ---------------------------------------------------------------- fitting
fit = fit_effect(conc, y)
check("fit converges on noise-free data", fit.success)
check("fit recovers E0, Emax and EC50 to 0.1%", max(rel_errors(fit, truth).values()) < 1e-3, f"({rel_errors(fit, truth)})")
check("fixed Hill slope is left at 1", fit.effect.hill == 1.0)
check("fit reports uncertainty with spare points", fit.cv_percent is not None and set(fit.cv_percent) == set(fit.free))
check("fit reports standard errors", fit.se is not None and all(v >= 0 for v in fit.se.values()))
check("noise-free fit has a tiny RMSE", fit.rmse < 1e-6)

truth_h = EffectParameters(e0=5.0, emax=40.0, ec50=3.0, hill=2.2)
c_h = np.geomspace(0.2, 60, 14)
fit = fit_effect(c_h, effect_from_conc(c_h, truth_h), free=("e0", "emax", "ec50", "hill"))
check("fit recovers a Hill slope", max(rel_errors(fit, truth_h).values()) < 1e-2, f"({rel_errors(fit, truth_h)})")

truth_i = EffectParameters(e0=100.0, emax=-80.0, ec50=3.0, hill=1.0)
fit = fit_effect(conc, effect_from_conc(conc, truth_i))
check("fit recovers an inhibitory effect", max(rel_errors(fit, truth_i).values()) < 1e-3, f"({rel_errors(fit, truth_i)})")

rng = np.random.default_rng(3)
noisy = y + 1.0 * rng.standard_normal(y.size)
fit = fit_effect(conc, noisy)
errs = rel_errors(fit, truth)
check("noisy fit stays within 25% of the true values", max(errs.values()) < 0.25, f"({errs})")
check("noisy fit has an RMSE of about the noise level", 0.3 < fit.rmse < 2.0, f"({fit.rmse:.3f})")

fit = fit_effect(conc, y, start=EffectParameters(e0=0.0, emax=100.0, ec50=1.0))
check("fit also converges from a poor explicit start", max(rel_errors(fit, truth).values()) < 1e-3)

part = fit_effect(conc, y, start=replace(truth, e0=5.0), free=("emax", "ec50"))
check("fixed E0 is left unchanged", part.effect.e0 == 5.0)

spread = [0, 3, 6]
three = fit_effect(conc[spread], y[spread], free=("e0", "emax", "ec50"))
check("no uncertainty when observations equal parameters", three.cv_percent is None and three.se is None)

tight = fit_effect(c_h, effect_from_conc(c_h, truth_h), free=("e0", "emax", "ec50", "hill"), bounds={"hill": (0.1, 1.0)})
check("a limit that binds is reported in at_bound", "hill" in tight.at_bound)
check("the fitted Hill slope respects its limit", tight.effect.hill <= 1.0 + 1e-6)
check("an unconstrained fit reports nothing at a limit", fit_effect(conc, y).at_bound == ())

# data confined to the linear part of the curve: EC50 and Emax are poorly determined
low = np.linspace(0.05, 1.0, 12)
far = EffectParameters(e0=0.0, emax=100.0, ec50=50.0)
rng = np.random.default_rng(5)
flat = effect_from_conc(low, far) + 0.2 * rng.standard_normal(low.size)
fit = fit_effect(low, flat)
check("unsaturated data give a large CV for EC50", fit.cv_percent["ec50"] is None or fit.cv_percent["ec50"] > 50, f"({fit.cv_percent})")

# ---------------------------------------------------------------- input validation
check("rejects mismatched lengths", raises(lambda: fit_effect(conc, y[:-1])))
check("rejects negative concentrations", raises(lambda: fit_effect(-conc, y)))
check("rejects NaN effects", raises(lambda: fit_effect(conc, np.where(conc > 1, np.nan, y))))
check("rejects an unknown parameter", raises(lambda: fit_effect(conc, y, free=("foo",))))
check("rejects an empty free list", raises(lambda: fit_effect(conc, y, free=())))
check("rejects a repeated free parameter", raises(lambda: fit_effect(conc, y, free=("e0", "e0"))))
check("rejects too few observations", raises(lambda: fit_effect(conc[:2], y[:2])))
check("rejects all-zero concentrations", raises(lambda: fit_effect(np.zeros(5), np.arange(5.0))))
check("rejects reversed bounds", raises(lambda: fit_effect(conc, y, bounds={"ec50": (5.0, 1.0)})))
check(
    "rejects non-positive log-scale bounds",
    raises(lambda: fit_effect(conc, y, free=("e0", "emax", "ec50", "hill"), bounds={"hill": (0.0, 5.0)})),
)

# ---------------------------------------------------------------- glue with the PK model
inputs = ModelInputs(dose=200.0, n_doses=1, ka=1.2, volume=30.0, clearance=5.0)
t_pd = np.array([0.0, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24], dtype=float)
pd_truth = EffectParameters(e0=10.0, emax=60.0, ec50=2.0, hill=1.0)
c_at = concentrations_at(inputs, parse_effect_observations("time,effect\n" + "\n".join(f"{t},0" for t in t_pd)))
effect_obs = parse_effect_observations(
    "time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(t_pd, effect_from_conc(c_at, pd_truth)))
)
check("concentrations come from the PK model", c_at[0] == 0.0 and c_at.max() > 1.0)

result, conc_used = run_effect_fit(inputs, effect_obs)
check("glue fit recovers E0, Emax and EC50", max(rel_errors(result, pd_truth).values()) < 1e-3, f"({rel_errors(result, pd_truth)})")
check("glue returns the concentrations it used", np.allclose(conc_used, c_at))

keep_hill = replace(inputs, hill=1.0)
result2, _ = run_effect_fit(replace(keep_hill, e0=0.0), effect_obs, free=("e0", "emax", "ec50"))
check("fixed parameters keep the form values", result2.effect.hill == 1.0)

result3, _ = run_effect_fit(replace(inputs, e0=9.0, emax=55.0, ec50=2.5), effect_obs, use_form_start=True)
check("fit can start from the form values", result3.start.ec50 == 2.5 and result3.start.e0 == 9.0)

applied = apply_effect_fit(inputs, result)
check("apply copies the fit into the form and switches the effect on", applied.use_effect and abs(applied.ec50 - 2.0) < 1e-2)
check("apply leaves the PK values alone", applied.ka == inputs.ka and applied.clearance == inputs.clearance)

conc_p, pred = predict_effects(applied, effect_obs)
check("predictions match the data after applying the fit", np.max(np.abs(pred - effect_obs.effect)) < 1e-3)

rows = effect_fit_table_rows(result)
check("table rows have parameter, estimate, unit and uncertainty", len(rows) == 3 and all(len(r) == 4 for r in rows))
check("table uncertainty is a standard error", all(r[3].startswith("\u00b1") for r in rows))
check("table row labels and units are right", rows[0][0] == "E0" and rows[2][2] == "mg/L")

notes = effect_fit_notes(result, conc_used)
check("notes always remind that the effect is direct", any("no delay" in n for n in notes))
check("a clean fit has no warning about convergence", not any("did not converge" in n for n in notes))

# only low concentrations measured -> EC50 above the data range
low_obs = parse_effect_observations("time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(t_pd[:4], effect_from_conc(c_at[:4], EffectParameters(e0=0.0, emax=500.0, ec50=80.0)))))
low_result, low_conc = run_effect_fit(inputs, low_obs, free=("emax", "ec50"))
low_notes = effect_fit_notes(low_result, low_conc)
check("notes warn when EC50 is above every measured concentration", any("never approaches saturation" in n for n in low_notes), f"({low_notes})")

sat = EffectParameters(e0=0.0, emax=50.0, ec50=0.001)
sat_obs = parse_effect_observations("time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(t_pd, effect_from_conc(c_at, sat))))
sat_result, sat_conc = run_effect_fit(inputs, sat_obs, free=("e0", "emax", "ec50"))
sat_notes = effect_fit_notes(sat_result, sat_conc)
check("notes warn when EC50 is below every non-zero concentration", any("saturated everywhere" in n for n in sat_notes), f"({sat_notes})")
check("a baseline of zero is not flagged as poorly determined", not any("Poorly determined" in n for n in sat_notes), f"({sat_notes})")
check("table shows a plain standard error for a baseline near zero", effect_fit_table_rows(sat_result)[0][3].startswith("\u00b1"))

rng = np.random.default_rng(11)
very_noisy = parse_effect_observations(
    "time,effect\n" + "\n".join(
        f"{t},{e:.10g}" for t, e in zip(t_pd, effect_from_conc(c_at, pd_truth) + 40.0 * rng.standard_normal(t_pd.size))
    )
)
vn_result, vn_conc = run_effect_fit(inputs, very_noisy)
check("very noisy data are flagged as poorly determined", any("Poorly determined" in n for n in effect_fit_notes(vn_result, vn_conc)), f"({effect_fit_notes(vn_result, vn_conc)})")

before_drug = parse_effect_observations("time,effect\n0,1\n0,1.1\n0,0.9\n")
check("fit refuses data taken before any drug is present", raises(lambda: run_effect_fit(inputs, before_drug, free=("e0", "emax", "ec50"))))

print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
