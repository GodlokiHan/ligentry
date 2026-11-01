"""Tests for biosim.pkpd.effect_residuals (no Qt). Run: python scripts/test_effect_residuals.py"""
from __future__ import annotations

import numpy as np

from biosim.pkpd.effect_data import EffectObservations
from biosim.pkpd.effect_residuals import effect_residuals_for
from biosim.pkpd.effect_workflow import concentrations_at, form_effect
from biosim.pkpd.model import effect_from_conc
from biosim.pkpd.workflow import ModelInputs

failures = 0


def check(condition: bool, label: str) -> None:
    global failures
    print(("PASS " if condition else "FAIL ") + label)
    if not condition:
        failures += 1


inputs = ModelInputs(e0=5.0, emax=80.0, ec50=0.8, hill=1.0)  # use_effect stays False
times = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 24.0])
blank = EffectObservations(times=times, effect=np.zeros_like(times))
conc = concentrations_at(inputs, blank)
truth = effect_from_conc(conc, form_effect(inputs))

# data exactly on the model -> zero residuals
res = effect_residuals_for(inputs, EffectObservations(times=times, effect=truth.copy()))
check(len(res) == times.size, "one residual per observation")
check(np.allclose(res.times, times), "times kept")
check(np.allclose(res.predicted, truth), "predicted effect comes from the form's model")
check(np.allclose(res.residual, 0.0, atol=1e-9), "zero residuals for data on the model")
check(res.rmse < 1e-9, "rmse is zero for data on the model")

# works whether or not the effect box is ticked
on = ModelInputs(**{**inputs.__dict__, "use_effect": True})
check(np.allclose(effect_residuals_for(on, EffectObservations(times, truth)).predicted, truth),
      "same prediction with the effect switched on")

# constant offset -> constant residual, observed minus predicted
res = effect_residuals_for(inputs, EffectObservations(times=times, effect=truth + 3.0))
check(np.allclose(res.residual, 3.0), "residual is observed minus predicted")
check(abs(res.rmse - 3.0) < 1e-9, "rmse equals the offset")

# a wrong baseline in the form shows up in the residuals
off = ModelInputs(**{**inputs.__dict__, "e0": 15.0})
res = effect_residuals_for(off, EffectObservations(times=times, effect=truth))
check(np.allclose(res.residual, -10.0), "baseline error of 10 gives residuals of -10")

# the returned arrays are copies, not views of the observations
obs = EffectObservations(times=times.copy(), effect=truth.copy())
res = effect_residuals_for(inputs, obs)
res.observed[0] = 999.0
res.times[0] = 999.0
check(obs.effect[0] == truth[0] and obs.times[0] == times[0], "observations are not modified")

# invalid form values raise ValueError instead of returning nonsense
bad = ModelInputs(**{**inputs.__dict__, "dose": -1.0})
try:
    effect_residuals_for(bad, blank)
    check(False, "a negative dose is rejected")
except ValueError:
    check(True, "a negative dose is rejected")

print("\nALL PASSED" if failures == 0 else f"\n{failures} FAILED")
raise SystemExit(1 if failures else 0)
