"""Checks for the effect-compartment (delayed effect) model, its fitting and workflow glue.

    python scripts/test_effect_compartment.py

Needs only numpy and scipy. No network, no RDKit, no Qt.
"""
from __future__ import annotations

import csv
import math
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np

from ligentry.pkpd import workflow as wf
from ligentry.pkpd.effect_data import parse_effect_observations
from ligentry.pkpd.effect_fit import DEFAULT_KE0_BOUNDS, fit_effect
from ligentry.pkpd.effect_residuals import effect_residuals_for
from ligentry.pkpd.effect_workflow import (
    apply_effect_fit,
    driving_concentrations_at,
    effect_fit_notes,
    effect_fit_table_rows,
    effect_free_choices,
    form_effect,
    predict_effects,
    run_effect_fit,
)
from ligentry.pkpd.model import (
    EffectParameters,
    PKParameters,
    Regimen,
    effect_from_conc,
    effect_site_concentration,
    equilibration_half_life,
    simulate,
)

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


def rel(a, b) -> float:
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))) / np.max(np.abs(b)))


def has(notes, text) -> bool:
    return any(text in n for n in notes)


# ------------------------------------------------------------------ closed forms
V, CL, D = 20.0, 4.0, 100.0
K = CL / V
ONE = PKParameters(ka=1.0, volume=V, clearance=CL)
IV = Regimen(dose=D, route="iv_bolus")
t = np.linspace(0, 40, 401)

for ke0 in (0.05, 0.7, 5.0):
    exact = D / V * ke0 / (ke0 - K) * (np.exp(-K * t) - np.exp(-ke0 * t))
    err = rel(effect_site_concentration(ONE, IV, t, ke0), exact)
    check(f"IV bolus effect-site concentration matches the closed form (ke0 {ke0})", err < 1e-6, f"(rel err {err:.2e})")

exact = D / V * K * t * np.exp(-K * t)
check("ke0 equal to the elimination rate matches the limit t*exp(-kt)",
      rel(effect_site_concentration(ONE, IV, t, K), exact) < 1e-6)

# ------------------------------------------------------------------ simulation behaviour
eff_direct = EffectParameters(e0=5.0, emax=50.0, ec50=2.0, hill=1.5)
eff_delay = replace(eff_direct, ke0=0.4)
plain = simulate(ONE, Regimen(dose=D), t)
direct = simulate(ONE, Regimen(dose=D), t, eff_direct)
delayed = simulate(ONE, Regimen(dose=D), t, eff_delay)

check("a direct effect has no effect-site concentration", direct.effect_conc is None)
check("a delayed effect carries the effect-site concentration", delayed.effect_conc is not None and delayed.effect_conc.shape == t.shape)
check("the effect compartment does not change the PK",
      rel(delayed.conc, plain.conc) < 1e-7 and rel(delayed.central, plain.central) < 1e-7)
check("direct effect is still computed from plasma concentration",
      np.allclose(direct.effect, effect_from_conc(direct.conc, eff_direct)))
check("delayed effect is computed from the effect-site concentration",
      np.allclose(delayed.effect, effect_from_conc(delayed.effect_conc, eff_delay)))
check("effect-site concentration starts at zero", delayed.effect_conc[0] == 0.0)
check("effect-site concentration is never negative", np.all(delayed.effect_conc >= 0.0))
check("the effect site peaks after plasma", t[np.argmax(delayed.effect_conc)] > t[np.argmax(delayed.conc)])
check("the effect site peaks lower than plasma", delayed.effect_conc.max() < delayed.conc.max())
check("the delayed effect peaks later than the direct effect", t[np.argmax(delayed.effect)] > t[np.argmax(direct.effect)])

t_long = np.linspace(0, 3000, 300001)
sim_long = simulate(ONE, Regimen(dose=D), t_long, eff_delay)
trap = getattr(np, "trapezoid", None) or np.trapz
auc_c = float(trap(sim_long.conc, t_long))
auc_ce = float(trap(sim_long.effect_conc, t_long))
check("effect-site AUC equals plasma AUC (the compartment has unit gain)", abs(auc_ce - auc_c) / auc_c < 1e-4, f"({auc_ce:.4f} vs {auc_c:.4f})")

fast = simulate(ONE, IV, np.array([1.0, 5.0, 20.0]), replace(eff_direct, ke0=1e4))
check("a very large ke0 gives back the direct effect",
      np.allclose(fast.effect, effect_from_conc(fast.conc, eff_direct), rtol=1e-3, atol=1e-3))

# rising and falling plasma concentration at the same level give different effects: hysteresis
sim_h = simulate(ONE, Regimen(dose=D), t, eff_delay)
i_peak = int(np.argmax(sim_h.conc))
target = 0.5 * sim_h.conc.max()
up = int(np.argmin(np.abs(sim_h.conc[:i_peak] - target)))
down = i_peak + int(np.argmin(np.abs(sim_h.conc[i_peak:] - target)))
check("the same plasma concentration gives different effects on the way up and down (hysteresis)",
      abs(sim_h.effect[up] - sim_h.effect[down]) > 1.0)
sim_d = simulate(ONE, Regimen(dose=D), t, eff_direct)
check("a direct effect shows far less hysteresis than a delayed one",
      abs(sim_d.effect[up] - sim_d.effect[down]) < 0.1 * abs(sim_h.effect[up] - sim_h.effect[down]))

# repeated IV doses: the effect site is a linear filter of plasma, so doses add
reg3 = Regimen(dose=60.0, interval=8.0, n_doses=3, route="iv_bolus")
tt = np.linspace(0, 40, 801)
total = np.zeros_like(tt)
for d in reg3.dose_times():
    single = effect_site_concentration(ONE, Regimen(dose=60.0, route="iv_bolus"), np.clip(tt - d, 0, None), 0.4)
    total += np.where(tt >= d, single, 0.0)
check("repeated doses obey superposition in the effect site",
      rel(effect_site_concentration(ONE, reg3, tt, 0.4), total) < 1e-6)

end = simulate(ONE, Regimen(dose=100.0, interval=8.0, n_doses=2, route="iv_bolus"), np.array([7.999999, 8.0]), eff_delay)
check("a dose at the end of the window still works with an effect compartment", end.conc[1] > end.conc[0] and end.effect_conc is not None)

shuffled = np.array([20.0, 0.0, 3.0, 12.0, 1.0])
a = simulate(ONE, IV, shuffled, eff_delay).effect_conc
b = simulate(ONE, IV, np.sort(shuffled), eff_delay).effect_conc
check("unsorted times work with an effect compartment", np.allclose(a[np.argsort(shuffled, kind="stable")], b))

TWO = PKParameters(ka=1.3, volume=20.0, clearance=4.0, peripheral_volume=35.0, intercompartmental_clearance=9.0, bioavailability=0.7)
s2 = simulate(TWO, Regimen(dose=150.0), t, eff_delay)
check("two-compartment PK with an effect compartment keeps its peripheral amount", s2.peripheral is not None and s2.effect_conc is not None)
check("two-compartment effect site peaks after plasma", t[np.argmax(s2.effect_conc)] > t[np.argmax(s2.conc)])
mm = PKParameters(ka=1.0, volume=20.0, elimination="michaelis_menten", vmax=10.0, km=2.0)
s3 = simulate(mm, Regimen(dose=200.0), t, eff_delay)
check("Michaelis-Menten elimination works with an effect compartment", np.all(np.isfinite(s3.effect)) and s3.effect_conc.max() > 0)

# ------------------------------------------------------------------ validation
check("ke0 must be positive", raises(lambda: EffectParameters(ke0=0.0)))
check("a negative ke0 is rejected", raises(lambda: EffectParameters(ke0=-1.0)))
check("ke0 None means direct, a number means delayed", not EffectParameters().delayed and EffectParameters(ke0=0.3).delayed)
check("equilibration half-life is ln2/ke0", abs(equilibration_half_life(0.5) - math.log(2) / 0.5) < 1e-12)
check("equilibration half-life rejects ke0 <= 0", raises(lambda: equilibration_half_life(0.0)))

# ------------------------------------------------------------------ fitting
inputs = wf.ModelInputs(
    dose=200.0, ka=1.2, volume=30.0, clearance=5.0,
    use_effect=True, delayed_effect=True, ke0=0.3, e0=10.0, emax=60.0, ec50=2.0, hill=1.0,
)
t_pd = np.array([0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24], dtype=float)
pk, reg, truth = wf.build_model(inputs)
y_true = simulate(pk, reg, t_pd, truth).effect


def make_obs(y):
    return parse_effect_observations("time,effect\n" + "\n".join(f"{x},{v:.10g}" for x, v in zip(t_pd, y)))


obs = make_obs(y_true)
free4 = ("e0", "emax", "ec50", "ke0")


def errs(result):
    return {n: abs(v - getattr(truth, n)) / abs(getattr(truth, n)) for n, v in result.estimates.items()}


res, drive = run_effect_fit(replace(inputs, ke0=1.0, e0=0.0, emax=1.0, ec50=1.0), obs, free4)
check("delayed fit recovers E0, Emax, EC50 and ke0 to 0.1%", res.success and max(errs(res).values()) < 1e-3, f"({errs(res)})")
check("the fit returns the effect-site concentrations it used", np.allclose(drive, effect_site_concentration(pk, reg, t_pd, res.effect.ke0), rtol=1e-6))
check("fitted result carries ke0", res.effect.ke0 is not None and abs(res.effect.ke0 - 0.3) < 1e-3)

for ke_start in (0.02, 0.3, 8.0):
    r, _ = run_effect_fit(replace(inputs, ke0=ke_start, e0=0.0, emax=1.0, ec50=1.0), obs, free4, use_form_start=True)
    check(f"delayed fit converges from a form start of ke0 = {ke_start}", max(errs(r).values()) < 1e-3, f"({errs(r)})")

rng = np.random.default_rng(2)
noisy = make_obs(y_true + 1.5 * rng.standard_normal(y_true.size))
rn, _ = run_effect_fit(replace(inputs, ke0=1.0), noisy, free4)
check("noisy delayed fit stays within 25% of the true values", max(errs(rn).values()) < 0.25, f"({errs(rn)})")
check("noisy delayed fit reports uncertainty for ke0", rn.cv_percent is not None and rn.cv_percent["ke0"] is not None and rn.se["ke0"] > 0)

only = run_effect_fit(inputs, obs, ("ke0",), use_form_start=True)[0]
check("ke0 alone can be fitted with the other parameters fixed", abs(only.effect.ke0 - 0.3) < 1e-3 and only.effect.emax == 60.0)

fixed_ke0, _ = run_effect_fit(replace(inputs, ke0=0.3), obs, ("e0", "emax", "ec50"))
check("with ke0 not fitted the form's ke0 is kept", fixed_ke0.effect.ke0 == 0.3 and fixed_ke0.free == ("e0", "emax", "ec50"))
check("fitting the other parameters at the right ke0 recovers them", max(errs(fixed_ke0).values()) < 1e-3, f"({errs(fixed_ke0)})")

direct_obs = make_obs(effect_from_conc(simulate(pk, reg, t_pd).conc, EffectParameters(e0=10.0, emax=60.0, ec50=2.0)))
rd, cd = run_effect_fit(replace(inputs, delayed_effect=False), direct_obs, free4)
check("direct data push ke0 to its upper limit", "ke0" in rd.at_bound and rd.effect.ke0 >= 0.99 * DEFAULT_KE0_BOUNDS[1])
check("direct data still give the right EC50 and Emax", abs(rd.effect.ec50 - 2.0) / 2.0 < 0.02 and abs(rd.effect.emax - 60.0) / 60.0 < 0.02, f"({rd.estimates})")
check("a delayed fit works when the form's delayed box is off", rd.effect.ke0 is not None)

conc_plasma = simulate(pk, reg, t_pd).conc
check("a delayed fit needs a site function", raises(lambda: fit_effect(conc_plasma, y_true, free=("e0", "emax", "ec50", "ke0"))))
check("a start with ke0 also needs a site function", raises(lambda: fit_effect(conc_plasma, y_true, start=truth, free=("e0", "emax"))))
check("ke0 is limited to its bounds", fit_effect(
    conc_plasma, y_true, start=replace(truth, ke0=0.3), free=("e0", "emax", "ec50", "ke0"),
    bounds={"ke0": (0.5, 5.0)}, site=lambda k: effect_site_concentration(pk, reg, t_pd, k),
).effect.ke0 >= 0.5 - 1e-9)
check("ke0 bounds must be positive", raises(lambda: fit_effect(
    conc_plasma, y_true, start=replace(truth, ke0=0.3), free=("ke0",), bounds={"ke0": (0.0, 1.0)},
    site=lambda k: effect_site_concentration(pk, reg, t_pd, k))))
check("a direct fit still returns no ke0", fit_effect(conc_plasma, direct_obs.effect).effect.ke0 is None)

# ------------------------------------------------------------------ workflow glue
check("direct choices are the four original parameters", effect_free_choices() == ("e0", "emax", "ec50", "hill"))
check("a delayed effect adds ke0 to the choices", effect_free_choices(True) == ("e0", "emax", "ec50", "hill", "ke0"))
check("default form inputs stay direct", not wf.ModelInputs().delayed_effect and wf.build_model(wf.ModelInputs(use_effect=True))[2].ke0 is None)
check("the form reaches the model as ke0", wf.build_model(inputs)[2].ke0 == 0.3)
check("the delayed box off leaves the effect direct", wf.build_model(replace(inputs, delayed_effect=False))[2].ke0 is None)
check("form_effect carries ke0 only when delayed", form_effect(inputs).ke0 == 0.3 and form_effect(replace(inputs, delayed_effect=False)).ke0 is None)

c_drive = driving_concentrations_at(inputs, obs)
check("driving concentrations are effect-site values for a delayed effect", np.allclose(c_drive, effect_site_concentration(pk, reg, t_pd, 0.3)))
check("driving concentrations are plasma values for a direct effect", np.allclose(driving_concentrations_at(replace(inputs, delayed_effect=False), obs), conc_plasma))
c_p, pred = predict_effects(inputs, obs)
check("predict_effects reproduces the simulated delayed effect", np.max(np.abs(pred - y_true)) < 1e-6)

applied = apply_effect_fit(replace(inputs, delayed_effect=False, ke0=0.9), res)
check("applying a delayed fit switches the delayed effect on and copies ke0", applied.delayed_effect and abs(applied.ke0 - 0.3) < 1e-3 and applied.use_effect)
check("applying a direct fit leaves the delayed setting alone", not apply_effect_fit(replace(inputs, delayed_effect=False), fit_effect(conc_plasma, direct_obs.effect)).delayed_effect)

rows = effect_fit_table_rows(res)
check("fit table has a ke0 row in 1/h", [r[0] for r in rows] == ["E0", "Emax", "EC50", "ke0"] and rows[3][2] == "1/h")

notes = effect_fit_notes(res, drive)
check("delayed fit notes describe the effect compartment and half-life", has(notes, "Effect-compartment model") and has(notes, "equilibration half-life"))
check("delayed fit notes do not claim there is no delay", not has(notes, "no delay"))
check("delayed fit notes warn that a wrong PK model can look like a delay", has(notes, "wrong PK model"))
dnotes = effect_fit_notes(rd, cd)
check("a ke0 on its upper limit is explained as no detectable delay", has(dnotes, "no detectable delay"))
check("direct fit notes are unchanged", has(effect_fit_notes(fit_effect(conc_plasma, direct_obs.effect), conc_plasma), "no delay"))

# residuals: a direct model against delayed data leaves a pattern, the right model does not
good = effect_residuals_for(inputs, obs)
bad = effect_residuals_for(replace(inputs, delayed_effect=False), obs)
check("residuals against the true delayed model are tiny", good.rmse < 1e-6)
check("residuals against a direct model show the delay", bad.rmse > 100 * max(good.rmse, 1e-6) and bad.rmse > 1.0, f"({bad.rmse:.3f})")
signs = np.sign(bad.residual[np.abs(bad.residual) > 0.5])
changes = int(np.sum(signs[1:] != signs[:-1]))
check("the direct-model residuals are systematic, not scattered", changes <= 4, f"({changes} sign changes)")

# ------------------------------------------------------------------ window time, text, export
text = wf.effect_half_life_text(wf.build_model(inputs)[2])
check("half-life text for a delayed effect", text.startswith("Effect-site equilibration half-life") and "ke0" in text)
check("half-life text is empty for a direct effect", wf.effect_half_life_text(None) == "" and wf.effect_half_life_text(EffectParameters()) == "")

sim = wf.simulate_inputs(inputs, n_points=20)
check("simulate_inputs carries the effect-site curve", sim.sim.effect_conc is not None)
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "sim.csv"
    wf.simulation_to_csv(path, sim)
    header = list(csv.reader(open(path)))[0]
    check("delayed export has an effect-site column before the effect",
          header == ["time_h", "conc_mg_per_L", "gut_mg", "central_mg", "effect_site_mg_per_L", "effect"])
    wf.simulation_to_csv(path, wf.simulate_inputs(replace(inputs, delayed_effect=False), n_points=20))
    check("direct export is unchanged", list(csv.reader(open(path)))[0] == ["time_h", "conc_mg_per_L", "gut_mg", "central_mg", "effect"])
    wf.simulation_to_csv(path, wf.simulate_inputs(replace(inputs, compartments=2), n_points=20))
    check("two-compartment delayed export has peripheral, effect site and effect",
          list(csv.reader(open(path)))[0] == ["time_h", "conc_mg_per_L", "gut_mg", "central_mg", "peripheral_mg", "effect_site_mg_per_L", "effect"])

# end to end: two-compartment IV model with a delayed effect
two_in = replace(inputs, compartments=2, route="iv_bolus", dose=100.0, volume=20.0, clearance=4.0,
                 peripheral_volume=35.0, intercompartmental_clearance=9.0, ke0=0.6)
pk2, reg2, eff2 = wf.build_model(two_in)
y2 = simulate(pk2, reg2, t_pd, eff2).effect
r2, _ = run_effect_fit(replace(two_in, ke0=2.0, e0=0.0, emax=1.0, ec50=1.0), make_obs(y2), free4)
check("end to end: IV two-compartment delayed effect is recovered",
      abs(r2.effect.ke0 - 0.6) / 0.6 < 1e-3 and abs(r2.effect.ec50 - 2.0) / 2.0 < 1e-3, f"({r2.estimates})")

print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
