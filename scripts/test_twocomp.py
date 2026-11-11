"""Checks for the two-compartment model, IV bolus dosing and their workflow glue.

    python scripts/test_twocomp.py

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
from ligentry.pkpd.effect_workflow import run_effect_fit
from ligentry.pkpd.fit import FitResult, fit_concentration
from ligentry.pkpd.model import (
    EffectParameters,
    PKParameters,
    Regimen,
    effect_from_conc,
    half_lives,
    hybrid_rate_constants,
    simulate,
)
from ligentry.pkpd.observations import parse_observations
from ligentry.pkpd.residuals import residuals_for

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


# ------------------------------------------------------------------ closed forms
V1, V2, Q, CL, KA, F, D = 20.0, 35.0, 9.0, 4.0, 1.3, 0.7, 150.0
TWO = PKParameters(
    ka=KA, volume=V1, clearance=CL, bioavailability=F,
    peripheral_volume=V2, intercompartmental_clearance=Q,
)
K10, K12, K21 = CL / V1, Q / V1, Q / V2
ALPHA, BETA = hybrid_rate_constants(TWO)


def iv_exact(t, dose=D):
    t = np.asarray(t, dtype=float)
    return dose / V1 * (
        (ALPHA - K21) / (ALPHA - BETA) * np.exp(-ALPHA * t)
        + (K21 - BETA) / (ALPHA - BETA) * np.exp(-BETA * t)
    )


def oral_exact(t, dose=D):
    t = np.asarray(t, dtype=float)
    return F * dose * KA / V1 * (
        (K21 - ALPHA) / ((KA - ALPHA) * (BETA - ALPHA)) * np.exp(-ALPHA * t)
        + (K21 - BETA) / ((KA - BETA) * (ALPHA - BETA)) * np.exp(-BETA * t)
        + (K21 - KA) / ((ALPHA - KA) * (BETA - KA)) * np.exp(-KA * t)
    )


t = np.linspace(0, 60, 601)
iv = Regimen(dose=D, route="iv_bolus")

check("hybrid constants satisfy the sum rule", abs(ALPHA + BETA - (K10 + K12 + K21)) < 1e-12)
check("hybrid constants satisfy the product rule", abs(ALPHA * BETA - K10 * K21) < 1e-12)
check("alpha is larger than beta", ALPHA > BETA > 0)

err = rel(simulate(TWO, iv, t).conc, iv_exact(t))
check("IV bolus matches the biexponential closed form", err < 1e-7, f"(rel err {err:.2e})")
err = rel(simulate(TWO, Regimen(dose=D), t).conc, oral_exact(t))
check("oral dose matches the three-exponential closed form", err < 1e-7, f"(rel err {err:.2e})")

one = PKParameters(ka=KA, volume=V1, clearance=CL)
s1 = simulate(one, iv, t)
check("one-compartment IV bolus is a single exponential", rel(s1.conc, D / V1 * np.exp(-CL / V1 * t)) < 1e-7)
check("IV dose is in the central compartment at time zero", abs(simulate(TWO, iv, np.array([0.0])).central[0] - D) < 1e-9)
check("one-compartment results have no peripheral amount", s1.peripheral is None)
check("two-compartment results carry a peripheral amount", simulate(TWO, iv, t).peripheral is not None)

# ------------------------------------------------------------------ exposure and mass balance
t_long = np.linspace(0, 4000, 400001)
auc_oral = simulate(TWO, Regimen(dose=D), t_long).metrics()["auc"]
check("oral AUC equals F*dose/CL", abs(auc_oral - F * D / CL) / (F * D / CL) < 1e-4, f"({auc_oral:.4f})")
auc_iv = simulate(TWO, iv, t_long).metrics()["auc"]
check("IV AUC equals dose/CL", abs(auc_iv - D / CL) / (D / CL) < 1e-4, f"({auc_iv:.4f})")

tt = np.linspace(0, 30, 30001)
s = simulate(TWO, iv, tt)
auc_30 = float(np.trapezoid(s.conc, tt) if hasattr(np, "trapezoid") else np.trapz(s.conc, tt))
balance = s.central[-1] + s.peripheral[-1] + CL * auc_30
check("dose = central + peripheral + eliminated (mass balance)", abs(balance - D) / D < 1e-5, f"({balance:.5f})")
check("amounts never go negative", np.all(s.central >= -1e-9) and np.all(s.peripheral >= -1e-9))

# ------------------------------------------------------------------ limiting cases
tiny_q = replace(TWO, intercompartmental_clearance=1e-9)
check("Q -> 0 gives the one-compartment oral result",
      rel(simulate(tiny_q, Regimen(dose=D), t).conc, simulate(PKParameters(ka=KA, volume=V1, clearance=CL, bioavailability=F), Regimen(dose=D), t).conc) < 1e-6)
fast = replace(TWO, intercompartmental_clearance=1e4)
late = np.array([2.0, 5.0, 10.0, 20.0])
v_total = V1 + V2
check("very fast distribution acts as one compartment of volume V1 + V2",
      rel(simulate(fast, iv, late).conc, D / v_total * np.exp(-CL / v_total * late)) < 2e-3)

# ------------------------------------------------------------------ route
a = simulate(TWO, iv, t).conc
b = simulate(replace(TWO, ka=9.0, bioavailability=0.2), iv, t).conc
check("IV bolus ignores ka and bioavailability", np.allclose(a, b, rtol=1e-12, atol=1e-12))
check("IV bolus leaves nothing in the gut", np.all(simulate(TWO, iv, t).gut == 0.0))

reg_iv = Regimen(dose=100.0, interval=8.0, n_doses=3, route="iv_bolus")
tt = np.linspace(0, 40, 801)
sim = simulate(TWO, reg_iv, tt)
total = np.zeros_like(tt)
for d in reg_iv.dose_times():
    total += np.where(tt >= d, iv_exact(np.clip(tt - d, 0, None), 100.0), 0.0)
check("repeated IV doses obey superposition", rel(sim.conc, total) < 1e-7)
before = simulate(TWO, reg_iv, np.array([8.0 - 1e-6, 8.0]))
check("IV concentration jumps by dose/V1 at a dose time",
      abs((before.conc[1] - before.conc[0]) - 100.0 / V1) < 1e-3)

end_dose = simulate(one, Regimen(dose=100.0, interval=8.0, n_doses=2), np.array([0.0, 8.0]))
check("a dose exactly at the end of the window is included (oral gut amount)", end_dose.gut[1] > 0.99 * 100.0)
end_iv = simulate(TWO, Regimen(dose=100.0, interval=8.0, n_doses=2, route="iv_bolus"), np.array([7.999999, 8.0]))
check("a dose exactly at the end of the window is included (IV concentration)", end_iv.conc[1] - end_iv.conc[0] > 0.99 * 100.0 / V1)

# ------------------------------------------------------------------ Michaelis-Menten with a peripheral compartment
mm_lin = PKParameters(
    ka=KA, volume=V1, elimination="michaelis_menten", vmax=CL * 1e6, km=1e6,
    peripheral_volume=V2, intercompartmental_clearance=Q, bioavailability=F,
)
check("two-compartment Michaelis-Menten tends to linear for Km >> C",
      rel(simulate(mm_lin, Regimen(dose=D), t).conc, simulate(TWO, Regimen(dose=D), t).conc) < 1e-4)
sat = PKParameters(
    ka=1.0, volume=V1, elimination="michaelis_menten", vmax=5.0, km=1e-6,
    peripheral_volume=V2, intercompartmental_clearance=Q,
)
s = simulate(sat, Regimen(dose=500.0, route="iv_bolus"), np.array([6.0, 12.0]))
slope = ((s.central + s.peripheral)[1] - (s.central + s.peripheral)[0]) / 6.0
check("saturated elimination removes Vmax mg/h from the whole body", abs(slope + 5.0) / 5.0 < 5e-3, f"({slope:.4f})")

# ------------------------------------------------------------------ half-lives
lives = half_lives(TWO)
check("terminal half-life is ln2/beta", abs(lives["terminal"] - math.log(2) / BETA) < 1e-12)
check("distribution half-life is ln2/alpha", abs(lives["distribution"] - math.log(2) / ALPHA) < 1e-12)
check("one-compartment half-life is ln2*V/CL", abs(half_lives(one)["terminal"] - math.log(2) * V1 / CL) < 1e-12)
check("one-compartment model has no distribution half-life", half_lives(one)["distribution"] is None)
check("hybrid constants need two compartments", raises(lambda: hybrid_rate_constants(one)))

# ------------------------------------------------------------------ validation
check("a peripheral volume needs a clearance", raises(lambda: PKParameters(ka=1, volume=1, clearance=1, peripheral_volume=2.0)))
check("a clearance needs a peripheral volume", raises(lambda: PKParameters(ka=1, volume=1, clearance=1, intercompartmental_clearance=2.0)))
check("peripheral parameters must be positive", raises(lambda: replace(TWO, peripheral_volume=0.0)))
check("rejects an unknown route", raises(lambda: Regimen(dose=1.0, route="iv_drip")))
check("compartments property", one.compartments == 1 and TWO.compartments == 2)
shuffled = np.array([20.0, 0.0, 3.0, 48.0, 1.0, 3.0])
a = simulate(TWO, iv, shuffled).conc
b = simulate(TWO, iv, np.sort(shuffled)).conc
check("unsorted times work for the two-compartment model", np.allclose(a[np.argsort(shuffled, kind="stable")], b))

# ------------------------------------------------------------------ fitting
t_obs = np.array([0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48], dtype=float)
c_iv = simulate(TWO, iv, t_obs).conc
start = replace(TWO, volume=V1 * 1.4, clearance=CL * 0.7, peripheral_volume=V2 * 0.6, intercompartmental_clearance=Q * 1.5)
names = ("volume", "clearance", "peripheral_volume", "intercompartmental_clearance")
fit = fit_concentration(t_obs, c_iv, iv, start, free=names)
errs = {n: abs(v - getattr(TWO, n)) / getattr(TWO, n) for n, v in fit.estimates.items()}
check("IV two-compartment fit converges", fit.success)
check("IV fit recovers V, CL, V2 and Q to 0.1%", max(errs.values()) < 1e-3, f"({errs})")
check("fit result remembers the route", fit.route == "iv_bolus")
check("fit reports uncertainty for the peripheral parameters", fit.cv_percent is not None and "peripheral_volume" in fit.cv_percent)

rng = np.random.default_rng(7)
c_noisy = c_iv * (1 + 0.03 * rng.standard_normal(c_iv.size))
fit = fit_concentration(t_obs, c_noisy, iv, start, free=names)
errs = {n: abs(v - getattr(TWO, n)) / getattr(TWO, n) for n, v in fit.estimates.items()}
check("noisy IV fit recovers all four parameters to 25%", max(errs.values()) < 0.25, f"({errs})")

check("ka cannot be fitted with IV bolus dosing", raises(lambda: fit_concentration(t_obs, c_iv, iv, start, free=("ka", "volume"))))
check("peripheral parameters cannot be fitted without them", raises(
    lambda: fit_concentration(t_obs, c_iv, iv, one, free=("peripheral_volume",))))

oral = Regimen(dose=D)
t_oral = np.array([0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48, 72], dtype=float)
c_oral = simulate(TWO, oral, t_oral).conc
ostart = replace(TWO, ka=KA * 1.3, volume=V1 * 1.2, clearance=CL * 0.85,
                 peripheral_volume=V2 * 0.8, intercompartmental_clearance=Q * 1.2)
ofit = fit_concentration(t_oral, c_oral, oral, ostart, free=("ka", *names))
oerrs = {n: abs(v - getattr(TWO, n)) / getattr(TWO, n) for n, v in ofit.estimates.items()}
check("oral two-compartment fit recovers all five parameters to 1%", max(oerrs.values()) < 1e-2, f"({oerrs})")

# ------------------------------------------------------------------ workflow glue
check("free_choices keeps the one-compartment oral list", wf.free_choices("linear") == ("ka", "volume", "clearance"))
check("free_choices keeps the Michaelis-Menten list", wf.free_choices("michaelis_menten") == ("ka", "volume", "vmax", "km"))
check("two compartments add V2 and Q", wf.free_choices("linear", 2) == ("ka", "volume", "clearance", "peripheral_volume", "intercompartmental_clearance"))
check("IV bolus removes ka", "ka" not in wf.free_choices("linear", 2, "iv_bolus"))
check("free_choices rejects a bad route", raises(lambda: wf.free_choices("linear", 1, "sublingual")))
check("free_choices rejects three compartments", raises(lambda: wf.free_choices("linear", 3)))

inputs = wf.ModelInputs()
pk, reg, eff = wf.build_model(inputs)
check("default inputs still build a one-compartment oral model", pk.compartments == 1 and reg.route == "oral" and eff is None)
two_in = replace(inputs, compartments=2, route="iv_bolus")
pk2, reg2, _ = wf.build_model(two_in)
check("two-compartment inputs reach the model", pk2.compartments == 2 and pk2.peripheral_volume == inputs.peripheral_volume and reg2.route == "iv_bolus")
check("build_model rejects a zero peripheral volume", raises(lambda: wf.build_model(replace(two_in, peripheral_volume=0.0))))
check("build_model rejects three compartments", raises(lambda: wf.build_model(replace(inputs, compartments=3))))


def old_span(pkx, regx):
    last = float(regx.dose_times()[-1])
    return wf._nice_ceiling(last + 3.0 / pkx.ka + 5.0 * math.log(2.0) / pkx.ke)


pk1, reg1, _ = wf.build_model(inputs)
check("one-compartment oral time window is unchanged", wf.default_time_span(pk1, reg1) == old_span(pk1, reg1))
oral_two = wf.default_time_span(*wf.build_model(replace(inputs, compartments=2))[:2])
iv_two = wf.default_time_span(*wf.build_model(replace(inputs, compartments=2, route="iv_bolus"))[:2])
check("a peripheral compartment lengthens the time window", oral_two >= wf.default_time_span(pk1, reg1))
check("IV dosing does not need an absorption allowance", iv_two <= oral_two)

for label, extra in (
    ("linear oral", {}),
    ("linear IV two-compartment", {"compartments": 2, "route": "iv_bolus"}),
    ("linear oral two-compartment", {"compartments": 2}),
    ("Michaelis-Menten IV two-compartment", {"compartments": 2, "route": "iv_bolus", "elimination": "michaelis_menten"}),
    ("Michaelis-Menten oral", {"elimination": "michaelis_menten"}),
):
    res = wf.simulate_inputs(replace(inputs, **extra))
    c = res.sim.conc
    check(f"automatic window shows most of the elimination ({label})", c[-1] < 0.1 * c.max(), f"(C_end/Cmax = {c[-1] / c.max():.3f})")

check("half-life text for one compartment", wf.half_life_text(pk1).startswith("Terminal half-life"))
check("half-life text for two compartments", wf.half_life_text(TWO).startswith("Half-lives: distribution"))
check("half-life text flags the Michaelis-Menten tail", "low concentration" in wf.half_life_text(sat))

# ------------------------------------------------------------------ fit notes
def fake(params, free, route="oral", cv=None):
    return FitResult(
        params=params, free=tuple(free), estimates={n: getattr(params, n) for n in free},
        cv_percent=cv if cv is not None else {n: 5.0 for n in free},
        rmse=0.0, n_obs=20, success=True, message="ok", route=route,
    )


def has(notes, text):
    return any(text in n for n in notes)


n = wf.fit_notes(fake(TWO, names, "iv_bolus"))
check("IV fit notes ask for early samples", has(n, "early, rapidly falling"))
check("IV fit notes do not claim apparent values", not has(n, "apparent"))
n = wf.fit_notes(fake(TWO, ("ka", *names)))
check("oral two-compartment notes recommend IV data", has(n, "IV bolus data"))
check("oral two-compartment notes mention V2/F and Q/F", has(n, "V2/F"))
blurred = replace(TWO, clearance=2.0, intercompartmental_clearance=0.4, peripheral_volume=4.0)
a_, b_ = hybrid_rate_constants(blurred)
check("test setup: phases really are close together", a_ / b_ < 3.0, f"({a_ / b_:.2f})")
n = wf.fit_notes(fake(blurred, names, "iv_bolus"))
check("notes warn when the phases are not separated", has(n, "not well separated"))
check("well separated phases give no such warning", not has(wf.fit_notes(fake(TWO, names, "iv_bolus")), "not well separated"))
flip = replace(TWO, ka=0.01)
check("flip-flop warning compares ka with the terminal rate for two compartments", has(wf.fit_notes(fake(flip, ("ka", "volume"))), "terminal rate constant"))
n = wf.fit_notes(fake(replace(one, ka=0.01), ("ka", "volume")))
check("one-compartment flip-flop note is unchanged", has(n, "elimination rate constant CL/V") and has(n, "flip-flop"))
n = wf.fit_notes(fake(one, ("volume", "clearance")))
check("one-compartment apparent-value note is unchanged", has(n, "V and CL are apparent values (V/F and CL/F)"))
check("one-compartment fit gets no two-compartment notes", not has(n, "distribution phase"))
real = wf.run_fit(replace(inputs, route="iv_bolus", compartments=2, dose=D, volume=V1 * 1.2, clearance=CL * 0.8, peripheral_volume=V2 * 0.7, intercompartmental_clearance=Q * 1.3),
                  parse_observations("time,conc\n" + "\n".join(f"{x},{y:.10g}" for x, y in zip(t_obs, c_iv))),
                  wf.free_choices("linear", 2, "iv_bolus"))
check("run_fit works end to end for an IV two-compartment model", real.success and max(abs(v - getattr(TWO, k)) / getattr(TWO, k) for k, v in real.estimates.items()) < 1e-3)
applied = wf.apply_fit(replace(inputs, compartments=2, route="iv_bolus"), real)
check("apply_fit copies V2 and Q into the form", abs(applied.peripheral_volume - V2) / V2 < 1e-3 and abs(applied.intercompartmental_clearance - Q) / Q < 1e-3)
rows = wf.fit_table_rows(real)
check("fit table labels the new parameters", {r[0] for r in rows} == {"V", "CL", "V2", "Q"} and {r[2] for r in rows} == {"L", "L/h"})

# ------------------------------------------------------------------ export
res = wf.simulate_inputs(replace(inputs, compartments=2, use_effect=True), n_points=20)
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "sim.csv"
    wf.simulation_to_csv(path, res)
    rows = list(csv.reader(open(path)))
    check("two-compartment export has a peripheral column before the effect",
          rows[0] == ["time_h", "conc_mg_per_L", "gut_mg", "central_mg", "peripheral_mg", "effect"])
    res1 = wf.simulate_inputs(inputs, n_points=20)
    wf.simulation_to_csv(path, res1)
    check("one-compartment export is unchanged", list(csv.reader(open(path)))[0] == ["time_h", "conc_mg_per_L", "gut_mg", "central_mg"])

# ------------------------------------------------------------------ other modules still work with two compartments
ins = replace(inputs, compartments=2, route="iv_bolus", dose=D, volume=V1, clearance=CL, peripheral_volume=V2, intercompartmental_clearance=Q)
obs = parse_observations("time,conc\n" + "\n".join(f"{x},{y:.10g}" for x, y in zip(t_obs, c_iv)))
check("residuals work with a two-compartment model", residuals_for(ins, obs).rmse < 1e-6)
pd_truth = EffectParameters(e0=10.0, emax=60.0, ec50=2.0, hill=1.0)
c_pd = simulate(wf.build_model(ins)[0], wf.build_model(ins)[1], t_obs).conc
eobs = parse_effect_observations("time,effect\n" + "\n".join(f"{x},{y:.10g}" for x, y in zip(t_obs, effect_from_conc(c_pd, pd_truth))))
eres, _ = run_effect_fit(ins, eobs)
check("effect fitting takes its concentrations from the two-compartment model",
      abs(eres.effect.ec50 - 2.0) / 2.0 < 1e-3 and abs(eres.effect.emax - 60.0) / 60.0 < 1e-3)

print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
