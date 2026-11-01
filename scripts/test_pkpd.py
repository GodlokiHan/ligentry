"""Checks for ligentry.pkpd against closed-form solutions and parameter recovery.

    python scripts/test_pkpd.py

Needs only numpy and scipy. No network, no RDKit.
"""
from __future__ import annotations

import sys
from dataclasses import replace

import numpy as np

from ligentry.pkpd.fit import fit_concentration
from ligentry.pkpd.model import (
    EffectParameters,
    PKParameters,
    Regimen,
    effect_from_conc,
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


def bateman(t, dose, f, ka, v, cl):
    """Closed-form plasma concentration for a single oral dose, linear elimination."""
    ke = cl / v
    t = np.asarray(t, dtype=float)
    return f * dose * ka / (v * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t))


BASE = PKParameters(ka=1.2, volume=30.0, clearance=5.0, bioavailability=0.8)

# 1. single dose against the closed form
t = np.linspace(0, 48, 400)
sim = simulate(BASE, Regimen(dose=200.0), t)
exact = bateman(t, 200.0, 0.8, 1.2, 30.0, 5.0)
err = np.max(np.abs(sim.conc - exact)) / exact.max()
check("single oral dose matches the Bateman equation", err < 1e-6, f"(rel err {err:.2e})")

# 2. AUC to infinity equals F * dose / CL
t_long = np.linspace(0, 400, 40001)
auc = simulate(BASE, Regimen(dose=200.0), t_long).metrics()["auc"]
auc_exact = 0.8 * 200.0 / 5.0
check("AUC equals F*dose/CL", abs(auc - auc_exact) / auc_exact < 1e-4, f"({auc:.4f} vs {auc_exact:.4f})")

# 3. repeated dosing is the sum of shifted single-dose curves
reg = Regimen(dose=100.0, interval=8.0, n_doses=4)
t = np.linspace(0, 60, 601)
sim = simulate(BASE, reg, t)
total = np.zeros_like(t)
for d in reg.dose_times():
    shifted = np.clip(t - d, 0, None)
    total += np.where(t >= d, bateman(shifted, 100.0, 0.8, 1.2, 30.0, 5.0), 0.0)
err = np.max(np.abs(sim.conc - total)) / total.max()
check("repeated dosing obeys superposition", err < 1e-6, f"(rel err {err:.2e})")

# 4. concentration is continuous at a dose time, the gut amount jumps
at = np.array([7.999999, 8.0, 8.000001])
s = simulate(BASE, reg, at)
check("concentration is continuous across a dose", np.ptp(s.conc) < 1e-4)
check("gut amount jumps at a dose time", s.gut[1] - s.gut[0] > 0.5 * 0.8 * 100.0)

# 5. late first dose: nothing before it
late = simulate(BASE, Regimen(dose=100.0, start=5.0), np.array([0.0, 4.9, 5.0, 6.0, 12.0]))
check("nothing in the body before the first dose", np.all(late.conc[:2] == 0.0))
check("drug appears after a late first dose", late.conc[3] > 0 and late.conc[4] > 0)

# 6. times in any order give the same answer
shuffled = np.array([12.0, 0.0, 3.0, 48.0, 1.0, 3.0])
a = simulate(BASE, Regimen(dose=200.0), shuffled).conc
b = simulate(BASE, Regimen(dose=200.0), np.sort(shuffled)).conc
check("unsorted times are handled", np.allclose(a[np.argsort(shuffled, kind="stable")], b))

# 7. Michaelis-Menten reduces to linear when Km is far above any concentration
mm_lin = PKParameters(
    ka=1.2, volume=30.0, elimination="michaelis_menten",
    vmax=5.0 * 1e6, km=1e6, bioavailability=0.8,
)
t = np.linspace(0, 48, 400)
lin = simulate(BASE, Regimen(dose=200.0), t).conc
mm = simulate(mm_lin, Regimen(dose=200.0), t).conc
err = np.max(np.abs(lin - mm)) / lin.max()
check("Michaelis-Menten tends to linear for Km >> C", err < 1e-4, f"(rel err {err:.2e})")

# 8. saturated Michaelis-Menten eliminates at a constant rate Vmax
sat = PKParameters(ka=50.0, volume=10.0, elimination="michaelis_menten", vmax=5.0, km=1e-6)
s = simulate(sat, Regimen(dose=100.0), np.array([5.0, 10.0]))
slope = (s.central[1] - s.central[0]) / 5.0
check("saturated elimination removes Vmax mg/h", abs(slope + 5.0) / 5.0 < 5e-3, f"(slope {slope:.4f})")

# 9. effect model
eff = EffectParameters(e0=2.0, emax=10.0, ec50=3.0, hill=1.0)
check("effect is E0 at zero concentration", float(effect_from_conc(0.0, eff)) == 2.0)
check("effect is E0 + Emax/2 at EC50", abs(float(effect_from_conc(3.0, eff)) - 7.0) < 1e-12)
eff2 = replace(eff, hill=2.5)
check("EC50 gives half effect for any Hill slope", abs(float(effect_from_conc(3.0, eff2)) - 7.0) < 1e-12)
check("effect saturates at E0 + Emax", abs(float(effect_from_conc(1e6, eff)) - 12.0) < 1e-3)
inh = EffectParameters(e0=100.0, emax=-80.0, ec50=2.0)
check("negative Emax models inhibition", abs(float(effect_from_conc(2.0, inh)) - 60.0) < 1e-12)
s = simulate(BASE, Regimen(dose=200.0), np.linspace(0, 24, 49), effect=eff)
check("simulate attaches an effect curve", s.effect is not None and s.effect.shape == s.conc.shape)
check("no effect curve without effect parameters", simulate(BASE, Regimen(dose=200.0), t).effect is None)

# 10. metrics
m = simulate(BASE, Regimen(dose=200.0), np.linspace(0, 48, 4801)).metrics()
ke = 5.0 / 30.0
tmax_exact = np.log(1.2 / ke) / (1.2 - ke)
check("Tmax matches the closed form", abs(m["tmax"] - tmax_exact) < 0.01, f"({m['tmax']:.3f} vs {tmax_exact:.3f})")

# 11. input validation
check("rejects a non-positive dose", raises(lambda: Regimen(dose=0.0)))
check("rejects multiple doses without an interval", raises(lambda: Regimen(dose=1.0, n_doses=3)))
check("rejects linear elimination without clearance", raises(lambda: PKParameters(ka=1.0, volume=1.0)))
check(
    "rejects Michaelis-Menten without vmax/km",
    raises(lambda: PKParameters(ka=1.0, volume=1.0, elimination="michaelis_menten")),
)
check("rejects negative times", raises(lambda: simulate(BASE, Regimen(dose=1.0), np.array([-1.0, 2.0]))))
check("rejects bad bioavailability", raises(lambda: replace(BASE, bioavailability=1.5)))

# 12. fit recovers known parameters from noise-free data
truth = PKParameters(ka=1.2, volume=30.0, clearance=5.0)
reg1 = Regimen(dose=200.0)
t_obs = np.array([0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48], dtype=float)
c_obs = simulate(truth, reg1, t_obs).conc
guess = PKParameters(ka=0.8, volume=40.0, clearance=3.5)
fit = fit_concentration(t_obs, c_obs, reg1, guess)
rel = {k: abs(v - getattr(truth, k)) / getattr(truth, k) for k, v in fit.estimates.items()}
check("fit converges on noise-free data", fit.success)
check("fit recovers ka, V and CL to 0.1%", max(rel.values()) < 1e-3, f"({rel})")
check("fit reports parameter uncertainty", fit.cv_percent is not None and set(fit.cv_percent) == set(fit.free))

# 13. fit with 5% proportional noise
rng = np.random.default_rng(1)
c_noisy = c_obs * (1 + 0.05 * rng.standard_normal(c_obs.size))
fit = fit_concentration(t_obs, c_noisy, reg1, guess)
rel = {k: abs(v - getattr(truth, k)) / getattr(truth, k) for k, v in fit.estimates.items()}
check("noisy fit recovers parameters to 15%", max(rel.values()) < 0.15, f"({rel})")

# 14. bounds are respected
bounded = fit_concentration(t_obs, c_obs, reg1, guess, bounds={"ka": (0.5, 0.6)})
check("bounds keep ka inside its range", bounded.estimates["ka"] <= 0.6 * (1 + 1e-6))

# 15. fit with a subset of free parameters leaves the rest alone
part = fit_concentration(t_obs, c_obs, reg1, replace(truth, ka=1.2), free=("volume", "clearance"))
check("fixed parameters are left unchanged", part.params.ka == 1.2)

# 16. Michaelis-Menten fit (ka and V known)
mm_truth = PKParameters(ka=1.5, volume=20.0, elimination="michaelis_menten", vmax=8.0, km=2.0)
reg2 = Regimen(dose=300.0)
t2 = np.array([0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24], dtype=float)
c2 = simulate(mm_truth, reg2, t2).conc
mm_guess = replace(mm_truth, vmax=5.0, km=4.0)
mm_fit = fit_concentration(t2, c2, reg2, mm_guess, free=("vmax", "km"))
rel = {k: abs(v - getattr(mm_truth, k)) / getattr(mm_truth, k) for k, v in mm_fit.estimates.items()}
check("fit recovers Michaelis-Menten Vmax and Km", max(rel.values()) < 1e-2, f"({rel})")

# 17. fit input validation
check("fit rejects negative concentrations", raises(lambda: fit_concentration(t_obs, -c_obs, reg1, guess)))
check("fit rejects an unknown parameter", raises(lambda: fit_concentration(t_obs, c_obs, reg1, guess, free=("foo",))))
check("fit rejects vmax on linear elimination", raises(lambda: fit_concentration(t_obs, c_obs, reg1, guess, free=("vmax",))))
check("fit rejects too few observations", raises(lambda: fit_concentration(t_obs[:2], c_obs[:2], reg1, guess)))

print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
