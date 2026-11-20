"""Tests for IV infusion dosing and oral lag time (Tlag) across model, fit, and GUI.

Run: python scripts/test_infusion_tlag.py
"""
from __future__ import annotations

import os
import sys
from dataclasses import replace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
from PySide6.QtWidgets import QApplication

from ligentry.gui.pkpd_window import PKPDWindow
from ligentry.pkpd import joint_workflow, workflow as wf
from ligentry.pkpd.effect_data import parse_effect_observations
from ligentry.pkpd.model import (
    PKParameters,
    Regimen,
    hybrid_rate_constants,
    simulate,
)
from ligentry.pkpd.observations import parse_observations

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


# ------------------------------------------------------------------ 1. IV infusion analytic checks
V, CL, D, TINF = 25.0, 5.0, 100.0, 2.0
K = CL / V
R0 = D / TINF
PK1 = PKParameters(ka=1.0, volume=V, clearance=CL)
REG_INF = Regimen(dose=D, route="iv_infusion", infusion_duration=TINF)

t = np.linspace(0.0, 24.0, 241)
sim1 = simulate(PK1, REG_INF, t)
exact1 = np.where(
    t <= TINF,
    (R0 / CL) * (1.0 - np.exp(-K * t)),
    (R0 / CL) * (1.0 - np.exp(-K * TINF)) * np.exp(-K * (t - TINF)),
)
err1 = float(np.max(np.abs(sim1.conc - exact1)) / np.max(exact1))
check("1-comp IV infusion matches closed form", err1 < 1e-7, f"(rel err {err1:.2e})")
check("IV infusion leaves gut at zero", np.all(sim1.gut == 0.0))
check("IV infusion peaks at end of infusion", abs(t[np.argmax(sim1.conc)] - TINF) < 1e-9)

# Short infusion approaches IV bolus
short_inf = Regimen(dose=D, route="iv_infusion", infusion_duration=1e-4)
bolus = Regimen(dose=D, route="iv_bolus")
t_post = np.array([0.1, 0.5, 1.0, 4.0, 12.0])
check(
    "very short IV infusion matches IV bolus",
    np.allclose(simulate(PK1, short_inf, t_post).conc, simulate(PK1, bolus, t_post).conc, rtol=1e-3),
)

# 2-comp IV infusion closed form
V2, Q = 40.0, 8.0
PK2 = PKParameters(ka=1.0, volume=V, clearance=CL, peripheral_volume=V2, intercompartmental_clearance=Q)
ALPHA, BETA = hybrid_rate_constants(PK2)
K21 = Q / V2
A_coef = (ALPHA - K21) / (V * (ALPHA - BETA))
B_coef = (K21 - BETA) / (V * (ALPHA - BETA))


def inf2_exact(tt: np.ndarray) -> np.ndarray:
    t_in = np.minimum(tt, TINF)
    t_post_inf = np.maximum(0.0, tt - TINF)
    term_a = (A_coef / ALPHA) * (1.0 - np.exp(-ALPHA * t_in)) * np.exp(-ALPHA * t_post_inf)
    term_b = (B_coef / BETA) * (1.0 - np.exp(-BETA * t_in)) * np.exp(-BETA * t_post_inf)
    return R0 * (term_a + term_b)


sim2 = simulate(PK2, REG_INF, t)
err2 = float(np.max(np.abs(sim2.conc - inf2_exact(t))) / np.max(inf2_exact(t)))
check("2-comp IV infusion matches biexponential closed form", err2 < 1e-7, f"(rel err {err2:.2e})")

# Overlapping multiple IV infusions (interval 3h, duration 5h)
reg_overlap = Regimen(dose=150.0, interval=3.0, n_doses=3, route="iv_infusion", infusion_duration=5.0)
t_ov = np.linspace(0.0, 24.0, 481)
sim_ov = simulate(PK1, reg_overlap, t_ov)
superpos = np.zeros_like(t_ov)
single_inf = Regimen(dose=150.0, route="iv_infusion", infusion_duration=5.0)
for d_t in reg_overlap.dose_times():
    superpos += np.where(
        t_ov >= d_t, simulate(PK1, single_inf, np.clip(t_ov - d_t, 0.0, None)).conc, 0.0
    )
check("overlapping IV infusions obey superposition", np.allclose(sim_ov.conc, superpos, rtol=1e-7, atol=1e-9))

# ------------------------------------------------------------------ 2. Oral lag time (Tlag)
TLAG = 1.25
reg_lag = Regimen(dose=100.0, route="oral", tlag=TLAG)
reg_nolag = Regimen(dose=100.0, route="oral", tlag=0.0)
t_oral = np.linspace(0.0, 20.0, 401)
sim_lag = simulate(PK1, reg_lag, t_oral)
check("concentration is zero before Tlag", np.all(sim_lag.conc[t_oral < TLAG] == 0.0))
t_shifted = np.array([2.0, 3.0, 5.0, 10.0])
check(
    "concentration after Tlag matches unlagged curve shifted by Tlag",
    np.allclose(
        simulate(PK1, reg_lag, t_shifted).conc,
        simulate(PK1, reg_nolag, t_shifted - TLAG).conc,
        rtol=1e-7,
    ),
)

# Validation errors
check("rejects non-positive infusion_duration", raises(lambda: Regimen(dose=100.0, route="iv_infusion", infusion_duration=0.0)))
check("rejects negative tlag", raises(lambda: Regimen(dose=100.0, route="oral", tlag=-0.1)))
check("free_choices excludes ka and tlag for iv_infusion",
      "ka" not in wf.free_choices("linear", 1, "iv_infusion")
      and "tlag" not in wf.free_choices("linear", 1, "iv_infusion"))
check("free_choices includes tlag for oral when requested",
      "tlag" not in wf.free_choices("linear", 1, "oral")
      and "tlag" in wf.free_choices("linear", 1, "oral", include_tlag=True))

# ------------------------------------------------------------------ 3. Fitting Tlag (sequential & joint)
truth_in = wf.ModelInputs(
    dose=100.0, ka=1.5, volume=25.0, clearance=5.0, route="oral", tlag=0.8,
    use_effect=True, e0=0.0, emax=100.0, ec50=1.5,
)
pk_t, reg_t, eff_t = wf.build_model(truth_in)
t_obs = np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0])
c_obs = parse_observations(
    "time,conc\n" + "\n".join(f"{x},{y:.10g}" for x, y in zip(t_obs, simulate(pk_t, reg_t, t_obs).conc))
)
e_obs = parse_effect_observations(
    "time,effect\n" + "\n".join(f"{x},{y:.10g}" for x, y in zip(t_obs, simulate(pk_t, reg_t, t_obs, eff_t).effect))
)

start_in = replace(truth_in, ka=1.0, volume=30.0, clearance=4.0, tlag=0.2, emax=70.0, ec50=2.0)
seq_fit = wf.run_fit(start_in, c_obs, ("ka", "volume", "clearance", "tlag"))
check("sequential fit recovers Tlag within 0.02 h", abs(seq_fit.estimates["tlag"] - 0.8) < 0.02, str(seq_fit.estimates))
applied_seq = wf.apply_fit(start_in, seq_fit)
check("apply_fit copies Tlag into ModelInputs", abs(applied_seq.tlag - 0.8) < 0.02)

j_fit, _ = joint_workflow.run_joint_fit(
    start_in, c_obs, e_obs,
    free_pk=("ka", "volume", "clearance", "tlag"),
    free_pd=("emax", "ec50"),
)
check("joint fit recovers Tlag and PD together",
      abs(j_fit.pk_estimates["tlag"] - 0.8) < 0.02 and abs(j_fit.effect.ec50 - 1.5) < 0.05)

# ------------------------------------------------------------------ 4. GUI controls
app = QApplication.instance() or QApplication(sys.argv)
win = PKPDWindow()
check("oral starts with tlag enabled and infusion disabled",
      win.tlag_spin.isEnabled() and not win.infusion_spin.isEnabled() and win.free_boxes["tlag"].isEnabled())

win.route_combo.setCurrentIndex(win.route_combo.findData("iv_infusion"))
check("iv_infusion enables infusion_spin and disables ka, F, tlag",
      win.infusion_spin.isEnabled()
      and not win.ka_spin.isEnabled()
      and not win.f_spin.isEnabled()
      and not win.tlag_spin.isEnabled()
      and not win.free_boxes["tlag"].isEnabled())

win.infusion_spin.setValue(2.5)
win.update_plot()
g = win.gather_inputs()
check("gather_inputs reads iv_infusion and infusion_duration",
      g.route == "iv_infusion" and abs(g.infusion_duration - 2.5) < 1e-9)

win.set_inputs(replace(wf.ModelInputs(), route="oral", tlag=0.65))
g2 = win.gather_inputs()
check("set_inputs round-trips oral tlag", g2.route == "oral" and abs(g2.tlag - 0.65) < 1e-9)

win.close()
print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
