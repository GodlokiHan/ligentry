"""Tests for joint (simultaneous) PK/PD fitting and GUI controls (no display needed).

Run: python scripts/test_joint_fit.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

from ligentry.gui.pkpd_window import PKPDWindow
from ligentry.pkpd import joint_workflow, workflow as wf
from ligentry.pkpd.effect_data import parse_effect_observations
from ligentry.pkpd.joint_fit import fit_joint_pkpd
from ligentry.pkpd.model import EffectParameters, PKParameters, Regimen, simulate
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


# ------------------------------------------------------------------ 1. Direct joint fit (oral 1-comp)
truth_in = wf.ModelInputs(
    dose=100.0, ka=1.2, volume=30.0, clearance=5.0,
    use_effect=True, e0=5.0, emax=80.0, ec50=1.2, hill=1.0,
)
pk_t, reg_t, eff_t = wf.build_model(truth_in)
t_pk = np.array([0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0])
t_pd = np.array([0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 18.0, 24.0])

c_obs = parse_observations(
    "time,conc\n" + "\n".join(f"{t},{c:.10g}" for t, c in zip(t_pk, simulate(pk_t, reg_t, t_pk).conc))
)
e_obs = parse_effect_observations(
    "time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(t_pd, simulate(pk_t, reg_t, t_pd, eff_t).effect))
)

perturbed = replace(
    truth_in, ka=0.7, volume=45.0, clearance=3.5, e0=0.0, emax=50.0, ec50=2.5
)
res, drive = joint_workflow.run_joint_fit(
    perturbed, c_obs, e_obs,
    free_pk=("ka", "volume", "clearance"),
    free_pd=("e0", "emax", "ec50"),
)
check("direct joint fit converges", res.success, res.message)
check("direct joint fit recovers ka, V, CL to 1%",
      all(abs(res.pk_estimates[k] - getattr(pk_t, k)) / getattr(pk_t, k) < 0.01
          for k in ("ka", "volume", "clearance")), str(res.pk_estimates))
check("direct joint fit recovers E0, Emax, EC50 to 1%",
      abs(res.effect.e0 - 5.0) < 0.1
      and abs(res.effect.emax - 80.0) / 80.0 < 0.01
      and abs(res.effect.ec50 - 1.2) / 1.2 < 0.01, str(res.pd_estimates))
check(
    "direct joint fit has tiny RMSEs",
    res.pk_rmse < 1e-3 and res.pd_rmse < 1e-2,
    f"(pk_rmse={res.pk_rmse:.2e}, pd_rmse={res.pd_rmse:.2e})",
)

# ------------------------------------------------------------------ 2. Delayed joint fit (ke0)
delay_in = replace(truth_in, delayed_effect=True, ke0=0.35)
pk_d, reg_d, eff_d = wf.build_model(delay_in)
e_obs_d = parse_effect_observations(
    "time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(t_pd, simulate(pk_d, reg_d, t_pd, eff_d).effect))
)
res_d, drive_d = joint_workflow.run_joint_fit(
    replace(perturbed, delayed_effect=False, ke0=1.5),
    c_obs, e_obs_d,
    free_pk=("ka", "volume", "clearance"),
    free_pd=("e0", "emax", "ec50", "ke0"),
)
check("delayed joint fit converges", res_d.success, res_d.message)
check("delayed joint fit recovers ke0 to 1%",
      res_d.effect.ke0 is not None and abs(res_d.effect.ke0 - 0.35) / 0.35 < 0.01,
      f"got {res_d.effect.ke0}")
check("delayed joint fit recovers PK and PD together",
      abs(res_d.pk.clearance - 5.0) / 5.0 < 0.01 and abs(res_d.effect.ec50 - 1.2) / 1.2 < 0.01)

applied = joint_workflow.apply_joint_fit(perturbed, res_d)
check("apply_joint_fit updates PK, PD, use_effect and delayed_effect",
      applied.use_effect
      and applied.delayed_effect
      and abs(applied.clearance - 5.0) < 0.05
      and abs(applied.ke0 - 0.35) < 0.01)

# ------------------------------------------------------------------ 3. Noisy data & uncertainty
rng = np.random.default_rng(42)
c_noisy = parse_observations(
    "time,conc\n" + "\n".join(
        f"{t},{max(1e-4, c * (1 + 0.04 * rng.standard_normal())):.10g}"
        for t, c in zip(t_pk, simulate(pk_t, reg_t, t_pk).conc)
    )
)
e_noisy = parse_effect_observations(
    "time,effect\n" + "\n".join(
        f"{t},{e + 1.5 * rng.standard_normal():.10g}"
        for t, e in zip(t_pd, simulate(pk_d, reg_d, t_pd, eff_d).effect)
    )
)
res_n, drive_n = joint_workflow.run_joint_fit(
    perturbed, c_noisy, e_noisy,
    free_pk=("ka", "volume", "clearance"),
    free_pd=("e0", "emax", "ec50", "ke0"),
)
check("noisy joint fit stays within 25% of true parameters",
      abs(res_n.pk.clearance - 5.0) / 5.0 < 0.25
      and abs(res_n.effect.ec50 - 1.2) / 1.2 < 0.25
      and abs(res_n.effect.ke0 - 0.35) / 0.35 < 0.25)
check("joint fit reports SE and CV% for both PK and PD parameters",
      res_n.se is not None and res_n.cv_percent is not None
      and "clearance" in res_n.se and "ec50" in res_n.se and "ke0" in res_n.se)
rows = joint_workflow.joint_fit_table_rows(res_n)
check("table formats 7 unified rows", len(rows) == 7 and rows[0][0] == "ka" and rows[-1][0] == "ke0")
notes = joint_workflow.joint_fit_notes(res_n, drive_n)
check("joint notes mention effect-compartment fit",
      any("Joint effect-compartment fit" in n for n in notes))

# ------------------------------------------------------------------ 4. Validation errors
check("rejects empty free_pk", raises(lambda: fit_joint_pkpd(
    c_obs.times, c_obs.conc, e_obs.times, e_obs.effect, reg_t, pk_t, free_pk=(), free_pd=("e0",),
)))
check("rejects empty free_pd", raises(lambda: fit_joint_pkpd(
    c_obs.times, c_obs.conc, e_obs.times, e_obs.effect, reg_t, pk_t, free_pk=("clearance",), free_pd=(),
)))
check("rejects ka for IV bolus", raises(lambda: fit_joint_pkpd(
    c_obs.times, c_obs.conc, e_obs.times, e_obs.effect,
    Regimen(dose=100.0, route="iv_bolus"), pk_t, free_pk=("ka",), free_pd=("e0",),
)))

# ------------------------------------------------------------------ 5. GUI offscreen test
app = QApplication.instance() or QApplication(sys.argv)
popups: list[str] = []
QMessageBox.warning = staticmethod(lambda *a, **k: popups.append(str(a[2]) if len(a) > 2 else ""))
QMessageBox.information = staticmethod(lambda *a, **k: popups.append(str(a[2]) if len(a) > 2 else ""))

win = PKPDWindow()
check("joint button starts disabled until both datasets are loaded", not win.joint_fit_btn.isEnabled())

with tempfile.TemporaryDirectory() as tmp:
    p_pk = Path(tmp) / "pk.csv"
    p_pd = Path(tmp) / "pd.csv"
    p_pk.write_text("time,conc\n" + "\n".join(f"{t},{c:.10g}" for t, c in zip(c_obs.times, c_obs.conc)))
    p_pd.write_text("time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(e_obs_d.times, e_obs_d.effect)))

    chosen = {"path": str(p_pk)}
    QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (chosen["path"], ""))
    win.load_data()
    check("joint button stays disabled with only PK data", not win.joint_fit_btn.isEnabled())

    chosen["path"] = str(p_pd)
    win.load_effect_data()
    check("joint button enables when both PK and PD data are loaded", win.joint_fit_btn.isEnabled())

    win.set_inputs(replace(perturbed, delayed_effect=False))
    win.eff_ke0_box.setChecked(True)
    win.fit_joint_clicked()
    check("GUI joint fit runs without popups", not popups, str(popups))
    check("GUI joint table populates 7 rows", win.joint_table.rowCount() == 7)
    check("GUI joint apply button enables", win.joint_apply_btn.isEnabled())

    win.apply_joint_fit_clicked()
    out_in = win.gather_inputs()
    check("GUI apply copies PK + PD + delayed_effect into form",
          out_in.use_effect and out_in.delayed_effect and abs(out_in.ke0 - 0.35) < 0.01)

    win.clear_effect_data()
    check("clearing PD data disables joint button again", not win.joint_fit_btn.isEnabled())

win.close()
print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
