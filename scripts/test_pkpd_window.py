"""Drives the PK/PD window without showing it (Qt "offscreen" platform).

    python scripts/test_pkpd_window.py

File dialogs and message boxes are replaced by stand-ins, so nothing opens and
nothing waits for a click. Covers the model-structure controls (route,
compartments, elimination), simulating, loading data, fitting, applying a fit and
exporting. Needs PySide6 and matplotlib, but no display.
"""
from __future__ import annotations

import csv
import os
import sys
import tempfile
from pathlib import Path

os.environ["QT_QPA_PLATFORM"] = "offscreen"

import numpy as np

from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

from ligentry.gui.pkpd_window import PKPDWindow
from ligentry.pkpd import workflow
from ligentry.pkpd.model import PKParameters, Regimen, simulate

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


app = QApplication.instance() or QApplication(sys.argv)

messages: list[str] = []
QMessageBox.warning = staticmethod(lambda *a, **k: messages.append(f"warning: {a[2] if len(a) > 2 else ''}"))
QMessageBox.information = staticmethod(lambda *a, **k: messages.append(f"info: {a[2] if len(a) > 2 else ''}"))

chosen: dict[str, str] = {}
QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (chosen.get("open", ""), ""))
QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (chosen.get("save", ""), ""))

win = PKPDWindow()

# ---------------------------------------------------------------- defaults
check("window starts one-compartment and oral", win.comp_combo.currentData() == 1 and win.route_combo.currentData() == "oral")
check("peripheral boxes are disabled for one compartment", not win.volume2_spin.isEnabled() and not win.q_spin.isEnabled())
check("peripheral fit boxes are disabled for one compartment",
      not win.free_boxes["peripheral_volume"].isEnabled() and not win.free_boxes["intercompartmental_clearance"].isEnabled())
check("ka and F are enabled for oral dosing", win.ka_spin.isEnabled() and win.f_spin.isEnabled())
win.update_plot()
check("one-compartment readout shows the terminal half-life", "Terminal half-life" in win.metrics_label.text())
check("default simulation is one-compartment", win.last_result.pk.compartments == 1 and win.last_result.sim.peripheral is None)

# ---------------------------------------------------------------- structure controls
win.comp_combo.setCurrentIndex(win.comp_combo.findData(2))
check("two compartments enable V2 and Q", win.volume2_spin.isEnabled() and win.q_spin.isEnabled())
check("two compartments enable their fit boxes",
      win.free_boxes["peripheral_volume"].isEnabled() and win.free_boxes["intercompartmental_clearance"].isEnabled())
win.update_plot()
check("two-compartment readout shows both half-lives", "Half-lives: distribution" in win.metrics_label.text())
check("two-compartment simulation has a peripheral amount", win.last_result.sim.peripheral is not None)

win.route_combo.setCurrentIndex(win.route_combo.findData("iv_bolus"))
check("IV bolus disables ka and F", not win.ka_spin.isEnabled() and not win.f_spin.isEnabled())
check("IV bolus disables the ka fit box", not win.free_boxes["ka"].isEnabled())
win.update_plot()
check("IV simulation uses the IV route", win.last_result.regimen.route == "iv_bolus")
check("IV concentration at time zero is dose / V1",
      abs(win.last_result.sim.conc[0] - win.dose_spin.value() / win.volume_spin.value()) < 1e-6)

win.elim_combo.setCurrentIndex(win.elim_combo.findData("michaelis_menten"))
check("Michaelis-Menten enables Vmax and Km, disables CL", win.vmax_spin.isEnabled() and win.km_spin.isEnabled() and not win.clearance_spin.isEnabled())
win.update_plot()
check("Michaelis-Menten readout marks the low-concentration tail", "low concentration" in win.metrics_label.text())
check("Michaelis-Menten two-compartment simulation runs", win.status_label.text() == "")
win.elim_combo.setCurrentIndex(win.elim_combo.findData("linear"))

# ---------------------------------------------------------------- inputs round trip
inputs = win.gather_inputs()
check("gather_inputs reads the route and compartments", inputs.route == "iv_bolus" and inputs.compartments == 2)
changed = workflow.ModelInputs(
    dose=250.0, route="oral", compartments=1, peripheral_volume=77.0, intercompartmental_clearance=3.5,
)
win.set_inputs(changed)
back = win.gather_inputs()
check("set_inputs sets route and compartments", back.route == "oral" and back.compartments == 1)
check("set_inputs sets V2 and Q", abs(back.peripheral_volume - 77.0) < 1e-9 and abs(back.intercompartmental_clearance - 3.5) < 1e-9)

# ---------------------------------------------------------------- load, fit, apply
truth = PKParameters(ka=1.0, volume=20.0, clearance=4.0, peripheral_volume=35.0, intercompartmental_clearance=9.0)
reg = Regimen(dose=150.0, route="iv_bolus")
times = np.array([0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48], dtype=float)
conc = simulate(truth, reg, times).conc
with tempfile.TemporaryDirectory() as tmp:
    data_path = Path(tmp) / "iv.csv"
    data_path.write_text("time,conc\n" + "\n".join(f"{t},{c:.10g}" for t, c in zip(times, conc)))

    win.set_inputs(workflow.ModelInputs(
        dose=150.0, route="iv_bolus", compartments=2, volume=26.0, clearance=3.0,
        peripheral_volume=22.0, intercompartmental_clearance=12.0,
    ))
    chosen["open"] = str(data_path)
    win.load_data()
    check("loading data enables the fit button", win.fit_btn.isEnabled() and "16 points" in win.data_label.text())
    win.update_plot()
    check("the readout reports the fit against the data", "Against the concentration data" in win.metrics_label.text())

    messages.clear()
    win.fit_clicked()
    check("the fit ran without a warning dialog", not messages, f"({messages})")
    check("the fit table lists V, CL, V2 and Q", win.fit_table.rowCount() == 4)
    labels = {win.fit_table.item(i, 0).text() for i in range(win.fit_table.rowCount())}
    check("ka is not fitted for IV dosing", labels == {"V", "CL", "V2", "Q"}, f"({labels})")
    check("the fit recovers the true values", win.last_fit is not None and all(
        abs(v - getattr(truth, k)) / getattr(truth, k) < 1e-3 for k, v in win.last_fit.estimates.items()))
    check("the notes mention the distribution phase", "distribution phase" in win.fit_info.text())
    check("applying becomes possible", win.apply_btn.isEnabled())
    win.apply_fit_clicked()
    applied = win.gather_inputs()
    check("applying copies V2 and Q into the form",
          abs(applied.peripheral_volume - 35.0) < 0.05 and abs(applied.intercompartmental_clearance - 9.0) < 0.05)
    win.update_plot()
    check("after applying, the model matches the data", "RMSE 0" in win.metrics_label.text() or "e-" in win.metrics_label.text())

    # ---------------------------------------------------------------- effect data with the two-compartment model
    eff_path = Path(tmp) / "pd.csv"
    from ligentry.pkpd.model import EffectParameters, effect_from_conc
    pd_truth = EffectParameters(e0=10.0, emax=60.0, ec50=2.0)
    eff_path.write_text("time,effect\n" + "\n".join(f"{t},{e:.10g}" for t, e in zip(times, effect_from_conc(conc, pd_truth))))
    chosen["open"] = str(eff_path)
    win.load_effect_data()
    messages.clear()
    win.fit_effect_clicked()
    check("the effect fit ran without a warning dialog", not messages, f"({messages})")
    check("the effect fit recovers EC50 and Emax", win.last_effect_fit is not None
          and abs(win.last_effect_fit.effect.ec50 - 2.0) < 0.01 and abs(win.last_effect_fit.effect.emax - 60.0) < 0.1)

    # ---------------------------------------------------------------- export
    win.effect_box.setChecked(True)
    win.update_plot()
    out_path = Path(tmp) / "sim.csv"
    chosen["save"] = str(out_path)
    win.export_csv()
    header = next(csv.reader(open(out_path)))
    check("the export has a peripheral column", "peripheral_mg" in header and header[-1] == "effect", f"({header})")

    # ---------------------------------------------------------------- errors are shown, not raised
    win.set_inputs(workflow.ModelInputs(compartments=2, peripheral_volume=60.0, intercompartmental_clearance=8.0))
    win.volume2_spin.setValue(win.volume2_spin.minimum())
    win.update_plot()
    check("a tiny peripheral volume still simulates", win.last_result is not None)

print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
