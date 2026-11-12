"""Drives the PK/PD window's delayed-effect (effect compartment) controls without a screen.

    python scripts/test_pkpd_window_delayed.py

Runs Qt offscreen (QT_QPA_PLATFORM=offscreen), with file dialogs and message boxes
replaced so nothing waits for a click.
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
from ligentry.pkpd import workflow as wf
from ligentry.pkpd.model import simulate

PASSED = 0
FAILED = 0
POPUPS: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASSED, FAILED
    if condition:
        PASSED += 1
        print(f"ok   - {name}")
    else:
        FAILED += 1
        print(f"FAIL - {name} {detail}")


app = QApplication.instance() or QApplication([])
QMessageBox.warning = staticmethod(lambda *a, **k: POPUPS.append(str(a[2]) if len(a) > 2 else ""))
QMessageBox.information = staticmethod(lambda *a, **k: POPUPS.append(str(a[2]) if len(a) > 2 else ""))

win = PKPDWindow()
win.show()
app.processEvents()


def labels(ax) -> list[str]:
    return [line.get_label() for line in ax.get_lines()]


# ------------------------------------------------------------------ defaults
d = win.gather_inputs()
check("the window starts with a direct effect", not d.delayed_effect)
check("the default ke0 matches the model default", d.ke0 == wf.ModelInputs().ke0)
check("the default inputs equal ModelInputs()", d == wf.ModelInputs())
check("the ke0 box starts disabled", not win.ke0_spin.isEnabled())
check("the effect box title says direct", "direct" in win.effect_box.title())
check("the ke0 fit box starts unticked", not win.eff_ke0_box.isChecked())

# ------------------------------------------------------------------ the checkbox
win.delayed_check.setChecked(True)
check("ticking the delayed box enables ke0", win.ke0_spin.isEnabled())
check("the effect box title mentions the effect compartment", "effect compartment" in win.effect_box.title())
win.ke0_spin.setValue(0.8)
g = win.gather_inputs()
check("gather_inputs reads the delayed box and ke0", g.delayed_effect and abs(g.ke0 - 0.8) < 1e-9)
win.delayed_check.setChecked(False)
check("unticking disables ke0 again and restores the title", not win.ke0_spin.isEnabled() and "direct" in win.effect_box.title())

win.set_inputs(replace(wf.ModelInputs(), use_effect=True, delayed_effect=True, ke0=0.25))
g = win.gather_inputs()
check("set_inputs round-trips the delayed effect", g.delayed_effect and abs(g.ke0 - 0.25) < 1e-9 and g.use_effect)

# ------------------------------------------------------------------ drawing
win.update_plot()
res = win.last_result
check("a delayed simulation carries the effect-site curve", res is not None and res.sim.effect_conc is not None)
check("the concentration panel draws the effect site", "effect site" in labels(win.fig.axes[0]))
check("the readout shows the equilibration half-life", "equilibration half-life" in win.metrics_label.text())

win.set_inputs(replace(win.gather_inputs(), delayed_effect=False))
win.update_plot()
check("a direct effect draws no effect-site curve", "effect site" not in labels(win.fig.axes[0]))
check("a direct effect shows no equilibration half-life", "equilibration" not in win.metrics_label.text())

win.set_inputs(replace(win.gather_inputs(), use_effect=False, delayed_effect=True))
win.update_plot()
check("with the effect box off nothing is delayed", win.last_result.sim.effect_conc is None)

# ------------------------------------------------------------------ fitting ke0 from effect data
truth_in = wf.ModelInputs(
    dose=200.0, ka=1.2, volume=30.0, clearance=5.0,
    use_effect=True, delayed_effect=True, ke0=0.3, e0=10.0, emax=60.0, ec50=2.0, hill=1.0,
)
t_pd = np.array([0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24], dtype=float)
pk, reg, eff = wf.build_model(truth_in)
y = simulate(pk, reg, t_pd, eff).effect
tmp = tempfile.TemporaryDirectory()
path = Path(tmp.name) / "effect.csv"
path.write_text("time,effect\n" + "\n".join(f"{a},{b:.10g}" for a, b in zip(t_pd, y)))
QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (str(path), ""))

win.set_inputs(replace(truth_in, delayed_effect=False, ke0=1.0, e0=0.0, emax=1.0, ec50=1.0))
win.load_effect_data()
app.processEvents()
win.update_plot()
check("effect data are loaded", win.eff_obs is not None and len(win.eff_obs) == len(t_pd))
check("the effect residual panel appears", len(win.fig.axes) == 3, f"({len(win.fig.axes)} panels)")

win.eff_ke0_box.setChecked(True)
win.fit_effect_clicked()
fit = win.last_effect_fit
check("no warning popup during the fit", not POPUPS, f"({POPUPS})")
check("the window fits ke0 when its box is ticked", fit is not None and "ke0" in fit.free)
check("the fitted ke0 is right", fit is not None and abs(fit.effect.ke0 - 0.3) / 0.3 < 1e-3, f"({fit.effect if fit else None})")
check("the fit table lists ke0 in 1/h", win.eff_table.rowCount() == 4 and win.eff_table.item(3, 0).text() == "ke0" and win.eff_table.item(3, 2).text() == "1/h")
check("the fit notes describe the effect compartment", "Effect-compartment model" in win.eff_info.text())

win.apply_effect_fit_clicked()
g = win.gather_inputs()
check("using the fit switches the delayed effect on", win.delayed_check.isChecked() and g.delayed_effect)
check("using the fit copies ke0 into the box", abs(g.ke0 - 0.3) / 0.3 < 1e-3, f"({g.ke0})")
win.update_plot()
check("after using the fit the model matches the data", win.last_result is not None and "effect site" in labels(win.fig.axes[0]))

# a direct fit (ke0 box unticked) still works and leaves the delayed setting alone
win.set_inputs(replace(win.gather_inputs(), delayed_effect=False))
win.eff_ke0_box.setChecked(False)
win.fit_effect_clicked()
check("a fit without ke0 does not include it", win.last_effect_fit is not None and "ke0" not in win.last_effect_fit.free and win.last_effect_fit.effect.ke0 is None)
win.apply_effect_fit_clicked()
check("using a direct fit leaves the delayed box off", not win.delayed_check.isChecked())

win.close()
tmp.cleanup()
print(f"\n{PASSED}/{PASSED + FAILED} checks passed")
sys.exit(0 if FAILED == 0 else 1)
