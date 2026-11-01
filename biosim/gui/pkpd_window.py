"""PK/PD simulator window: simulate a dosing regimen, compare it with measured
concentrations and fit the model parameters.

This module only gathers numbers from the widgets and draws results. All the
modelling lives in biosim.pkpd (workflow.py, model.py, fit.py, observations.py),
which is tested without Qt.

Run it on its own with:  python -m biosim.gui.pkpd_window
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_API", "pyside6")  # before matplotlib chooses a Qt binding

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractSpinBox,
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from biosim.pkpd import workflow
from biosim.pkpd.observations import Observations, ObservationError, load_observations

_INTRO = (
    "One-compartment oral model with a direct Emax effect. It is a teaching and "
    "exploration tool: results are only as good as the parameters you enter, and "
    "nothing here is dosing advice."
)
_PARAMETER_BOXES = (  # fit checkbox label, parameter name
    ("ka", "ka"),
    ("V", "volume"),
    ("CL", "clearance"),
    ("Vmax", "vmax"),
    ("Km", "km"),
)


class PKPDWindow(QMainWindow):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("biosim - PK/PD simulator")
        self.resize(1300, 820)

        self.obs: Observations | None = None
        self.last_result: workflow.SimulationResult | None = None
        self.last_fit = None
        self.last_dir = Path.home()

        d = workflow.ModelInputs()

        # ---- dosing ----------------------------------------------------------------
        self.dose_spin = self._spin(d.dose, 0.001, 1e6, 3, " mg")
        self.interval_spin = self._spin(d.interval, 0.01, 1e4, 3, " h")
        self.n_spin = QSpinBox()
        self.n_spin.setRange(1, 500)
        self.n_spin.setValue(d.n_doses)
        self.n_spin.valueChanged.connect(self.schedule_update)
        dosing = QGroupBox("Dosing")
        form = QFormLayout(dosing)
        form.addRow("Dose", self.dose_spin)
        form.addRow("Number of doses", self.n_spin)
        form.addRow("Interval between doses", self.interval_spin)

        # ---- pharmacokinetics ----------------------------------------------------------
        self.ka_spin = self._spin(d.ka, 0.001, 100.0, 4, " 1/h")
        self.volume_spin = self._spin(d.volume, 0.01, 1e5, 3, " L")
        self.elim_combo = QComboBox()
        self.elim_combo.addItem("Linear (clearance)", workflow.LINEAR)
        self.elim_combo.addItem("Michaelis-Menten (Vmax, Km)", workflow.MICHAELIS_MENTEN)
        self.elim_combo.currentIndexChanged.connect(self._elimination_changed)
        self.clearance_spin = self._spin(d.clearance, 0.001, 1e5, 4, " L/h")
        self.vmax_spin = self._spin(d.vmax, 0.001, 1e6, 4, " mg/h")
        self.km_spin = self._spin(d.km, 0.0001, 1e5, 4, " mg/L")
        self.f_spin = self._spin(d.bioavailability, 0.01, 1.0, 3, "")
        self.f_spin.setToolTip("Fraction of each dose that reaches the gut (bioavailability F)")
        pk_box = QGroupBox("Pharmacokinetics")
        form = QFormLayout(pk_box)
        form.addRow("Absorption rate ka", self.ka_spin)
        form.addRow("Volume V", self.volume_spin)
        form.addRow("Elimination", self.elim_combo)
        form.addRow("Clearance CL", self.clearance_spin)
        form.addRow("Vmax", self.vmax_spin)
        form.addRow("Km", self.km_spin)
        form.addRow("Bioavailability F", self.f_spin)

        # ---- pharmacodynamics ------------------------------------------------------------
        self.e0_spin = self._spin(d.e0, -1e6, 1e6, 4, "")
        self.emax_spin = self._spin(d.emax, -1e6, 1e6, 4, "")
        self.ec50_spin = self._spin(d.ec50, 0.0001, 1e5, 4, " mg/L")
        self.hill_spin = self._spin(d.hill, 0.1, 10.0, 2, "")
        self.effect_box = QGroupBox("Effect (direct Emax model)")
        self.effect_box.setCheckable(True)
        self.effect_box.setChecked(False)
        self.effect_box.toggled.connect(self.schedule_update)
        form = QFormLayout(self.effect_box)
        form.addRow("Baseline E0", self.e0_spin)
        form.addRow("Emax (negative = inhibition)", self.emax_spin)
        form.addRow("EC50", self.ec50_spin)
        form.addRow("Hill coefficient", self.hill_spin)

        # ---- time window -----------------------------------------------------------------
        self.auto_time = QCheckBox("Choose the time window automatically")
        self.auto_time.setChecked(True)
        self.auto_time.toggled.connect(self._auto_time_changed)
        self.tend_spin = self._spin(36.0, 0.1, workflow.MAX_TIME_H, 2, " h")
        self.tend_spin.setEnabled(False)
        time_box = QGroupBox("Time window")
        form = QFormLayout(time_box)
        form.addRow(self.auto_time)
        form.addRow("End time", self.tend_spin)

        # ---- readout ---------------------------------------------------------------------
        self.metrics_label = QLabel("")
        self.metrics_label.setWordWrap(True)
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)

        # ---- measured data and fitting -----------------------------------------------------
        self.data_label = QLabel("No data loaded")
        self.data_label.setWordWrap(True)
        self.load_btn = QPushButton("Load data...")
        self.load_btn.setToolTip("CSV with time (h) and concentration (mg/L) columns")
        self.load_btn.clicked.connect(self.load_data)
        self.clear_btn = QPushButton("Clear data")
        self.clear_btn.clicked.connect(self.clear_data)
        data_row = QHBoxLayout()
        data_row.addWidget(self.load_btn)
        data_row.addWidget(self.clear_btn)

        self.free_boxes: dict[str, QCheckBox] = {}
        free_row = QHBoxLayout()
        for label, name in _PARAMETER_BOXES:
            box = QCheckBox(label)
            box.setChecked(True)
            self.free_boxes[name] = box
            free_row.addWidget(box)
        self.weight_combo = QComboBox()
        self.weight_combo.addItem("Proportional (relative error)", "proportional")
        self.weight_combo.addItem("None (absolute error)", "none")
        self.fit_btn = QPushButton("Fit parameters")
        self.fit_btn.setEnabled(False)
        self.fit_btn.clicked.connect(self.fit_clicked)
        self.apply_btn = QPushButton("Use fitted values")
        self.apply_btn.setEnabled(False)
        self.apply_btn.clicked.connect(self.apply_fit_clicked)
        fit_row = QHBoxLayout()
        fit_row.addWidget(self.fit_btn)
        fit_row.addWidget(self.apply_btn)

        self.fit_table = QTableWidget(0, 4)
        self.fit_table.setHorizontalHeaderLabels(["Parameter", "Estimate", "Unit", "CV %"])
        self.fit_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.fit_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.fit_table.setMaximumHeight(140)
        self.fit_info = QLabel("")
        self.fit_info.setWordWrap(True)

        data_box = QGroupBox("Measured data and fitting")
        lay = QVBoxLayout(data_box)
        lay.addWidget(self.data_label)
        lay.addLayout(data_row)
        lay.addWidget(QLabel("Parameters to fit (starting from the values above):"))
        lay.addLayout(free_row)
        form = QFormLayout()
        form.addRow("Weighting", self.weight_combo)
        lay.addLayout(form)
        lay.addLayout(fit_row)
        lay.addWidget(self.fit_table)
        lay.addWidget(self.fit_info)

        self.export_btn = QPushButton("Export simulation CSV...")
        self.export_btn.clicked.connect(self.export_csv)

        # ---- left panel ------------------------------------------------------------------
        intro = QLabel(_INTRO)
        intro.setWordWrap(True)
        left = QWidget()
        left_lay = QVBoxLayout(left)
        for w in (intro, dosing, pk_box, self.effect_box, time_box, self.metrics_label,
                  self.status_label, data_box, self.export_btn):
            left_lay.addWidget(w)
        left_lay.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(left)
        scroll.setMinimumWidth(440)

        # ---- plot --------------------------------------------------------------------------
        self.fig = Figure(figsize=(7, 6), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.addWidget(self.toolbar)
        right_lay.addWidget(self.canvas, 1)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(scroll)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        self.setCentralWidget(splitter)

        # updates are debounced so dragging a value does not simulate on every step
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(120)
        self._timer.timeout.connect(self.update_plot)

        self._elimination_changed()
        self.update_plot()

    # ---- helpers ------------------------------------------------------------------------

    def _spin(self, value, lo, hi, decimals, suffix) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setDecimals(decimals)
        spin.setRange(lo, hi)
        spin.setValue(value)
        spin.setStepType(QAbstractSpinBox.StepType.AdaptiveDecimalStepType)
        spin.setKeyboardTracking(False)  # update when editing finishes, not on every key
        if suffix:
            spin.setSuffix(suffix)
        spin.valueChanged.connect(self.schedule_update)
        return spin

    def schedule_update(self, *_args) -> None:
        self._timer.start()

    def _elimination_changed(self, *_args) -> None:
        linear = self.elim_combo.currentData() == workflow.LINEAR
        self.clearance_spin.setEnabled(linear)
        self.vmax_spin.setEnabled(not linear)
        self.km_spin.setEnabled(not linear)
        self.free_boxes["clearance"].setEnabled(linear)
        self.free_boxes["vmax"].setEnabled(not linear)
        self.free_boxes["km"].setEnabled(not linear)
        self.schedule_update()

    def _auto_time_changed(self, checked: bool) -> None:
        self.tend_spin.setEnabled(not checked)
        self.schedule_update()

    def gather_inputs(self) -> workflow.ModelInputs:
        return workflow.ModelInputs(
            dose=self.dose_spin.value(),
            interval=self.interval_spin.value(),
            n_doses=self.n_spin.value(),
            ka=self.ka_spin.value(),
            volume=self.volume_spin.value(),
            elimination=self.elim_combo.currentData(),
            clearance=self.clearance_spin.value(),
            vmax=self.vmax_spin.value(),
            km=self.km_spin.value(),
            bioavailability=self.f_spin.value(),
            use_effect=self.effect_box.isChecked(),
            e0=self.e0_spin.value(),
            emax=self.emax_spin.value(),
            ec50=self.ec50_spin.value(),
            hill=self.hill_spin.value(),
            t_end=None if self.auto_time.isChecked() else self.tend_spin.value(),
        )

    def set_inputs(self, inputs: workflow.ModelInputs) -> None:
        self.dose_spin.setValue(inputs.dose)
        self.interval_spin.setValue(inputs.interval)
        self.n_spin.setValue(inputs.n_doses)
        self.ka_spin.setValue(inputs.ka)
        self.volume_spin.setValue(inputs.volume)
        self.elim_combo.setCurrentIndex(max(self.elim_combo.findData(inputs.elimination), 0))
        self.clearance_spin.setValue(inputs.clearance)
        self.vmax_spin.setValue(inputs.vmax)
        self.km_spin.setValue(inputs.km)
        self.f_spin.setValue(inputs.bioavailability)
        self.effect_box.setChecked(inputs.use_effect)
        self.e0_spin.setValue(inputs.e0)
        self.emax_spin.setValue(inputs.emax)
        self.ec50_spin.setValue(inputs.ec50)
        self.hill_spin.setValue(inputs.hill)
        self.schedule_update()

    # ---- simulate and draw ---------------------------------------------------------------

    def update_plot(self) -> None:
        inputs = self.gather_inputs()
        min_t = 0.0
        if self.obs is not None and len(self.obs):
            min_t = float(self.obs.times.max()) * 1.05
        try:
            result = workflow.simulate_inputs(inputs, min_t_end=min_t)
        except (ValueError, RuntimeError) as exc:
            self.status_label.setStyleSheet("color: #c0392b;")
            self.status_label.setText(f"Cannot simulate: {exc}")
            self.metrics_label.setText("")
            return
        self.status_label.setStyleSheet("")
        self.status_label.setText("")
        self.last_result = result

        if self.auto_time.isChecked():
            self.tend_spin.blockSignals(True)  # showing the automatic value must not re-trigger
            self.tend_spin.setValue(result.t_end)
            self.tend_spin.blockSignals(False)

        m = result.metrics
        self.metrics_label.setText(
            f"Cmax {m['cmax']:.4g} mg/L at {m['tmax']:.3g} h\n"
            f"AUC(0-{result.t_end:g} h) {m['auc']:.4g} mg*h/L"
        )
        self.redraw(result)

    def redraw(self, result: workflow.SimulationResult) -> None:
        sim = result.sim
        self.fig.clear()
        if sim.effect is not None:
            ax_c, ax_e = self.fig.subplots(
                2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]}
            )
        else:
            ax_c, ax_e = self.fig.subplots(), None

        ax_c.plot(sim.t, sim.conc, color="#1f77b4", lw=1.8, label="model")
        if result.regimen.n_doses <= 40:
            for t in result.regimen.dose_times():
                ax_c.axvline(t, color="gray", lw=0.7, ls=":")
        if self.obs is not None and len(self.obs):
            ax_c.plot(self.obs.times, self.obs.conc, "o", color="#d62728", ms=5,
                      label=f"data ({len(self.obs)} points)")
        ax_c.set_ylabel("Concentration (mg/L)")
        ax_c.set_ylim(bottom=0)
        ax_c.set_xlim(0, result.t_end)
        ax_c.grid(True, alpha=0.3)
        ax_c.legend(loc="upper right")

        if ax_e is not None:
            ax_e.plot(sim.t, sim.effect, color="#2ca02c", lw=1.8)
            ax_e.set_ylabel("Effect")
            ax_e.set_xlabel("Time (h)")
            ax_e.grid(True, alpha=0.3)
        else:
            ax_c.set_xlabel("Time (h)")
        self.canvas.draw_idle()

    # ---- measured data and fitting ----------------------------------------------------------

    def load_data(self) -> None:
        path, _filter = QFileDialog.getOpenFileName(
            self, "Load measured concentrations", str(self.last_dir),
            "Data files (*.csv *.tsv *.txt);;All files (*)",
        )
        if not path:
            return
        try:
            obs = load_observations(path)
        except ObservationError as exc:
            QMessageBox.warning(self, "biosim", str(exc))
            return
        self.obs = obs
        self.last_dir = Path(path).parent
        self.data_label.setText(f"{Path(path).name}: {len(obs)} points")
        self.fit_btn.setEnabled(True)
        self._reset_fit_display()
        self.schedule_update()

    def clear_data(self) -> None:
        self.obs = None
        self.data_label.setText("No data loaded")
        self.fit_btn.setEnabled(False)
        self._reset_fit_display()
        self.schedule_update()

    def _reset_fit_display(self) -> None:
        self.last_fit = None
        self.fit_table.setRowCount(0)
        self.fit_info.setText("")
        self.apply_btn.setEnabled(False)

    def fit_clicked(self) -> None:
        if self.obs is None:
            return
        free = [n for n, box in self.free_boxes.items() if box.isEnabled() and box.isChecked()]
        if not free:
            QMessageBox.information(self, "biosim", "Tick at least one parameter to fit.")
            return
        inputs = self.gather_inputs()
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            result = workflow.run_fit(
                inputs, self.obs, free, weighting=self.weight_combo.currentData()
            )
        except (ValueError, RuntimeError) as exc:
            QMessageBox.warning(self, "biosim", f"The fit could not be done:\n{exc}")
            return
        finally:
            QApplication.restoreOverrideCursor()

        self.last_fit = result
        rows = workflow.fit_table_rows(result)
        self.fit_table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, text in enumerate(row):
                self.fit_table.setItem(i, j, QTableWidgetItem(text))
        notes = workflow.fit_notes(result)
        self.fit_info.setText(
            f"RMSE {result.rmse:.4g} mg/L over {result.n_obs} points."
            + ("\n" + "\n".join(notes) if notes else "")
        )
        self.apply_btn.setEnabled(True)

    def apply_fit_clicked(self) -> None:
        if self.last_fit is None:
            return
        self.set_inputs(workflow.apply_fit(self.gather_inputs(), self.last_fit))

    def export_csv(self) -> None:
        if self.last_result is None:
            QMessageBox.information(self, "biosim", "Nothing to export yet.")
            return
        path, _filter = QFileDialog.getSaveFileName(
            self, "Export simulation", str(self.last_dir / "pkpd_simulation.csv"),
            "CSV files (*.csv)",
        )
        if not path:
            return
        try:
            workflow.simulation_to_csv(path, self.last_result)
        except OSError as exc:
            QMessageBox.warning(self, "biosim", f"Could not write the file:\n{exc}")


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    window = PKPDWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
