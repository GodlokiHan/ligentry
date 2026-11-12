"""PK/PD simulator window: simulate a dosing regimen, compare it with measured
concentrations and effects, and fit the model parameters.

This module only gathers numbers from the widgets and draws results. All the
modelling lives in ligentry.pkpd (workflow.py, model.py, fit.py, observations.py,
effect_*.py), which is tested without Qt.

Run it on its own with:  python -m ligentry.gui.pkpd_window
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

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

from ligentry.pkpd import effect_workflow, workflow
from ligentry.pkpd.effect_data import EffectObservations, load_effect_observations
from ligentry.pkpd.effect_residuals import effect_residuals_for
from ligentry.pkpd.observations import Observations, ObservationError, load_observations
from ligentry.pkpd.residuals import Residuals, residuals_for

_INTRO = (
    "One- or two-compartment model with oral or IV bolus dosing and an Emax effect, "
    "direct or delayed through an effect compartment. It is a teaching and exploration "
    "tool: results are only as good as the parameters you enter, and nothing here is "
    "dosing advice."
)
_PARAMETER_BOXES = (  # fit checkbox label, parameter name
    ("ka", "ka"),
    ("V", "volume"),
    ("CL", "clearance"),
    ("Vmax", "vmax"),
    ("Km", "km"),
    ("V2", "peripheral_volume"),
    ("Q", "intercompartmental_clearance"),
)
_EFFECT_BOXES = (  # fit checkbox label, parameter name, ticked by default
    ("E0", "e0", True),
    ("Emax", "emax", True),
    ("EC50", "ec50", True),
    ("Hill", "hill", False),
)
_FIT_ERRORS = (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError)


class PKPDWindow(QMainWindow):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("ligentry - PK/PD simulator")
        self.resize(1300, 820)

        self.obs: Observations | None = None
        self.eff_obs: EffectObservations | None = None
        self.last_result: workflow.SimulationResult | None = None
        self.last_fit = None
        self.last_effect_fit = None
        self.last_dir = Path.home()

        d = workflow.ModelInputs()

        # ---- dosing ----------------------------------------------------------------
        self.dose_spin = self._spin(d.dose, 0.001, 1e6, 3, " mg")
        self.interval_spin = self._spin(d.interval, 0.01, 1e4, 3, " h")
        self.n_spin = QSpinBox()
        self.n_spin.setRange(1, 500)
        self.n_spin.setValue(d.n_doses)
        self.n_spin.valueChanged.connect(self.schedule_update)
        self.route_combo = QComboBox()
        self.route_combo.addItem("Oral", workflow.ORAL)
        self.route_combo.addItem("IV bolus", workflow.IV_BOLUS)
        self.route_combo.setToolTip(
            "Oral doses are absorbed with rate ka. An IV bolus goes straight into the central "
            "compartment: ka and bioavailability then do not apply."
        )
        self.route_combo.currentIndexChanged.connect(self._structure_changed)
        dosing = QGroupBox("Dosing")
        form = QFormLayout(dosing)
        form.addRow("Route", self.route_combo)
        form.addRow("Dose", self.dose_spin)
        form.addRow("Number of doses", self.n_spin)
        form.addRow("Interval between doses", self.interval_spin)

        # ---- pharmacokinetics ----------------------------------------------------------
        self.ka_spin = self._spin(d.ka, 0.001, 100.0, 4, " 1/h")
        self.volume_spin = self._spin(d.volume, 0.01, 1e5, 3, " L")
        self.comp_combo = QComboBox()
        self.comp_combo.addItem("One compartment", 1)
        self.comp_combo.addItem("Two compartments", 2)
        self.comp_combo.currentIndexChanged.connect(self._structure_changed)
        self.volume2_spin = self._spin(d.peripheral_volume, 0.01, 1e5, 3, " L")
        self.q_spin = self._spin(d.intercompartmental_clearance, 0.001, 1e5, 4, " L/h")
        self.elim_combo = QComboBox()
        self.elim_combo.addItem("Linear (clearance)", workflow.LINEAR)
        self.elim_combo.addItem("Michaelis-Menten (Vmax, Km)", workflow.MICHAELIS_MENTEN)
        self.elim_combo.currentIndexChanged.connect(self._structure_changed)
        self.clearance_spin = self._spin(d.clearance, 0.001, 1e5, 4, " L/h")
        self.vmax_spin = self._spin(d.vmax, 0.001, 1e6, 4, " mg/h")
        self.km_spin = self._spin(d.km, 0.0001, 1e5, 4, " mg/L")
        self.f_spin = self._spin(d.bioavailability, 0.01, 1.0, 3, "")
        self.f_spin.setToolTip("Fraction of each dose that reaches the gut (bioavailability F)")
        pk_box = QGroupBox("Pharmacokinetics")
        form = QFormLayout(pk_box)
        form.addRow("Compartments", self.comp_combo)
        form.addRow("Absorption rate ka", self.ka_spin)
        form.addRow("Volume V (central)", self.volume_spin)
        form.addRow("Peripheral volume V2", self.volume2_spin)
        form.addRow("Intercompartmental clearance Q", self.q_spin)
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
        self._effect_spins = {
            "e0": self.e0_spin,
            "emax": self.emax_spin,
            "ec50": self.ec50_spin,
            "hill": self.hill_spin,
        }
        self.delayed_check = QCheckBox("Delayed effect (effect compartment)")
        self.delayed_check.setToolTip(
            "The effect follows an effect-site concentration that chases the plasma "
            "concentration with rate constant ke0, instead of following plasma directly. "
            "Use it when the effect lags behind the concentration."
        )
        self.ke0_spin = self._spin(d.ke0, 0.001, 100.0, 4, " 1/h")
        self.ke0_spin.setToolTip("Effect-site rate constant ke0. The equilibration half-life is ln2/ke0.")
        self.ke0_spin.setEnabled(False)
        self.delayed_check.toggled.connect(self._delayed_changed)
        self.effect_box = QGroupBox("Effect (direct Emax model)")
        self.effect_box.setCheckable(True)
        self.effect_box.setChecked(False)
        self.effect_box.toggled.connect(self.schedule_update)
        form = QFormLayout(self.effect_box)
        form.addRow("Baseline E0", self.e0_spin)
        form.addRow("Emax (negative = inhibition)", self.emax_spin)
        form.addRow("EC50", self.ec50_spin)
        form.addRow("Hill coefficient", self.hill_spin)
        form.addRow(self.delayed_check)
        form.addRow("Effect-site rate ke0", self.ke0_spin)

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

        # ---- measured concentrations and fitting ---------------------------------------------
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
        lay.addWidget(QLabel("Parameters to fit (uses the dosing and starting values above):"))
        lay.addLayout(free_row)
        form = QFormLayout()
        form.addRow("Weighting", self.weight_combo)
        lay.addLayout(form)
        lay.addLayout(fit_row)
        lay.addWidget(self.fit_table)
        lay.addWidget(self.fit_info)

        # ---- measured effects and effect fitting -----------------------------------------------
        self.eff_data_label = QLabel("No effect data loaded")
        self.eff_data_label.setWordWrap(True)
        self.eff_load_btn = QPushButton("Load effect data...")
        self.eff_load_btn.setToolTip(
            "CSV with time (h) and effect columns. The concentration at each time comes from "
            "the PK model entered above, so fit or enter the PK parameters first."
        )
        self.eff_load_btn.clicked.connect(self.load_effect_data)
        self.eff_clear_btn = QPushButton("Clear effect data")
        self.eff_clear_btn.clicked.connect(self.clear_effect_data)
        eff_data_row = QHBoxLayout()
        eff_data_row.addWidget(self.eff_load_btn)
        eff_data_row.addWidget(self.eff_clear_btn)

        self.eff_free_boxes: dict[str, QCheckBox] = {}
        eff_free_row = QHBoxLayout()
        for label, name, ticked in _EFFECT_BOXES:
            box = QCheckBox(label)
            box.setChecked(ticked)
            self.eff_free_boxes[name] = box
            eff_free_row.addWidget(box)
        self.eff_ke0_box = QCheckBox("ke0")
        self.eff_ke0_box.setChecked(False)
        self.eff_ke0_box.setToolTip(
            "Also fit the effect-site rate constant ke0 (a delayed effect). Unless you tick "
            "\"Start from the values above\", ke0 starts from a scan over many values. "
            "Using the fit then switches the delayed effect on."
        )
        eff_free_row.addWidget(self.eff_ke0_box)
        self.eff_form_start = QCheckBox("Start from the values above (not guessed from the data)")
        self.eff_form_start.setToolTip(
            "By default the fitted parameters start from values guessed from the data and "
            "the parameters you leave unticked keep the values in the Effect box."
        )
        self.eff_fit_btn = QPushButton("Fit effect")
        self.eff_fit_btn.setEnabled(False)
        self.eff_fit_btn.clicked.connect(self.fit_effect_clicked)
        self.eff_apply_btn = QPushButton("Use fitted effect values")
        self.eff_apply_btn.setEnabled(False)
        self.eff_apply_btn.clicked.connect(self.apply_effect_fit_clicked)
        eff_fit_row = QHBoxLayout()
        eff_fit_row.addWidget(self.eff_fit_btn)
        eff_fit_row.addWidget(self.eff_apply_btn)

        self.eff_table = QTableWidget(0, 4)
        self.eff_table.setHorizontalHeaderLabels(["Parameter", "Estimate", "Unit", "\u00b1 SE"])
        self.eff_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.eff_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.eff_table.setMaximumHeight(140)
        self.eff_info = QLabel("")
        self.eff_info.setWordWrap(True)

        eff_box = QGroupBox("Effect data and fitting")
        lay = QVBoxLayout(eff_box)
        lay.addWidget(self.eff_data_label)
        lay.addLayout(eff_data_row)
        lay.addWidget(QLabel("Effect parameters to fit (the PK model above supplies the concentrations):"))
        lay.addLayout(eff_free_row)
        lay.addWidget(self.eff_form_start)
        lay.addLayout(eff_fit_row)
        lay.addWidget(self.eff_table)
        lay.addWidget(self.eff_info)

        self.export_btn = QPushButton("Export simulation CSV...")
        self.export_btn.clicked.connect(self.export_csv)

        # ---- left panel ------------------------------------------------------------------
        intro = QLabel(_INTRO)
        intro.setWordWrap(True)
        left = QWidget()
        left_lay = QVBoxLayout(left)
        for w in (intro, dosing, pk_box, self.effect_box, time_box, self.metrics_label,
                  self.status_label, data_box, eff_box, self.export_btn):
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

        self._structure_changed()
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

    def _structure_changed(self, *_args) -> None:
        """Enable only the boxes that apply to the chosen route, compartments and elimination."""
        linear = self.elim_combo.currentData() == workflow.LINEAR
        oral = self.route_combo.currentData() == workflow.ORAL
        two = self.comp_combo.currentData() == 2
        self.clearance_spin.setEnabled(linear)
        self.vmax_spin.setEnabled(not linear)
        self.km_spin.setEnabled(not linear)
        self.ka_spin.setEnabled(oral)
        self.f_spin.setEnabled(oral)
        self.volume2_spin.setEnabled(two)
        self.q_spin.setEnabled(two)
        self.free_boxes["clearance"].setEnabled(linear)
        self.free_boxes["vmax"].setEnabled(not linear)
        self.free_boxes["km"].setEnabled(not linear)
        self.free_boxes["ka"].setEnabled(oral)
        self.free_boxes["peripheral_volume"].setEnabled(two)
        self.free_boxes["intercompartmental_clearance"].setEnabled(two)
        self.schedule_update()

    def _delayed_changed(self, checked: bool) -> None:
        self.ke0_spin.setEnabled(checked)
        self.effect_box.setTitle(
            "Effect (Emax model, effect compartment)" if checked else "Effect (direct Emax model)"
        )
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
            route=self.route_combo.currentData(),
            compartments=self.comp_combo.currentData(),
            peripheral_volume=self.volume2_spin.value(),
            intercompartmental_clearance=self.q_spin.value(),
            delayed_effect=self.delayed_check.isChecked(),
            ke0=self.ke0_spin.value(),
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
        self.route_combo.setCurrentIndex(max(self.route_combo.findData(inputs.route), 0))
        self.comp_combo.setCurrentIndex(max(self.comp_combo.findData(inputs.compartments), 0))
        self.volume2_spin.setValue(inputs.peripheral_volume)
        self.q_spin.setValue(inputs.intercompartmental_clearance)
        self.ke0_spin.setValue(inputs.ke0)
        self.delayed_check.setChecked(inputs.delayed_effect)
        self.schedule_update()

    # ---- simulate and draw ---------------------------------------------------------------

    def update_plot(self) -> None:
        inputs = self.gather_inputs()
        data_end = 0.0
        if self.obs is not None and len(self.obs):
            data_end = max(data_end, float(self.obs.times.max()))
        if self.eff_obs is not None and len(self.eff_obs):
            data_end = max(data_end, float(self.eff_obs.times.max()))
        min_t = data_end * 1.05
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

        resid = None
        if self.obs is not None and len(self.obs):
            try:
                resid = residuals_for(inputs, self.obs)
            except (ValueError, RuntimeError):
                resid = None

        eresid = None
        if self.eff_obs is not None and len(self.eff_obs):
            try:
                eresid = effect_residuals_for(inputs, self.eff_obs)
            except (ValueError, RuntimeError):
                eresid = None

        m = result.metrics
        text = (
            f"Cmax {m['cmax']:.4g} mg/L at {m['tmax']:.3g} h\n"
            f"AUC(0-{result.t_end:g} h) {m['auc']:.4g} mg*h/L\n"
            f"{workflow.half_life_text(result.pk)}"
        )
        delay_text = workflow.effect_half_life_text(result.effect)
        if delay_text:
            text += "\n" + delay_text
        if resid is not None:
            text += f"\nAgainst the concentration data: RMSE {resid.rmse:.4g} mg/L ({len(resid)} points)"
        if eresid is not None:
            text += f"\nAgainst the effect data: RMSE {eresid.rmse:.4g} ({len(eresid)} points)"
        self.metrics_label.setText(text)
        self.redraw(result, resid, eresid)

    def redraw(
        self,
        result: workflow.SimulationResult,
        resid: Residuals | None = None,
        eresid: Residuals | None = None,
    ) -> None:
        sim = result.sim
        has_eff_data = self.eff_obs is not None and len(self.eff_obs) > 0
        panels = ["conc"]
        if resid is not None:
            panels.append("resid")
        if sim.effect is not None or has_eff_data:
            panels.append("effect")
        if eresid is not None:
            panels.append("eresid")
        ratios = {"conc": 3.0, "resid": 1.0, "effect": 1.5, "eresid": 1.0}

        self.fig.clear()
        axes = self.fig.subplots(
            len(panels), 1, sharex=True, squeeze=False,
            gridspec_kw={"height_ratios": [ratios[p] for p in panels]},
        )[:, 0]
        ax = dict(zip(panels, axes))

        ax_c = ax["conc"]
        ax_c.plot(sim.t, sim.conc, color="#1f77b4", lw=1.8, label="model")
        if sim.effect_conc is not None:
            ax_c.plot(sim.t, sim.effect_conc, color="#2ca02c", lw=1.4, ls="--", label="effect site")
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

        if "resid" in ax:
            ax_r = ax["resid"]
            r = resid.residual
            ax_r.axhline(0.0, color="gray", lw=0.8)
            ax_r.vlines(resid.times, 0.0, r, color="#d62728", lw=0.8, alpha=0.6)
            ax_r.plot(resid.times, r, "o", color="#d62728", ms=4)
            limit = float(np.max(np.abs(r))) * 1.25
            if limit > 0:
                ax_r.set_ylim(-limit, limit)
            ax_r.set_ylabel("Residual\n(mg/L)")
            ax_r.grid(True, alpha=0.3)

        if "effect" in ax:
            ax_e = ax["effect"]
            if sim.effect is not None:
                ax_e.plot(sim.t, sim.effect, color="#2ca02c", lw=1.8, label="model")
            if has_eff_data:
                ax_e.plot(self.eff_obs.times, self.eff_obs.effect, "o", color="#9467bd", ms=5,
                          label=f"data ({len(self.eff_obs)} points)")
            if sim.effect is None:
                ax_e.text(0.02, 0.94, "Tick the Effect box to draw the model curve",
                          transform=ax_e.transAxes, va="top", fontsize=8, color="gray")
            ax_e.set_ylabel("Effect")
            ax_e.grid(True, alpha=0.3)
            handles, _labels = ax_e.get_legend_handles_labels()
            if handles:
                ax_e.legend(loc="upper right")

        if "eresid" in ax:
            ax_er = ax["eresid"]
            er = eresid.residual
            ax_er.axhline(0.0, color="gray", lw=0.8)
            ax_er.vlines(eresid.times, 0.0, er, color="#9467bd", lw=0.8, alpha=0.6)
            ax_er.plot(eresid.times, er, "o", color="#9467bd", ms=4)
            limit = float(np.max(np.abs(er))) * 1.25
            if limit > 0:
                ax_er.set_ylim(-limit, limit)
            ax_er.set_ylabel("Residual\n(effect)")
            ax_er.grid(True, alpha=0.3)

        axes[-1].set_xlabel("Time (h)")
        self.canvas.draw_idle()

    # ---- measured concentrations and fitting ------------------------------------------------

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
            QMessageBox.warning(self, "ligentry", str(exc))
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
            QMessageBox.information(self, "ligentry", "Tick at least one parameter to fit.")
            return
        inputs = self.gather_inputs()
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            result = workflow.run_fit(
                inputs, self.obs, free, weighting=self.weight_combo.currentData()
            )
        except (ValueError, RuntimeError) as exc:
            QMessageBox.warning(self, "ligentry", f"The fit could not be done:\n{exc}")
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

    # ---- measured effects and effect fitting -------------------------------------------------

    def load_effect_data(self) -> None:
        path, _filter = QFileDialog.getOpenFileName(
            self, "Load measured effects", str(self.last_dir),
            "Data files (*.csv *.tsv *.txt);;All files (*)",
        )
        if not path:
            return
        try:
            eff_obs = load_effect_observations(path)
        except ObservationError as exc:
            QMessageBox.warning(self, "ligentry", str(exc))
            return
        self.eff_obs = eff_obs
        self.last_dir = Path(path).parent
        self.eff_data_label.setText(f"{Path(path).name}: {len(eff_obs)} points")
        self.eff_fit_btn.setEnabled(True)
        self._reset_effect_fit_display()
        self.schedule_update()

    def clear_effect_data(self) -> None:
        self.eff_obs = None
        self.eff_data_label.setText("No effect data loaded")
        self.eff_fit_btn.setEnabled(False)
        self._reset_effect_fit_display()
        self.schedule_update()

    def _reset_effect_fit_display(self) -> None:
        self.last_effect_fit = None
        self.eff_table.setRowCount(0)
        self.eff_info.setText("")
        self.eff_apply_btn.setEnabled(False)

    def fit_effect_clicked(self) -> None:
        if self.eff_obs is None:
            return
        free = [n for n, box in self.eff_free_boxes.items() if box.isChecked()]
        if self.eff_ke0_box.isChecked():
            free.append("ke0")
        if not free:
            QMessageBox.information(self, "ligentry", "Tick at least one effect parameter to fit.")
            return
        inputs = self.gather_inputs()
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            result, conc = effect_workflow.run_effect_fit(
                inputs, self.eff_obs, free, use_form_start=self.eff_form_start.isChecked()
            )
        except _FIT_ERRORS as exc:
            QMessageBox.warning(self, "ligentry", f"The effect fit could not be done:\n{exc}")
            return
        finally:
            QApplication.restoreOverrideCursor()

        self.last_effect_fit = result
        rows = effect_workflow.effect_fit_table_rows(result)
        self.eff_table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, text in enumerate(row):
                self.eff_table.setItem(i, j, QTableWidgetItem(text))
        notes = effect_workflow.effect_fit_notes(result, conc)
        self.eff_info.setText(
            f"RMSE {result.rmse:.4g} over {result.n_obs} points."
            + ("\n" + "\n".join(notes) if notes else "")
        )
        self.eff_apply_btn.setEnabled(True)

    def apply_effect_fit_clicked(self) -> None:
        if self.last_effect_fit is None:
            return
        fitted = self.last_effect_fit.effect
        self.set_inputs(effect_workflow.apply_effect_fit(self.gather_inputs(), self.last_effect_fit))

        # the boxes have fixed ranges and decimals: say so if a fitted value did not fit in them
        shown = self.gather_inputs()
        clipped = []
        spins = dict(self._effect_spins)
        if fitted.ke0 is not None:
            spins["ke0"] = self.ke0_spin
        for name, spin in spins.items():
            want, got = getattr(fitted, name), getattr(shown, name)
            tolerance = 0.5 * 10.0 ** (-spin.decimals()) + 1e-3 * abs(want)
            if abs(want - got) > tolerance:
                label = effect_workflow.EFFECT_PARAMETER_INFO[name][0]
                clipped.append(f"{label} ({want:.4g} fitted, {got:.4g} shown)")
        if clipped:
            self.eff_info.setText(
                self.eff_info.text()
                + "\nThe box limits changed these values when they were copied: "
                + "; ".join(clipped) + "."
            )

    def export_csv(self) -> None:
        if self.last_result is None:
            QMessageBox.information(self, "ligentry", "Nothing to export yet.")
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
            QMessageBox.warning(self, "ligentry", f"Could not write the file:\n{exc}")


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    window = PKPDWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
