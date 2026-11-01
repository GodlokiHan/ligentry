"""biosim desktop application (Qt)."""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import psutil
from PySide6.QtCore import QProcess, QSettings, Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QFileDialog, QFormLayout, QHBoxLayout,
    QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit,
    QProgressBar, QPushButton, QSpinBox, QSplitter, QStyle, QTableWidget,
    QTableWidgetItem, QToolBar, QVBoxLayout, QWidget,
)

from biosim.gui import runio
from biosim.gui.viewer import MolViewer

_PROGRESS_RE = re.compile(r"\[(\d+)/(\d+)\]")  # "[3/20] name" from the docking log

COL_RANK, COL_NAME, COL_SCORE, COL_ISO, COL_STATUS, COL_NOTES = range(6)
HEADERS = ["Rank", "Name", "Score (kcal/mol)", "Isomers", "Status", "Notes"]
MAX_RECENT = 8


def fmt_time(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.resize(1400, 850)
        self.proc: QProcess | None = None  # docking runs in its own process (biosim.gui.job)
        self._buffer = ""  # partial output line from the job
        self._last_error = ""
        self.workdir = Path("runs/gui")
        self.pdb_id = ""
        self.rows: list[dict] = []  # table row i <-> self.rows[i]
        self.dirty = False
        self.running = False
        self.stop_requested = False
        self._populating = False

        # run-state used by the timers
        self.t_start: float | None = None
        self.t_first_ligand: float | None = None
        self.completed = 0
        self.total = 0
        self.failed = False

        self.pdb_edit = QLineEdit("1HSG")
        self.chembl_edit = QLineEdit("CHEMBL243")
        self.n_spin = QSpinBox()
        self.n_spin.setRange(1, 200)
        self.n_spin.setValue(5)
        self.smiles_edit = QPlainTextEdit()
        self.smiles_edit.setPlaceholderText(
            "Or paste ligands, one per line: SMILES [name]\n"
            "(used only when the ChEMBL field is empty)"
        )
        self.smiles_edit.setMaximumHeight(110)
        self.exh_spin = QSpinBox()
        self.exh_spin.setRange(1, 64)
        self.exh_spin.setValue(8)
        self.stereo_check = QCheckBox("Enumerate unspecified stereoisomers")
        self.stereo_check.setChecked(True)
        self.stereo_check.setToolTip(
            "Ligands whose SMILES leave stereocentres or double bonds unspecified are "
            "docked as every isomer (up to the cap) and the best-scoring one is kept. "
            "Slower: up to the cap times more docking per such ligand."
        )
        self.iso_spin = QSpinBox()
        self.iso_spin.setRange(1, 64)
        self.iso_spin.setValue(8)
        self.iso_spin.setToolTip("Most isomers docked per ligand; more are sampled at random.")
        self.stereo_check.toggled.connect(self.iso_spin.setEnabled)
        self.out_edit = QLineEdit("runs/gui")
        self.run_btn = QPushButton("Run docking")
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setToolTip(
            "Finish the current docking step, keep the results so far and stop")
        self.stop_btn.setEnabled(False)

        form = QFormLayout()
        form.addRow("PDB ID (with bound ligand)", self.pdb_edit)
        form.addRow("ChEMBL target", self.chembl_edit)
        form.addRow("Max ligands", self.n_spin)
        form.addRow("Ligand SMILES", self.smiles_edit)
        form.addRow("Exhaustiveness", self.exh_spin)
        form.addRow("Stereoisomers", self.stereo_check)
        form.addRow("Max isomers per ligand", self.iso_spin)
        form.addRow("Output folder", self.out_edit)

        buttons = QHBoxLayout()
        buttons.addWidget(self.run_btn, 3)
        buttons.addWidget(self.stop_btn, 1)

        # progress, timers, CPU meter
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.progress.setFormat("%v / %m ligands")
        self.progress.setAlignment(Qt.AlignCenter)

        self.elapsed_label = QLabel("Elapsed 00:00")
        self.remaining_label = QLabel("Remaining -")
        times = QHBoxLayout()
        times.addWidget(self.elapsed_label)
        times.addStretch(1)
        times.addWidget(self.remaining_label)

        self.cpu_bar = QProgressBar()
        self.cpu_bar.setRange(0, 100)
        self.cpu_bar.setFormat("CPU %p%")
        self.cpu_bar.setAlignment(Qt.AlignCenter)

        self.table = QTableWidget(0, len(HEADERS))
        self.table.setHorizontalHeaderLabels(HEADERS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        # only the Notes cells are editable (the other items drop the edit flag)
        self.table.setEditTriggers(
            QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed
        )
        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setSectionResizeMode(COL_NAME, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(COL_SCORE, QHeaderView.ResizeToContents)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setMaximumHeight(160)

        left = QWidget()
        lay = QVBoxLayout(left)
        lay.addLayout(form)
        lay.addLayout(buttons)
        lay.addWidget(self.progress)
        lay.addLayout(times)
        lay.addWidget(self.cpu_bar)
        lay.addWidget(self.table, 1)
        lay.addWidget(self.log)

        self.viewer = MolViewer(self)
        viewer_widget = self.viewer if isinstance(self.viewer, QWidget) else self.viewer.interactor

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left)
        splitter.addWidget(viewer_widget)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        self.setCentralWidget(splitter)

        self._build_actions()
        self._build_menus_and_toolbar()
        self.statusBar().showMessage("Ready")

        self.run_btn.clicked.connect(self.start)
        self.stop_btn.clicked.connect(self.stop)
        self.table.itemSelectionChanged.connect(self.on_selection_changed)
        self.table.itemChanged.connect(self.on_item_changed)

        psutil.cpu_percent(interval=None)  # prime: the first call always returns 0
        self.timer = QTimer(self)
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.tick)
        self.timer.start()

        self.update_title()
        self.update_actions()

    # ---- actions, menus, toolbar -----------------------------------------

    def _action(self, text, slot, shortcut=None, icon=None, tip="", checkable=False):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(self.style().standardIcon(icon))
        if shortcut:
            action.setShortcut(QKeySequence(shortcut))
        if tip:
            action.setStatusTip(tip)
            action.setToolTip(f"{text}: {tip}")
        action.setCheckable(checkable)
        if checkable:
            action.toggled.connect(slot)
        else:
            action.triggered.connect(lambda _checked=False, s=slot: s())
        return action

    def _build_actions(self) -> None:
        S = QStyle.StandardPixmap
        self.act_new = self._action(
            "New run", self.new_run, "Ctrl+N", S.SP_FileIcon, "Clear the results and start fresh")
        self.act_open = self._action(
            "Open run...", self.open_run, "Ctrl+O", S.SP_DirOpenIcon,
            "Load a finished run folder without docking again")
        self.act_save = self._action(
            "Save", self.save, "Ctrl+S", S.SP_DialogSaveButton,
            "Write the results (with your notes) back to the run folder")
        self.act_save_as = self._action(
            "Save as...", self.save_as, "Ctrl+Shift+S", S.SP_DriveHDIcon,
            "Copy the whole run folder to a new place and continue there")
        self.act_export = self._action(
            "Export CSV...", self.export_csv, "Ctrl+E", S.SP_ArrowDown,
            "Write the table to a CSV file anywhere")
        self.act_quit = self._action("Quit", self.close, "Ctrl+Q", S.SP_DialogCloseButton)

        self.act_copy_smiles = self._action(
            "Copy SMILES", self.copy_smiles, "Ctrl+Shift+C", S.SP_FileLinkIcon,
            "Copy the selected ligand's SMILES (best isomer) to the clipboard")
        self.act_edit_note = self._action(
            "Edit note", self.edit_note, "F2", S.SP_FileDialogDetailedView,
            "Edit the note of the selected ligand")
        self.act_remove = self._action(
            "Remove ligand", self.remove_selected, "Ctrl+Backspace", S.SP_TrashIcon,
            "Remove the selected ligand from the results table (files stay on disk)")

        self.act_run = self._action(
            "Run docking", self.start, "Ctrl+R", S.SP_MediaPlay, "Start docking with the form settings")
        self.act_stop = self._action(
            "Stop", self.stop, "Ctrl+.", S.SP_MediaStop,
            "Finish the current docking step, keep the results so far and stop")

        self.act_show_log = self._action(
            "Show log", self.toggle_log, None, None, "Show or hide the log box", checkable=True)
        self.act_show_log.setChecked(True)
        self.act_reset_view = self._action(
            "Reset 3D view", self.show_selected, None, S.SP_BrowserReload,
            "Re-centre the 3D view on the selected ligand")
        self.act_about = self._action("About biosim", self.about)
        self.act_pkpd = self._action(
            "PK/PD simulator...", self.open_pkpd, "Ctrl+Shift+P", S.SP_ComputerIcon,
            "Open the pharmacokinetics / pharmacodynamics simulator")

    def _build_menus_and_toolbar(self) -> None:
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)  # keep the menus inside the window (KDE global menu)

        file_menu = menubar.addMenu("&File")
        file_menu.addAction(self.act_new)
        file_menu.addAction(self.act_open)
        self.recent_menu = file_menu.addMenu("Open &recent")
        self.recent_menu.aboutToShow.connect(self._rebuild_recent_menu)
        file_menu.addSeparator()
        file_menu.addAction(self.act_save)
        file_menu.addAction(self.act_save_as)
        file_menu.addAction(self.act_export)
        file_menu.addSeparator()
        file_menu.addAction(self.act_quit)

        edit_menu = menubar.addMenu("&Edit")
        edit_menu.addAction(self.act_edit_note)
        edit_menu.addAction(self.act_copy_smiles)
        edit_menu.addSeparator()
        edit_menu.addAction(self.act_remove)

        run_menu = menubar.addMenu("&Run")
        run_menu.addAction(self.act_run)
        run_menu.addAction(self.act_stop)

        view_menu = menubar.addMenu("&View")
        view_menu.addAction(self.act_show_log)
        view_menu.addAction(self.act_reset_view)

        tools_menu = menubar.addMenu("&Tools")
        tools_menu.addAction(self.act_pkpd)

        help_menu = menubar.addMenu("&Help")
        help_menu.addAction(self.act_about)

        bar = QToolBar("Main toolbar", self)
        bar.setMovable(False)
        bar.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        for action in (self.act_new, self.act_open, self.act_save, self.act_save_as):
            bar.addAction(action)
        bar.addSeparator()
        bar.addAction(self.act_run)
        bar.addAction(self.act_stop)
        bar.addSeparator()
        for action in (self.act_edit_note, self.act_copy_smiles, self.act_remove, self.act_export):
            bar.addAction(action)
        bar.addSeparator()
        bar.addAction(self.act_reset_view)
        bar.addSeparator()
        bar.addAction(self.act_pkpd)
        self.addToolBar(bar)

    def open_pkpd(self) -> None:
        from biosim.gui.pkpd_window import PKPDWindow

        window = getattr(self, "pkpd_window", None)
        if window is None:
            window = PKPDWindow(self)
            self.pkpd_window = window
        window.show()
        window.raise_()
        window.activateWindow()

    def update_actions(self) -> None:
        has_rows = bool(self.rows)
        selected = 0 <= self.table.currentRow() < len(self.rows)
        idle = not self.running
        self.act_new.setEnabled(idle)
        self.act_open.setEnabled(idle)
        self.recent_menu.setEnabled(idle)
        self.act_run.setEnabled(idle)
        self.run_btn.setEnabled(idle)
        can_stop = self.running and not self.stop_requested
        self.act_stop.setEnabled(can_stop)
        self.stop_btn.setEnabled(can_stop)
        self.act_save.setEnabled(idle and has_rows)
        self.act_save_as.setEnabled(idle and has_rows)
        self.act_export.setEnabled(has_rows)
        self.act_copy_smiles.setEnabled(selected)
        self.act_edit_note.setEnabled(selected and idle)
        self.act_remove.setEnabled(selected and idle)
        self.act_reset_view.setEnabled(selected)

    def update_title(self) -> None:
        if self.rows:
            self.setWindowTitle(f"biosim - {self.workdir}{' *' if self.dirty else ''}")
        else:
            self.setWindowTitle("biosim")

    def set_dirty(self, flag: bool) -> None:
        self.dirty = flag
        self.update_title()

    # ---- recent folders ---------------------------------------------------

    @staticmethod
    def _settings() -> QSettings:
        return QSettings("biosim", "biosim")

    def _recent(self) -> list[str]:
        value = self._settings().value("recent_runs", [])
        if isinstance(value, str):
            value = [value]
        return [v for v in (value or []) if isinstance(v, str)]

    def _add_recent(self, folder: Path) -> None:
        folder_s = str(Path(folder).resolve())
        items = [folder_s] + [r for r in self._recent() if r != folder_s]
        self._settings().setValue("recent_runs", items[:MAX_RECENT])

    def _rebuild_recent_menu(self) -> None:
        self.recent_menu.clear()
        items = [r for r in self._recent() if Path(r).is_dir()]
        if not items:
            self.recent_menu.addAction("(none)").setEnabled(False)
            return
        for item in items:
            action = self.recent_menu.addAction(item)
            action.triggered.connect(lambda _c=False, p=item: self.open_run(Path(p)))

    # ---- table <-> rows ---------------------------------------------------

    def refresh_table(self, select: int | None = None) -> None:
        self._populating = True
        try:
            self.table.setRowCount(len(self.rows))
            read_only = Qt.ItemIsSelectable | Qt.ItemIsEnabled
            for i, row in enumerate(self.rows):
                score = row.get("vina_score_kcal_mol")
                values = {
                    COL_RANK: str(i + 1),
                    COL_NAME: row.get("name", ""),
                    COL_SCORE: f"{score:.2f}" if score is not None else "n/a",
                    COL_ISO: str(row.get("n_isomers", 1)),
                    COL_STATUS: row.get("status", ""),
                    COL_NOTES: row.get("notes", ""),
                }
                for col, text in values.items():
                    item = QTableWidgetItem(text)
                    if col != COL_NOTES:
                        item.setFlags(read_only)
                    self.table.setItem(i, col, item)
        finally:
            self._populating = False
        if select is not None and 0 <= select < len(self.rows):
            self.table.selectRow(select)
        self.update_actions()
        self.update_title()

    def on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._populating or item.column() != COL_NOTES:
            return
        if 0 <= item.row() < len(self.rows):
            self.rows[item.row()]["notes"] = item.text()
            self.set_dirty(True)

    def on_selection_changed(self) -> None:
        self.update_actions()
        self.show_selected()

    def _selected_index(self) -> int | None:
        row = self.table.currentRow()
        return row if 0 <= row < len(self.rows) else None

    # ---- file actions -----------------------------------------------------

    def maybe_save(self) -> bool:
        """True if it is fine to discard/replace the current results."""
        if not self.dirty:
            return True
        answer = QMessageBox.question(
            self, "biosim", "You have unsaved changes to the results. Save them first?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save,
        )
        if answer == QMessageBox.Save:
            return self.save()
        return answer == QMessageBox.Discard

    def new_run(self) -> None:
        if self.running or not self.maybe_save():
            return
        self.rows = []
        self.refresh_table()
        self.viewer.clear()
        self.log.clear()
        self.workdir = Path(self.out_edit.text().strip() or "runs/gui")
        self.set_dirty(False)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.elapsed_label.setText("Elapsed 00:00")
        self.remaining_label.setText("Remaining -")
        self.statusBar().showMessage("New run: fill in the form and press Run docking")

    def open_run(self, folder: Path | None = None) -> None:
        if self.running or not self.maybe_save():
            return
        if folder is None:
            start = self.workdir.parent if self.workdir.parent.is_dir() else Path.cwd()
            chosen = QFileDialog.getExistingDirectory(self, "Open run folder", str(start))
            if not chosen:
                return
            folder = Path(chosen)
        try:
            pdb_id, rows = runio.load_run(folder)
        except (runio.RunFolderError, OSError) as exc:
            QMessageBox.warning(self, "biosim", str(exc))
            return

        self.workdir = Path(folder)
        self.pdb_id = pdb_id
        self.rows = rows
        if pdb_id:
            self.pdb_edit.setText(pdb_id)
        self.out_edit.setText(str(folder))
        self.log.clear()
        self.log.appendPlainText(f"Loaded {len(rows)} ligands from {folder / 'results.csv'}")
        if not pdb_id:
            self.log.appendPlainText("No *_protein.pdb in this folder: the 3D viewer is unavailable.")
        self.viewer.clear()
        self.progress.setRange(0, max(len(rows), 1))
        self.progress.setValue(len(rows))
        self.elapsed_label.setText("Opened from disk")
        self.remaining_label.setText("Done")
        self.set_dirty(False)
        self.refresh_table(select=0 if rows else None)
        self._add_recent(folder)
        self.statusBar().showMessage(f"Opened {folder}")

    def save(self) -> bool:
        if not self.rows:
            return True
        try:
            path = runio.save_results(self.workdir, self.rows)
        except (runio.RunFolderError, OSError) as exc:
            QMessageBox.critical(self, "biosim", f"Could not save: {exc}")
            return False
        self.set_dirty(False)
        self.statusBar().showMessage(f"Saved {path}")
        return True

    def save_as(self) -> bool:
        if not self.rows:
            return False
        suggestion = self.workdir.parent / f"{self.workdir.name}_copy"
        chosen, _ = QFileDialog.getSaveFileName(
            self, "Save run as (creates a new folder)", str(suggestion)
        )
        if not chosen:
            return False
        dest = Path(chosen)
        try:
            runio.copy_run(self.workdir, dest)
            runio.save_results(dest, self.rows)
        except (runio.RunFolderError, OSError) as exc:
            QMessageBox.critical(self, "biosim", f"Could not save as: {exc}")
            return False
        self.workdir = dest
        self.out_edit.setText(str(dest))
        self.set_dirty(False)
        self._add_recent(dest)
        self.statusBar().showMessage(f"Saved run to {dest}")
        return True

    def export_csv(self) -> None:
        if not self.rows:
            return
        chosen, _ = QFileDialog.getSaveFileName(
            self, "Export results as CSV", str(self.workdir / "results_export.csv"),
            "CSV files (*.csv)",
        )
        if not chosen:
            return
        try:
            path = runio.export_csv(Path(chosen), self.rows)
        except OSError as exc:
            QMessageBox.critical(self, "biosim", f"Could not export: {exc}")
            return
        self.statusBar().showMessage(f"Exported {path}")

    # ---- edit actions -----------------------------------------------------

    def copy_smiles(self) -> None:
        idx = self._selected_index()
        if idx is None:
            return
        row = self.rows[idx]
        smiles = row.get("best_isomer_smiles") or row.get("smiles") or ""
        QApplication.clipboard().setText(smiles)
        self.statusBar().showMessage(f"Copied SMILES of {row.get('name', '')}")

    def edit_note(self) -> None:
        idx = self._selected_index()
        if idx is None:
            return
        item = self.table.item(idx, COL_NOTES)
        if item is not None:
            self.table.setCurrentItem(item)
            self.table.editItem(item)

    def remove_selected(self) -> None:
        idx = self._selected_index()
        if idx is None:
            return
        name = self.rows[idx].get("name", "")
        answer = QMessageBox.question(
            self, "biosim",
            f"Remove '{name}' from the results table?\n"
            "Its pose file stays on disk; Save rewrites results.csv without it.",
        )
        if answer != QMessageBox.Yes:
            return
        del self.rows[idx]
        self.set_dirty(True)
        self.refresh_table(select=min(idx, len(self.rows) - 1))
        if not self.rows:
            self.viewer.clear()

    # ---- view / help ------------------------------------------------------

    def toggle_log(self, visible: bool) -> None:
        self.log.setVisible(visible)

    def about(self) -> None:
        QMessageBox.about(
            self, "About biosim",
            "biosim: an open-source platform for docking, drug effectiveness and "
            "protein folding simulations.\n\nDocking engine: AutoDock Vina. "
            "Chemistry: RDKit, Meeko, Open Babel.",
        )

    def closeEvent(self, event) -> None:
        if self.running:
            QMessageBox.information(
                self, "biosim", "Docking is still running. Wait for it to finish before closing.")
            event.ignore()
            return
        if not self.maybe_save():
            event.ignore()
            return
        event.accept()

    # ---- run control ------------------------------------------------------

    def start(self) -> None:
        if self.running:
            return
        pdb_id = self.pdb_edit.text().strip().upper()
        chembl = self.chembl_edit.text().strip()
        smiles_text = self.smiles_edit.toPlainText()
        if len(pdb_id) != 4:
            QMessageBox.warning(self, "biosim", "Enter a 4-character PDB ID.")
            return
        if not chembl and not smiles_text.strip():
            QMessageBox.warning(self, "biosim", "Enter a ChEMBL target or paste SMILES.")
            return
        if not self.maybe_save():
            return

        self.pdb_id = pdb_id
        self.workdir = Path(self.out_edit.text().strip() or "runs/gui")
        self.rows = []
        self.refresh_table()
        self.set_dirty(False)
        self.log.clear()
        self.running = True
        self.stop_requested = False
        self.update_actions()

        self.t_start = time.time()
        self.t_first_ligand = None
        self.completed = 0
        self.total = 0
        self.failed = False
        self.progress.setRange(0, 0)  # busy animation until the first ligand starts
        self.elapsed_label.setText("Elapsed 00:00")
        self.remaining_label.setText("Setting up ...")

        if not self._launch_job(chembl, smiles_text):
            self.running = False
            self.t_start = None
            self.progress.setRange(0, 1)
            self.progress.setValue(0)
            self.remaining_label.setText("Remaining -")
            self.update_actions()

    def _launch_job(self, chembl: str, smiles_text: str) -> bool:
        """Start docking in a separate process so the window never shares a GIL with Vina."""
        spec = {
            "pdb_id": self.pdb_id,
            "chembl_target": chembl,
            "smiles_text": smiles_text,
            "n_ligands": self.n_spin.value(),
            "exhaustiveness": self.exh_spin.value(),
            "workdir": str(self.workdir.resolve()),
            "enumerate_stereo": self.stereo_check.isChecked(),
            "max_isomers": self.iso_spin.value(),
        }
        try:
            self.workdir.mkdir(parents=True, exist_ok=True)
            (self.workdir / "STOP").unlink(missing_ok=True)
            job_path = self.workdir / "job.json"  # also a record of the run's settings
            job_path.write_text(json.dumps(spec, indent=2))
        except OSError as exc:
            QMessageBox.critical(self, "biosim", f"Could not prepare the run folder: {exc}")
            return False

        self._buffer = ""
        self._last_error = ""
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self.on_proc_output)
        self.proc.finished.connect(self.on_proc_finished)
        self.proc.errorOccurred.connect(self.on_proc_error)
        self.proc.start(sys.executable, ["-u", "-m", "biosim.gui.job", str(job_path)])
        return True

    def on_proc_output(self) -> None:
        if self.proc is None:
            return
        self._buffer += bytes(self.proc.readAllStandardOutput()).decode("utf-8", errors="replace")
        *lines, self._buffer = self._buffer.split("\n")
        for line in lines:
            self._handle_line(line)

    def _handle_line(self, line: str) -> None:
        line = line.rstrip("\r")
        if not line.strip() or line.startswith("WARNING: At low exhaustiveness"):
            return  # harmless Vina notice
        if line.startswith("ERROR: "):
            self._last_error = line[len("ERROR: "):]
        self.on_log(line)

    def on_proc_finished(self, exit_code: int, exit_status) -> None:
        if self._buffer.strip():
            self._handle_line(self._buffer)
        self._buffer = ""
        crashed = exit_status == QProcess.ExitStatus.CrashExit
        if exit_code == 0 and not crashed:
            try:
                _pdb, rows = runio.load_run(self.workdir)
            except (runio.RunFolderError, OSError) as exc:
                self.on_failed(f"Docking finished but the results could not be read: {exc}")
            else:
                self.on_done(rows)
        else:
            self.on_failed(
                self._last_error
                or ("The docking process crashed." if crashed
                    else f"The docking process failed (exit code {exit_code}).")
            )
        self.on_finished()
        (self.workdir / "STOP").unlink(missing_ok=True)
        if self.proc is not None:
            self.proc.deleteLater()
            self.proc = None

    def on_proc_error(self, error) -> None:
        if error == QProcess.ProcessError.FailedToStart:  # 'finished' never fires in this case
            self.on_failed("Could not start the docking process.")
            self.on_finished()

    def stop(self) -> None:
        if not self.running or self.stop_requested or self.proc is None:
            return
        self.stop_requested = True
        try:
            (self.workdir / "STOP").write_text("stop\n")  # the job polls for this file
        except OSError as exc:
            self.stop_requested = False
            QMessageBox.warning(self, "biosim", f"Could not signal the job to stop: {exc}")
            return
        self.log.appendPlainText("Stop requested: finishing the current docking step ...")
        self.remaining_label.setText("Stopping ...")
        self.update_actions()

    def on_log(self, text: str) -> None:
        self.log.appendPlainText(text)
        match = _PROGRESS_RE.search(text)
        if not match:
            return
        index, total = int(match.group(1)), int(match.group(2))
        if self.t_first_ligand is None:
            self.t_first_ligand = time.time()
        self.total = total
        self.completed = index - 1  # the log line marks the START of ligand `index`
        self.progress.setRange(0, total)
        self.progress.setValue(self.completed)

    def on_failed(self, message: str) -> None:
        self.failed = True
        self.log.appendPlainText(f"FAILED: {message}")
        QMessageBox.critical(self, "biosim", message)

    def on_finished(self) -> None:
        elapsed = time.time() - self.t_start if self.t_start else 0.0
        self.t_start = None
        self.running = False
        if self.failed:
            self.progress.setRange(0, 1)
            self.progress.setValue(0)
            self.elapsed_label.setText(f"Elapsed {fmt_time(elapsed)}")
            self.remaining_label.setText("Stopped")
        elif self.stop_requested:
            self.progress.setRange(0, max(self.total, 1))
            self.progress.setValue(len(self.rows))
            self.elapsed_label.setText(f"Stopped after {fmt_time(elapsed)}")
            self.remaining_label.setText("Stopped: partial results kept")
        else:
            self.progress.setRange(0, max(self.total, 1))
            self.progress.setValue(max(self.total, 1))
            self.elapsed_label.setText(f"Total {fmt_time(elapsed)}")
            self.remaining_label.setText("Done")
        self.update_actions()

    def on_done(self, rows: list) -> None:
        self.rows = [dict(row, notes="") for row in rows]  # run_docking already wrote results.csv
        self.set_dirty(False)
        self.refresh_table(select=0 if self.rows else None)
        self._add_recent(self.workdir)
        label = "Stopped. Partial results" if self.stop_requested else "Done. Results"
        self.log.appendPlainText(f"{label}: {self.workdir / 'results.csv'}")

    # ---- timers -----------------------------------------------------------

    def tick(self) -> None:
        self.cpu_bar.setValue(int(psutil.cpu_percent(interval=None)))
        if self.t_start is None:
            return
        now = time.time()
        self.elapsed_label.setText(f"Elapsed {fmt_time(now - self.t_start)}")
        if self.t_first_ligand is None:
            self.remaining_label.setText("Setting up ...")
        elif self.completed < 1:
            self.remaining_label.setText("Remaining: estimating ...")
        else:
            per_ligand = (now - self.t_first_ligand) / self.completed
            remaining = per_ligand * (self.total - self.completed)
            self.remaining_label.setText(f"Remaining ~{fmt_time(remaining)}")

    # ---- viewer -----------------------------------------------------------

    def show_selected(self) -> None:
        idx = self._selected_index()
        if idx is None:
            return
        name = self.rows[idx].get("name", "")
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
        pose = self.workdir / "poses" / f"{safe}.pdbqt"
        receptor = self.workdir / f"{self.pdb_id}_protein.pdb"
        if receptor.exists():
            self.viewer.show_complex(receptor, pose if pose.exists() else None)


def main() -> int:
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()
