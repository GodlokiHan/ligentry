"""biosim desktop application (Qt)."""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

from biosim.gui.worker import DockWorker  # imports RDKit before Qt widgets
import psutil
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QFormLayout, QHBoxLayout, QLabel,
    QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar,
    QPushButton, QSpinBox, QSplitter, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from biosim.gui.viewer import MolViewer

_PROGRESS_RE = re.compile(r"\[(\d+)/(\d+)\]")  # "[3/20] name" from the docking log


def fmt_time(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("biosim")
        self.resize(1400, 850)
        self.worker: DockWorker | None = None
        self.workdir = Path("runs/gui")
        self.pdb_id = ""

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
        self.out_edit = QLineEdit("runs/gui")
        self.run_btn = QPushButton("Run docking")

        form = QFormLayout()
        form.addRow("PDB ID (with bound ligand)", self.pdb_edit)
        form.addRow("ChEMBL target", self.chembl_edit)
        form.addRow("Max ligands", self.n_spin)
        form.addRow("Ligand SMILES", self.smiles_edit)
        form.addRow("Exhaustiveness", self.exh_spin)
        form.addRow("Output folder", self.out_edit)

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

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Rank", "Name", "Score (kcal/mol)", "Status"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setMaximumHeight(160)

        left = QWidget()
        lay = QVBoxLayout(left)
        lay.addLayout(form)
        lay.addWidget(self.run_btn)
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

        self.run_btn.clicked.connect(self.start)
        self.table.itemSelectionChanged.connect(self.show_selected)

        psutil.cpu_percent(interval=None)  # prime: the first call always returns 0
        self.timer = QTimer(self)
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.tick)
        self.timer.start()

    # ---- run control ----------------------------------------------------

    def start(self) -> None:
        pdb_id = self.pdb_edit.text().strip().upper()
        chembl = self.chembl_edit.text().strip()
        smiles_text = self.smiles_edit.toPlainText()
        if len(pdb_id) != 4:
            QMessageBox.warning(self, "biosim", "Enter a 4-character PDB ID.")
            return
        if not chembl and not smiles_text.strip():
            QMessageBox.warning(self, "biosim", "Enter a ChEMBL target or paste SMILES.")
            return

        self.pdb_id = pdb_id
        self.workdir = Path(self.out_edit.text().strip() or "runs/gui")
        self.table.setRowCount(0)
        self.log.clear()
        self.run_btn.setEnabled(False)

        self.t_start = time.time()
        self.t_first_ligand = None
        self.completed = 0
        self.total = 0
        self.failed = False
        self.progress.setRange(0, 0)  # busy animation until the first ligand starts
        self.elapsed_label.setText("Elapsed 00:00")
        self.remaining_label.setText("Setting up ...")

        self.worker = DockWorker(
            pdb_id, chembl, smiles_text, self.n_spin.value(),
            self.exh_spin.value(), self.workdir,
        )
        self.worker.log.connect(self.on_log)
        self.worker.docked.connect(self.on_done)
        self.worker.failed.connect(self.on_failed)
        self.worker.finished.connect(self.on_finished)
        self.worker.start()

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
        self.run_btn.setEnabled(True)
        if self.failed:
            self.progress.setRange(0, 1)
            self.progress.setValue(0)
            self.elapsed_label.setText(f"Elapsed {fmt_time(elapsed)}")
            self.remaining_label.setText("Stopped")
        else:
            self.progress.setRange(0, max(self.total, 1))
            self.progress.setValue(max(self.total, 1))
            self.elapsed_label.setText(f"Total {fmt_time(elapsed)}")
            self.remaining_label.setText("Done")

    def on_done(self, rows: list) -> None:
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            score = row["vina_score_kcal_mol"]
            values = [
                str(i + 1), row["name"],
                f"{score:.2f}" if score is not None else "n/a", row["status"],
            ]
            for col, text in enumerate(values):
                self.table.setItem(i, col, QTableWidgetItem(text))
        self.log.appendPlainText(f"Done. Results: {self.workdir / 'results.csv'}")
        if rows:
            self.table.selectRow(0)

    # ---- timers ---------------------------------------------------------

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

    # ---- viewer ---------------------------------------------------------

    def show_selected(self) -> None:
        row = self.table.currentRow()
        if row < 0 or self.table.item(row, 1) is None:
            return
        name = self.table.item(row, 1).text()
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
