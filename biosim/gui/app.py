"""biosim desktop application (Qt)."""
from __future__ import annotations

import re
import sys
from pathlib import Path

from biosim.gui.worker import DockWorker  # imports RDKit before Qt widgets
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QFormLayout, QLineEdit, QMainWindow,
    QMessageBox, QPlainTextEdit, QPushButton, QSpinBox, QSplitter,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)
from PySide6.QtCore import Qt

from biosim.gui.viewer import MolViewer


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("biosim")
        self.resize(1400, 850)
        self.worker: DockWorker | None = None
        self.workdir = Path("runs/gui")
        self.pdb_id = ""

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

        self.worker = DockWorker(
            pdb_id, chembl, smiles_text, self.n_spin.value(),
            self.exh_spin.value(), self.workdir,
        )
        self.worker.log.connect(self.log.appendPlainText)
        self.worker.docked.connect(self.on_done)
        self.worker.failed.connect(self.on_failed)
        self.worker.finished.connect(lambda: self.run_btn.setEnabled(True))
        self.worker.start()

    def on_failed(self, message: str) -> None:
        self.log.appendPlainText(f"FAILED: {message}")
        QMessageBox.critical(self, "biosim", message)

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
