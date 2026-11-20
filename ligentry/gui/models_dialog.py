"""In-app AI Models manager dialog (CPU default + optional GPU model pack download)."""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QDialog, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QMessageBox, QProgressBar, QPushButton, QVBoxLayout,
)

from ligentry.core.model_manager import (
    DEFAULT_GPU_PACK_URL,
    download_gpu_pack,
    get_model_status,
    install_gpu_pack_from_archive,
    remove_gpu_pack,
)


class _DownloadWorker(QThread):
    progress = Signal(int, int)
    finished_ok = Signal(str)
    failed = Signal(str)

    def __init__(self, url: str) -> None:
        super().__init__()
        self.url = url

    def run(self) -> None:
        try:
            dest = download_gpu_pack(
                url=self.url,
                progress_cb=lambda done, total: self.progress.emit(done, total),
            )
            self.finished_ok.emit(str(dest))
        except Exception as exc:
            self.failed.emit(str(exc))


class AIModelsDialog(QDialog):
    """Manage CPU default and optional downloadable GPU pose-rescoring weights."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("ligentry - AI Pose Rescorer Models")
        self.resize(560, 340)
        self._worker: _DownloadWorker | None = None

        status_box = QGroupBox("Hardware & Active Model")
        form = QFormLayout(status_box)
        self.active_label = QLabel()
        self.gpu_hw_label = QLabel()
        self.cuda_rt_label = QLabel()
        self.cpu_pkg_label = QLabel()
        self.gpu_pkg_label = QLabel()
        form.addRow("Active rescorer tier:", self.active_label)
        form.addRow("NVIDIA GPU detected:", self.gpu_hw_label)
        form.addRow("ONNX CUDA runtime:", self.cuda_rt_label)
        form.addRow("CPU default model:", self.cpu_pkg_label)
        form.addRow("GPU optional pack:", self.gpu_pkg_label)

        pack_box = QGroupBox("Optional GPU Model Pack")
        pack_lay = QVBoxLayout(pack_box)
        info = QLabel(
            "The lightweight CPU model runs out of the box on any machine. "
            "If you have an NVIDIA GPU with CUDA-enabled ONNX Runtime, you can "
            "download the higher-capacity GPU model pack below."
        )
        info.setWordWrap(True)
        pack_lay.addWidget(info)

        url_row = QHBoxLayout()
        url_row.addWidget(QLabel("Pack URL:"))
        self.url_edit = QLineEdit(DEFAULT_GPU_PACK_URL)
        url_row.addWidget(self.url_edit, 1)
        pack_lay.addLayout(url_row)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        pack_lay.addWidget(self.progress)

        btn_row = QHBoxLayout()
        self.download_btn = QPushButton("Download GPU pack")
        self.import_btn = QPushButton("Install from .tar.gz...")
        self.remove_btn = QPushButton("Remove GPU pack")
        self.download_btn.clicked.connect(self._start_download)
        self.import_btn.clicked.connect(self._install_local)
        self.remove_btn.clicked.connect(self._remove_pack)
        btn_row.addWidget(self.download_btn)
        btn_row.addWidget(self.import_btn)
        btn_row.addWidget(self.remove_btn)
        pack_lay.addLayout(btn_row)

        close_row = QHBoxLayout()
        close_row.addStretch(1)
        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.accept)
        close_row.addWidget(self.close_btn)

        lay = QVBoxLayout(self)
        lay.addWidget(status_box)
        lay.addWidget(pack_box)
        lay.addLayout(close_row)

        self.refresh_status()

    def refresh_status(self) -> None:
        st = get_model_status(prefer_gpu=True)
        self.active_label.setText(
            f"<b>{st.active_tier.upper()}</b> ({st.model_path})"
            if st.model_path
            else "Unavailable (no weights found)"
        )
        self.gpu_hw_label.setText("Yes" if st.gpu_hardware_detected else "No")
        self.cuda_rt_label.setText(
            "Available" if st.cuda_provider_available else "Not available (CPU fallback)"
        )
        cpu_ok = (st.cpu_pack_dir / "model.onnx").is_file()
        self.cpu_pkg_label.setText("Installed (bundled)" if cpu_ok else "Missing")
        self.gpu_pkg_label.setText(
            f"Installed ({st.gpu_pack_dir})" if st.gpu_pack_installed else "Not installed"
        )
        self.remove_btn.setEnabled(st.gpu_pack_installed and self._worker is None)

    def _set_busy(self, busy: bool) -> None:
        self.download_btn.setEnabled(not busy)
        self.import_btn.setEnabled(not busy)
        self.remove_btn.setEnabled(not busy)
        self.close_btn.setEnabled(not busy)
        self.progress.setVisible(busy)

    def _start_download(self) -> None:
        url = self.url_edit.text().strip()
        if not url:
            QMessageBox.warning(self, "ligentry", "Enter a valid model pack URL.")
            return
        self._set_busy(True)
        self.progress.setRange(0, 0)
        self._worker = _DownloadWorker(url)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished_ok.connect(self._on_download_ok)
        self._worker.failed.connect(self._on_download_err)
        self._worker.start()

    def _on_progress(self, done: int, total: int) -> None:
        if total > 0:
            self.progress.setRange(0, total)
            self.progress.setValue(done)
        else:
            self.progress.setRange(0, 0)

    def _on_download_ok(self, dest: str) -> None:
        self._worker = None
        self._set_busy(False)
        self.refresh_status()
        QMessageBox.information(self, "ligentry", f"GPU model pack installed to:\n{dest}")

    def _on_download_err(self, msg: str) -> None:
        self._worker = None
        self._set_busy(False)
        self.refresh_status()
        QMessageBox.critical(self, "ligentry", f"Download failed:\n{msg}")

    def _install_local(self) -> None:
        chosen, _ = QFileDialog.getOpenFileName(
            self, "Select GPU model pack (.tar.gz)", str(Path.home()),
            "Tar Gzip Archives (*.tar.gz *.tgz);;All files (*)",
        )
        if not chosen:
            return
        try:
            dest = install_gpu_pack_from_archive(chosen)
        except Exception as exc:
            QMessageBox.critical(self, "ligentry", f"Could not install pack:\n{exc}")
            return
        self.refresh_status()
        QMessageBox.information(self, "ligentry", f"GPU model pack installed to:\n{dest}")

    def _remove_pack(self) -> None:
        remove_gpu_pack()
        self.refresh_status()
