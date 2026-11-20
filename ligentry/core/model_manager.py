"""AI model registry, hardware detection, and optional GPU model pack downloader."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

# Default hosted release URL for the optional GPU rescorer pack (override via env var)
DEFAULT_GPU_PACK_URL = os.environ.get(
    "LIGENTRY_GPU_MODEL_URL",
    "https://github.com/GodlokiHan/ligentry/releases/download/models-v1/ligentry-rescore-gpu-v1.tar.gz",
)
# Expected SHA-256 of the downloaded archive
DEFAULT_GPU_PACK_SHA256 = os.environ.get(
    "LIGENTRY_GPU_MODEL_SHA256",
    "b797513d97d3a1555d98be7b206aec219b602a49342c105ee7028c1bc9dec365",
)


@dataclass(frozen=True)
class ModelStatus:
    active_tier: str            # "cpu" or "gpu"
    model_path: Path | None
    gpu_hardware_detected: bool
    cuda_provider_available: bool
    gpu_pack_installed: bool
    gpu_pack_dir: Path
    cpu_pack_dir: Path


def user_models_dir() -> Path:
    xdg = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg) if xdg else Path.home() / ".local" / "share"
    return base / "ligentry" / "models"


def bundled_cpu_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "models" / "cpu_v1"


def gpu_pack_dir() -> Path:
    return user_models_dir() / "gpu_v1"


def detect_nvidia_gpu() -> bool:
    """Return True if an NVIDIA GPU is present on the system."""
    if shutil.which("nvidia-smi") is not None:
        return True
    return Path("/dev/nvidia0").exists()


def available_ort_providers() -> list[str]:
    try:
        import onnxruntime as ort

        return list(ort.get_available_providers())
    except Exception:
        return []


def is_gpu_pack_installed() -> bool:
    d = gpu_pack_dir()
    model_file = d / "model.onnx"
    if not model_file.is_file():
        return False
    manifest_file = d / "manifest.json"
    if manifest_file.is_file():
        try:
            meta = json.loads(manifest_file.read_text())
            expected = meta.get("sha256", "")
            if expected:
                actual = hashlib.sha256(model_file.read_bytes()).hexdigest()
                if actual != expected:
                    return False
        except Exception:
            return False
    return True


def get_model_status(prefer_gpu: bool = True) -> ModelStatus:
    providers = available_ort_providers()
    cuda_ok = "CUDAExecutionProvider" in providers
    hw_gpu = detect_nvidia_gpu() or cuda_ok
    gpu_installed = is_gpu_pack_installed()

    cpu_dir = bundled_cpu_dir()
    g_dir = gpu_pack_dir()

    if prefer_gpu and gpu_installed and cuda_ok:
        active = "gpu"
        m_path = g_dir / "model.onnx"
    elif (cpu_dir / "model.onnx").is_file():
        active = "cpu"
        m_path = cpu_dir / "model.onnx"
    elif gpu_installed:
        active = "gpu"
        m_path = g_dir / "model.onnx"
    else:
        active = "cpu"
        m_path = None

    return ModelStatus(
        active_tier=active,
        model_path=m_path,
        gpu_hardware_detected=hw_gpu,
        cuda_provider_available=cuda_ok,
        gpu_pack_installed=gpu_installed,
        gpu_pack_dir=g_dir,
        cpu_pack_dir=cpu_dir,
    )


def install_gpu_pack_from_archive(
    archive_path: str | Path,
    expected_sha256: str = DEFAULT_GPU_PACK_SHA256,
) -> Path:
    """Verify and extract a local `.tar.gz` GPU model pack into ~/.local/share/ligentry/models/gpu_v1."""
    src = Path(archive_path)
    if not src.is_file():
        raise FileNotFoundError(f"Archive not found: {src}")

    if expected_sha256:
        digest = hashlib.sha256(src.read_bytes()).hexdigest()
        if digest.lower() != expected_sha256.lower():
            raise ValueError(f"SHA-256 mismatch: expected {expected_sha256}, got {digest}")

    target = gpu_pack_dir()
    target.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="ligentry_gpu_pack_") as tmp:
        tmp_dir = Path(tmp)
        with tarfile.open(src, "r:gz") as tf:
            # Safe extraction: disallow path traversal
            for member in tf.getmembers():
                member_path = (tmp_dir / member.name).resolve()
                if not str(member_path).startswith(str(tmp_dir.resolve())):
                    raise ValueError(f"Unsafe path in archive: {member.name}")
            tf.extractall(tmp_dir)

        found_onnx = list(tmp_dir.rglob("model.onnx"))
        if not found_onnx:
            raise ValueError("Archive does not contain model.onnx")
        src_dir = found_onnx[0].parent

        for fname in ("model.onnx", "model.onnx.data", "manifest.json"):
            f = src_dir / fname
            if f.is_file():
                shutil.copy2(f, target / fname)

    if not is_gpu_pack_installed():
        raise RuntimeError("GPU model pack failed checksum verification after extraction")
    return target


def download_gpu_pack(
    url: str = DEFAULT_GPU_PACK_URL,
    expected_sha256: str = DEFAULT_GPU_PACK_SHA256,
    progress_cb: Callable[[int, int], None] | None = None,
) -> Path:
    """Download and install the optional GPU rescorer pack with optional byte progress callback."""
    if not url.startswith(("https://", "http://")):
        raise ValueError(f"Unsupported model URL scheme: {url}")

    with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Ligentry/0.1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp, open(tmp_path, "wb") as out:
            total = int(resp.headers.get("Content-Length") or 0)
            downloaded = 0
            while True:
                chunk = resp.read(65536)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)
                if progress_cb:
                    progress_cb(downloaded, total)
        return install_gpu_pack_from_archive(tmp_path, expected_sha256=expected_sha256)
    finally:
        tmp_path.unlink(missing_ok=True)


def remove_gpu_pack() -> None:
    d = gpu_pack_dir()
    if d.exists():
        shutil.rmtree(d)
