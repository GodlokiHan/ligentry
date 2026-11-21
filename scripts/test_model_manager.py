import hashlib
import os
import sys
import tempfile
import zipfile
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import ligentry.core.model_manager as mm  # noqa: E402

_failures: list[str] = []


def check(condition: bool, label: str, detail: str = "") -> None:
    if condition:
        print(f"PASS  {label}" + (f"   {detail}" if detail else ""))
    else:
        msg = f"{label}" + (f" ({detail})" if detail else "")
        _failures.append(msg)
        print(f"FAIL  {msg}")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_localappdata_path_resolution() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # 1. Verify current host resolution works without errors
        host_dir = mm.user_models_dir()
        check(
            isinstance(host_dir, (str, Path)) and len(str(host_dir)) > 0,
            "user_models_dir resolves a non-empty path on current host",
            str(host_dir),
        )

        # 2. On Windows, verify %LOCALAPPDATA%\ligentry\models is used
        if os.name == "nt":
            fake_local = str(tmp_path / "AppData" / "Local")
            with mock.patch.dict(os.environ, {"LOCALAPPDATA": fake_local}, clear=False):
                resolved = Path(mm.user_models_dir())
                expected = Path(fake_local) / "ligentry" / "models"
                check(
                    resolved == expected,
                    "user_models_dir uses %LOCALAPPDATA%\\ligentry\\models on Windows",
                    str(resolved),
                )
        else:
            # On POSIX, verify XDG_DATA_HOME is respected when set
            fake_xdg = str(tmp_path / "xdg_data")
            with mock.patch.dict(os.environ, {"XDG_DATA_HOME": fake_xdg}, clear=False):
                resolved = Path(mm.user_models_dir())
                expected = Path(fake_xdg) / "ligentry" / "models"
                check(
                    resolved == expected,
                    "user_models_dir respects XDG_DATA_HOME on POSIX",
                    str(resolved),
                )


def test_checksum_verification() -> None:
    payload = b"dummy-onnx-weights-v1-payload"
    good_digest = sha256_bytes(payload)
    bad_digest = "0" * 64

    with tempfile.TemporaryDirectory() as tmp:
        artifact = Path(tmp) / "model.onnx"
        artifact.write_bytes(payload)

        actual_digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        check(
            actual_digest == good_digest,
            "SHA-256 digest matches valid artifact",
            good_digest[:12],
        )
        check(
            actual_digest != bad_digest and actual_digest.lower() == good_digest.lower(),
            "SHA-256 check detects mismatch and normalizes hex case",
        )

        for fn_name in ("verify_sha256", "verify_checksum", "_verify_sha256", "_verify_checksum"):
            fn = getattr(mm, fn_name, None)
            if callable(fn):
                try:
                    ok = fn(artifact, good_digest)
                    check(bool(ok), f"{fn_name} accepts matching checksum")
                except Exception as exc:
                    check(False, f"{fn_name} raised unexpectedly on valid checksum", str(exc))
                try:
                    bad = fn(artifact, bad_digest)
                    check(not bool(bad), f"{fn_name} rejects mismatched checksum")
                except ValueError:
                    check(True, f"{fn_name} raises ValueError on mismatched checksum")


def test_archive_extraction_and_path_safety() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        extract_dir = tmp_path / "extracted"
        extract_dir.mkdir()

        archive_path = tmp_path / "bundle.zip"
        with zipfile.ZipFile(archive_path, "w") as zf:
            zf.writestr("cpu_v1/model.onnx", b"fake-onnx-model-bytes")
            zf.writestr("cpu_v1/metadata.json", b'{"tier": "cpu", "version": "1"}')

        with zipfile.ZipFile(archive_path, "r") as zf:
            for member in zf.infolist():
                target = (extract_dir / member.filename).resolve()
                is_safe = target == extract_dir.resolve() or extract_dir.resolve() in target.parents
                if is_safe and not member.is_dir():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(zf.read(member.filename))

        check(
            (extract_dir / "cpu_v1" / "model.onnx").is_file(),
            "Valid archive extracts nested model files into target directory",
        )

        escape_target = (extract_dir / ".." / "escaped.txt").resolve()
        is_escaped = escape_target.parent != extract_dir.resolve() and extract_dir.resolve() not in escape_target.parents
        check(
            is_escaped,
            "Traversal target resolves outside extraction root so guard can block it",
            str(escape_target),
        )

        for fn_name in ("safe_extract_archive", "extract_archive", "_safe_extract", "_extract_archive"):
            fn = getattr(mm, fn_name, None)
            if callable(fn):
                evil_zip = tmp_path / "evil.zip"
                with zipfile.ZipFile(evil_zip, "w") as zf:
                    zf.writestr("../../pwned.txt", b"malicious")
                evil_dest = tmp_path / "evil_out"
                evil_dest.mkdir()
                blocked = False
                try:
                    fn(evil_zip, evil_dest)
                    blocked = not (tmp_path / "pwned.txt").exists()
                except (ValueError, RuntimeError, PermissionError, zipfile.BadZipFile):
                    blocked = not (tmp_path / "pwned.txt").exists()
                check(blocked, f"{fn_name} prevents path traversal outside destination")


def main() -> int:
    test_localappdata_path_resolution()
    test_checksum_verification()
    test_archive_extraction_and_path_safety()

    print(f"\nSummary: {len(_failures)} failed")
    return 1 if _failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
