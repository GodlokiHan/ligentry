"""Entry point for the `ligentry-gui` command and the desktop launcher.

Works even when the conda environment is not activated (for example when the
app is started from the application menu): the environment's bin directory is
put on PATH so external tools such as Open Babel (obabel) are found.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _ensure_env_on_path() -> None:
    env_bin = str(Path(sys.executable).parent)
    parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    if env_bin not in parts:
        os.environ["PATH"] = os.pathsep.join([env_bin, *parts])


def main() -> int:
    _ensure_env_on_path()
    try:
        # lets Wayland/KDE match the window to the ligentry.desktop launcher and icon
        from PySide6.QtGui import QGuiApplication

        QGuiApplication.setDesktopFileName("ligentry")
    except Exception:
        pass
    from ligentry.gui.app import main as run_gui

    return run_gui()


if __name__ == "__main__":
    sys.exit(main())
