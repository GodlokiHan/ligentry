"""Install (or remove) a desktop launcher for biosim-gui on Linux.

    python scripts/install_launcher.py
    python scripts/install_launcher.py --uninstall

Run it with the Python of the environment that has biosim installed. It writes
two files under ~/.local/share (no root needed): the application-menu entry and
the icon.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ICON_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 128 128" width="128" height="128">
  <rect width="128" height="128" rx="24" fill="#12343b"/>
  <polygon points="64,24 98.64,44 98.64,84 64,104 29.36,84 29.36,44"
           fill="none" stroke="#8fd3c7" stroke-width="5" stroke-linejoin="round"/>
  <g fill="#e8f7f4">
    <circle cx="98.64" cy="44" r="8"/>
    <circle cx="98.64" cy="84" r="8"/>
    <circle cx="64" cy="104" r="8"/>
    <circle cx="29.36" cy="84" r="8"/>
    <circle cx="29.36" cy="44" r="8"/>
  </g>
  <circle cx="64" cy="24" r="8" fill="#f2a65a"/>
</svg>
"""

DESKTOP_TEMPLATE = """[Desktop Entry]
Type=Application
Name=biosim
GenericName=Molecular docking
Comment=Open-source molecular docking with a 3D viewer
Exec="{exec_path}"
Icon=biosim
Terminal=false
Categories=Science;Chemistry;Education;
StartupNotify=true
"""


def data_home() -> Path:
    return Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")


def refresh_caches(apps_dir: Path) -> None:
    for cmd in (["update-desktop-database", str(apps_dir)], ["kbuildsycoca6"]):
        if shutil.which(cmd[0]):
            subprocess.run(cmd, capture_output=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Install the biosim desktop launcher")
    parser.add_argument("--uninstall", action="store_true", help="remove the launcher and icon")
    args = parser.parse_args()

    base = data_home()
    apps_dir = base / "applications"
    desktop_path = apps_dir / "biosim.desktop"
    icon_path = base / "icons" / "hicolor" / "scalable" / "apps" / "biosim.svg"

    if args.uninstall:
        for path in (desktop_path, icon_path):
            if path.exists():
                path.unlink()
                print(f"removed {path}")
        refresh_caches(apps_dir)
        return 0

    launcher = Path(sys.executable).parent / "biosim-gui"
    if not launcher.exists():
        print(
            f"Could not find {launcher}.\n"
            "Activate the biosim environment and run: pip install -e . --no-deps",
            file=sys.stderr,
        )
        return 1

    apps_dir.mkdir(parents=True, exist_ok=True)
    icon_path.parent.mkdir(parents=True, exist_ok=True)
    icon_path.write_text(ICON_SVG)
    desktop_path.write_text(DESKTOP_TEMPLATE.format(exec_path=launcher))
    desktop_path.chmod(0o755)
    refresh_caches(apps_dir)

    print(f"wrote {desktop_path}")
    print(f"wrote {icon_path}")
    print("Look for 'biosim' in your application menu.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
