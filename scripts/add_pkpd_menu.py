"""Add Tools -> PK/PD simulator to biosim/gui/app.py (three small insertions).

Run once from anywhere:  python scripts/add_pkpd_menu.py
Safe to run twice: it does nothing if the menu entry is already there, and it
writes nothing unless all three insertion points are found exactly once and the
patched file still compiles.
"""
from __future__ import annotations

import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "biosim" / "gui" / "app.py"

ACTION_ANCHOR = '        self.act_about = self._action("About biosim", self.about)\n'
ACTION_NEW = (
    ACTION_ANCHOR
    + '        self.act_pkpd = self._action(\n'
    + '            "PK/PD simulator...", self.open_pkpd, "Ctrl+Shift+P", S.SP_ComputerIcon,\n'
    + '            "Open the pharmacokinetics / pharmacodynamics simulator")\n'
)

MENU_ANCHOR = '        help_menu = menubar.addMenu("&Help")\n'
MENU_NEW = (
    '        tools_menu = menubar.addMenu("&Tools")\n'
    + "        tools_menu.addAction(self.act_pkpd)\n\n"
    + MENU_ANCHOR
)

METHOD_ANCHOR = "    def update_actions(self) -> None:\n"
METHOD_NEW = (
    "    def open_pkpd(self) -> None:\n"
    "        from biosim.gui.pkpd_window import PKPDWindow\n\n"
    '        window = getattr(self, "pkpd_window", None)\n'
    "        if window is None:\n"
    "            window = PKPDWindow(self)\n"
    "            self.pkpd_window = window\n"
    "        window.show()\n"
    "        window.raise_()\n"
    "        window.activateWindow()\n\n"
    + METHOD_ANCHOR
)


def main() -> int:
    if not APP.is_file():
        print(f"ERROR: {APP} not found")
        return 1
    text = APP.read_text()
    if "open_pkpd" in text:
        print("app.py already has the PK/PD menu entry; nothing to do.")
        return 0
    for name, anchor in (
        ("action", ACTION_ANCHOR), ("menu", MENU_ANCHOR), ("method", METHOD_ANCHOR)
    ):
        found = text.count(anchor)
        if found != 1:
            print(f"ERROR: expected the {name} insertion point once, found {found}. Nothing changed.")
            return 1
    new = (
        text.replace(ACTION_ANCHOR, ACTION_NEW)
        .replace(MENU_ANCHOR, MENU_NEW)
        .replace(METHOD_ANCHOR, METHOD_NEW)
    )
    try:
        compile(new, str(APP), "exec")
    except SyntaxError as exc:
        print(f"ERROR: the patched file would not compile ({exc}). Nothing changed.")
        return 1
    APP.write_text(new)
    print("Added Tools -> PK/PD simulator to app.py.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
