"""Make the docking window's menus visible and add a PK/PD toolbar button.

Run once from anywhere:  python scripts/patch_app_menubar.py

1. menubar.setNativeMenuBar(False): on KDE the menu bar can be handed to a global
   menu that is not shown, which hides File/Edit/Run/View/Tools/Help entirely.
2. A "PK/PD simulator..." button at the end of the toolbar.

Safe to run twice; writes nothing unless both insertion points are found exactly
once and the patched file still compiles. Needs scripts/add_pkpd_menu.py to have
been run first.
"""
from __future__ import annotations

import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "biosim" / "gui" / "app.py"

MENUBAR_ANCHOR = "        menubar = self.menuBar()\n"
MENUBAR_NEW = (
    MENUBAR_ANCHOR
    + "        menubar.setNativeMenuBar(False)  # keep the menus inside the window (KDE global menu)\n"
)

TOOLBAR_ANCHOR = "        bar.addAction(self.act_reset_view)\n        self.addToolBar(bar)\n"
TOOLBAR_NEW = (
    "        bar.addAction(self.act_reset_view)\n"
    "        bar.addSeparator()\n"
    "        bar.addAction(self.act_pkpd)\n"
    "        self.addToolBar(bar)\n"
)


def main() -> int:
    if not APP.is_file():
        print(f"ERROR: {APP} not found")
        return 1
    text = APP.read_text()
    if "act_pkpd" not in text:
        print("ERROR: run scripts/add_pkpd_menu.py first. Nothing changed.")
        return 1
    if "setNativeMenuBar" in text:
        print("app.py already has this patch; nothing to do.")
        return 0
    for name, anchor in (("menu bar", MENUBAR_ANCHOR), ("toolbar", TOOLBAR_ANCHOR)):
        found = text.count(anchor)
        if found != 1:
            print(f"ERROR: expected the {name} insertion point once, found {found}. Nothing changed.")
            return 1
    new = text.replace(MENUBAR_ANCHOR, MENUBAR_NEW).replace(TOOLBAR_ANCHOR, TOOLBAR_NEW)
    try:
        compile(new, str(APP), "exec")
    except SyntaxError as exc:
        print(f"ERROR: the patched file would not compile ({exc}). Nothing changed.")
        return 1
    APP.write_text(new)
    print("Menu bar forced inside the window; PK/PD button added to the toolbar.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
