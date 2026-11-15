"""Offscreen tests for the Help window, the welcome guide and the Help menus.

Run from the project folder:  python scripts/test_help_gui.py

Settings go to a temporary folder and no real browser or web page is opened. The
embedded web view itself is only exercised when you set LIGENTRY_TEST_WEB=1.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
_tmp = tempfile.mkdtemp(prefix="ligentry-help-test-")
os.environ["XDG_CONFIG_HOME"] = _tmp
os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ.setdefault("LIGENTRY_NO_GUIDE", "1")

from PySide6.QtCore import QSettings, QUrl  # noqa: E402
from PySide6.QtWidgets import QApplication, QMainWindow  # noqa: E402

from ligentry.gui import help_menu, help_window  # noqa: E402
from ligentry.gui.guide import GuideDialog  # noqa: E402
from ligentry.help import guide as guide_logic  # noqa: E402
from ligentry.help.guide import GUIDE_PAGES  # noqa: E402

help_menu.prepare_application()
app = QApplication.instance() or QApplication(sys.argv)

passed = failed = 0


def check(name: str, cond: bool, detail: str = "") -> None:
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print(f"FAIL: {name} {detail}")


opened: list[str] = []
help_window.QDesktopServices.openUrl = staticmethod(lambda url: opened.append(url.toString()) or True)

# ---------------------------------------------------------------- help window: manual
w = help_window.HelpWindow()
check("starts on the contents page", w.history.current == "home" and "ligentry manual" in w.manual.toPlainText())
check("back disabled at start", not w.act_back.isEnabled() and not w.act_forward.isEnabled())

w.open_topic("pkpd-delayed")
check("open_topic shows the page", "Delayed effect" in w.manual.toPlainText() and w.history.current == "pkpd-delayed")
check("page has see-also text", "See also" in w.manual.toPlainText())
w.open_topic("no-such-topic")
check("unknown topic falls back to contents", w.history.current == "home")

w.show_page("glossary")
w.show_page("units")
w.go_back()
check("back works", w.history.current == "glossary" and w.act_forward.isEnabled())
w.go_forward()
check("forward works", w.history.current == "units")

w._on_anchor(QUrl("topic:docking-stereo"))
check("topic link navigates", w.history.current == "docking-stereo")
n_before = len(opened)
w._on_anchor(QUrl("http://localhost:8000/x"))
check("unsafe link refused", len(opened) == n_before and "cannot be opened" in w.statusBar().currentMessage())
w._on_anchor(QUrl("topic:nope"))
check("dead link ignored", w.history.current == "docking-stereo")

w.source_box.setCurrentIndex(w.source_box.findData(help_window.MANUAL_KEY))
w.search_box.setText("stereoisomers")
w.run_search()
check("manual search shows results", w.history.current == "search:stereoisomers" and "Search results" in w.manual.toPlainText())
w.run_search("zzzqqq")
check("empty manual search says so", "Nothing in the manual" in w.manual.toPlainText())
w.run_search("   ")
check("blank search does nothing harmful", "Type something" in w.statusBar().currentMessage())

# ---------------------------------------------------------------- help window: web (fallback path, no real browser)
w._web_failed = True
w._view = None
opened.clear()
w.run_search("ke0 effect compartment", "pubmed")
check("web search without engine opens the system browser",
      len(opened) == 1 and opened[0].startswith("https://pubmed.ncbi.nlm.nih.gov/?term=ke0"), str(opened))
check("web tab shown", w.tabs.currentIndex() == help_window.WEB_TAB)
check("fallback explains itself", "normal browser" in w.web_message.text())
opened.clear()
w.run_search("PMID 12345678")
check("PMID goes straight to PubMed", opened and opened[0] == "https://pubmed.ncbi.nlm.nih.gov/12345678/", str(opened))
opened.clear()
w.run_search("10.1000/xyz123", help_window.MANUAL_KEY)
check("DOI goes to the web even when searching the manual", opened and opened[0] == "https://doi.org/10.1000/xyz123", str(opened))
opened.clear()
w.open_web("file:///etc/passwd")
check("file address refused", not opened)
w.open_web("http://192.168.1.1/")
check("private address refused", not opened)
w.address.setText("who.int")
w._go_address()
check("address bar adds https", opened and opened[0] == "https://who.int", str(opened))
opened.clear()
w.open_in_browser()
check("open in browser uses last address", opened == ["https://who.int"], str(opened))
w._on_anchor(QUrl("websearch:youtube:pharmacokinetics%20lecture"))
check("manual web-search link opens the web",
      opened and opened[-1].startswith("https://www.youtube.com/results?search_query=pharmacokinetics"), str(opened))
w.focus_search()
w.zoom(+1)
w.zoom(-1)
w.zoom(0)
check("zoom does not break", w._zoom_steps == 0)
w.close()

if os.environ.get("LIGENTRY_TEST_WEB") == "1":
    w2 = help_window.HelpWindow()
    check("web engine importable", help_window.WEB_ENGINE_IMPORTED)
    view = w2._ensure_view()
    check("web view created", view is not None, w2.statusBar().currentMessage())
    if view is not None:
        w2.open_web("https://example.org/")
        check("web view given the address", "example.org" in w2.address.text())
        check("page refuses file://", view.page().acceptNavigationRequest(QUrl("file:///etc/passwd"), 0, True) is False)
    w2.close()
else:
    print("note: embedded web view not tested (set LIGENTRY_TEST_WEB=1 to try it)")

# ---------------------------------------------------------------- guide dialog
d = GuideDialog()
check("guide starts on page 1", d.index == 0 and not d.back_btn.isEnabled() and d.next_btn.text() == "Next")
check("guide shows title", d.title_label.text() == GUIDE_PAGES[0].title)
asked: list[str] = []
d.learn_more.connect(asked.append)
d._more()
check("learn more asks for the page's topic", asked == [GUIDE_PAGES[0].topic])
for _ in range(len(GUIDE_PAGES) - 1):
    d._next()
check("guide last page says Finish", d.index == len(GUIDE_PAGES) - 1 and d.next_btn.text() == "Finish")
d.go(-5)
check("guide index clamped", d.index == 0)
check("show again unticked by default", not d.show_again)
d.again.setChecked(True)
check("show again reads checkbox", d.show_again)
d.close()

# ---------------------------------------------------------------- first-launch rule with real QSettings
s = QSettings("ligentry", "ligentry")
s.clear()
check("fresh settings: guide due", guide_logic.should_show_guide(s, {}))
guide_logic.mark_guide_seen(s)
s.sync()
check("after seen: not due", not guide_logic.should_show_guide(QSettings("ligentry", "ligentry"), {}))
guide_logic.mark_guide_seen(s, show_again=True)
check("show again: due", guide_logic.should_show_guide(QSettings("ligentry", "ligentry"), {}))
guide_logic.mark_guide_seen(s)
check("env suppression honoured by schedule_guide", help_menu.schedule_guide(QMainWindow()) is False)
s.clear()
saved = os.environ.pop("LIGENTRY_NO_GUIDE", None)
host = QMainWindow()
check("fresh settings schedule the guide", help_menu.schedule_guide(host, delay_ms=10_000_000) is True)
if saved is not None:
    os.environ["LIGENTRY_NO_GUIDE"] = saved

# ---------------------------------------------------------------- Help menu on windows
plain = QMainWindow()
help_menu.add_help_menu(plain, "cli", about=lambda: None)
check("help menu added", plain.help_menu.title() == "&Help")
check("F1 bound", plain.act_help.shortcut().toString() == "F1")
check("search shortcut bound", plain.act_help_search.shortcut().toString() == "Ctrl+Shift+F")
names = [a.text() for a in plain.help_menu.actions() if a.text()]
check("menu entries", any("manual" in n for n in names) and any("Welcome" in n for n in names) and any("About" in n for n in names), str(names))
plain.act_help.trigger()
hw = help_menu._help_window
check("F1 opens the help window at the topic", hw is not None and hw.history.current == "cli")
plain.act_help_contents.trigger()
check("contents action reuses the window", help_menu._help_window is hw and hw.history.current == "home")
plain.act_help_search.trigger()
check("search action opens the window", help_menu._help_window is hw and hw.isVisible())
help_menu.close_help()
check("close_help forgets the window", help_menu._help_window is None)

try:
    from ligentry.gui.pkpd_window import PKPDWindow

    pk = PKPDWindow()
    check("PK/PD window has Help menu", hasattr(pk, "act_help") and pk.act_help.shortcut().toString() == "F1")
    pk.act_help.trigger()
    check("PK/PD F1 opens its manual page", help_menu._help_window.history.current == "pkpd-overview")
    pk.close()
    help_menu.close_help()
except Exception as exc:  # pragma: no cover
    check("PK/PD window help", False, repr(exc))

try:
    from ligentry.gui.app import MainWindow

    mw = MainWindow()
    check("main window has Help menu", hasattr(mw, "act_help") and mw.act_help.shortcut().toString() == "F1")
    check("main window keeps About", hasattr(mw, "act_about"))
    titles = [a.text() for a in mw.menuBar().actions()]
    check("Help menu is last", titles[-1] == "&Help", str(titles))
    mw.act_help.trigger()
    check("main F1 opens docking page", help_menu._help_window.history.current == "docking-overview")
    mw.close()
    check("closing main window closes help", help_menu._help_window is None)
except Exception as exc:  # the 3D viewer can fail offscreen on some machines
    print(f"SKIP main window checks: {exc!r}")

print(f"\n{passed} passed, {failed} failed")
sys.exit(1 if failed else 0)
