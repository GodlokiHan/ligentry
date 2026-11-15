"""Connects the Help window and the welcome guide to the application's windows.

* `prepare_application()` must be called before the QApplication is created.
* `add_help_menu(window, topic)` adds the Help menu (with F1 for that window's topic).
* `schedule_guide(window)` shows the welcome guide on the first launch.
* `open_help(...)` shows the single shared Help window.
"""
from __future__ import annotations

from PySide6.QtCore import QCoreApplication, QSettings, Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import QMessageBox

from ligentry.help import guide as guide_logic

_help_window = None
_guide_dialog = None


def prepare_application() -> None:
    """Needed by the embedded web browser; harmless when it is absent. Call before QApplication()."""
    try:
        QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts, True)
    except Exception:
        pass


def open_help(topic: str | None = None, query: str | None = None, source: str | None = None):
    """Show the Help window (creating it once). Optionally open a topic or run a search."""
    global _help_window
    from ligentry.gui.help_window import HelpWindow

    if _help_window is None:
        _help_window = HelpWindow()
        _help_window.destroyed.connect(_forget_help_window)
    if topic:
        _help_window.open_topic(topic)
    if query:
        _help_window.run_search(query, source)
    _help_window.show()
    _help_window.raise_()
    _help_window.activateWindow()
    return _help_window


def _forget_help_window(*_args) -> None:
    global _help_window
    _help_window = None


def close_help() -> None:
    """Close the Help window (so the application can quit when the main window closes)."""
    global _help_window
    if _help_window is not None:
        try:
            _help_window.close()
        except RuntimeError:
            pass
        _help_window = None


def show_guide(parent=None, mark_seen: bool = True):
    """Run the welcome guide as a dialog. Returns the dialog (already closed)."""
    global _guide_dialog
    from ligentry.gui.guide import GuideDialog

    dlg = GuideDialog(parent)
    _guide_dialog = dlg
    dlg.learn_more.connect(lambda topic: open_help(topic))
    dlg.exec()
    if mark_seen:
        guide_logic.mark_guide_seen(QSettings("ligentry", "ligentry"), show_again=dlg.show_again)
    _guide_dialog = None
    return dlg


def schedule_guide(window, delay_ms: int = 400) -> bool:
    """Show the guide shortly after the window appears, on the first launch only.

    Returns True when it was scheduled. LIGENTRY_NO_GUIDE=1 suppresses it.
    """
    if not guide_logic.should_show_guide(QSettings("ligentry", "ligentry")):
        return False
    QTimer.singleShot(delay_ms, lambda: show_guide(window))
    return True


def add_help_menu(window, topic: str = "welcome", about=None):
    """Add a Help menu to `window` (a QMainWindow). F1 opens the manual at `topic`."""
    menu = window.menuBar().addMenu("&Help")
    window.help_menu = menu

    def act(text, slot, shortcut=None, tip=""):
        a = QAction(text, window)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        if tip:
            a.setStatusTip(tip)
        a.triggered.connect(lambda _c=False, s=slot: s())
        return a

    window.act_help = act("ligentry &manual", lambda: open_help(topic), "F1",
                          "Open the manual at the page for this window")
    window.act_help_contents = act("Manual &contents", lambda: open_help("home"), None,
                                   "All manual pages")
    window.act_help_search = act("Search papers and the &web...",
                                 lambda: open_help(query=None).focus_search(), "Ctrl+Shift+F",
                                 "Search research papers, drug and protein databases, videos and the web")
    window.act_help_guide = act("&Welcome guide...", lambda: show_guide(window, mark_seen=False), None,
                                "A short tour of ligentry")
    menu.addAction(window.act_help)
    menu.addAction(window.act_help_contents)
    menu.addAction(window.act_help_search)
    menu.addSeparator()
    menu.addAction(window.act_help_guide)
    if about is not None:
        window.act_about = act("&About ligentry", about)
        menu.addSeparator()
        menu.addAction(window.act_about)
    return menu


def unavailable_message(parent, exc: Exception) -> None:
    QMessageBox.warning(parent, "ligentry", f"The Help window could not be opened: {exc}")
