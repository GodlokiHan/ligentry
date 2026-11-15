"""The Help window: the manual (offline) and a small web browser for research.

Two tabs:

* Manual: the pages from ligentry.help, shown in a QTextBrowser, with back/forward,
  contents and search.
* Web: an embedded browser (QtWebEngine) for papers, databases and videos. If the
  installed PySide6 has no web engine, or it cannot start, addresses open in the
  user's normal browser instead and everything else keeps working.

Nothing is loaded from the web until the user searches or clicks a web link.
"""
from __future__ import annotations

from PySide6.QtCore import QSettings, Qt, QUrl
from PySide6.QtGui import QAction, QDesktopServices, QKeySequence
from PySide6.QtWidgets import (
    QComboBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox,
    QPushButton, QStyle, QTabWidget, QTextBrowser, QToolBar, QVBoxLayout, QWidget,
)

from ligentry.help import pages, search
from ligentry.help.web import DEFAULT_SOURCE, SOURCES, is_safe_url, resolve_input

MANUAL_KEY = "manual"
MANUAL_TAB, WEB_TAB = 0, 1

try:  # the web engine is optional: PySide6-Addons provides it
    from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineProfile, QWebEngineSettings
    from PySide6.QtWebEngineWidgets import QWebEngineView

    class _SafePage(QWebEnginePage):
        """Refuses anything that is not an ordinary public web page."""

        def acceptNavigationRequest(self, url, nav_type, is_main_frame):
            scheme = url.scheme().lower()
            if scheme in ("http", "https"):
                return is_safe_url(url.toString())
            if scheme in ("file", "ftp", "ftps"):
                return False
            return super().acceptNavigationRequest(url, nav_type, is_main_frame)

    class _SafeView(QWebEngineView):
        def createWindow(self, _window_type):  # links that want a new window open here instead
            return self

    WEB_ENGINE_IMPORTED = True
except Exception:  # ImportError, or a broken Qt install
    WEB_ENGINE_IMPORTED = False


class HelpWindow(QMainWindow):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("ligentry Help")
        self.resize(1000, 760)
        self.history = pages.History()
        self._hits: list = []
        self._view = None  # created on first use
        self._web_failed = not WEB_ENGINE_IMPORTED
        self._last_url = ""

        self.tabs = QTabWidget()
        self.manual = QTextBrowser()
        self.manual.setOpenLinks(False)
        self.manual.setOpenExternalLinks(False)
        self.manual.anchorClicked.connect(self._on_anchor)
        self.tabs.addTab(self.manual, "Manual")
        self.web_panel = self._build_web_panel()
        self.tabs.addTab(self.web_panel, "Web")
        self.setCentralWidget(self.tabs)
        self.tabs.currentChanged.connect(lambda _i: self._update_nav())

        self._build_actions()
        self._build_toolbar()
        self.statusBar().showMessage("Ready")
        self._restore_settings()
        self.show_page(pages.HOME_ID, record=False)

    # ---- construction -----------------------------------------------------

    def _build_web_panel(self) -> QWidget:
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(4, 4, 4, 4)
        row = QHBoxLayout()
        self.address = QLineEdit()
        self.address.setPlaceholderText("Web address, or search words")
        self.address.returnPressed.connect(self._go_address)
        self.open_external_btn = QPushButton("Open in browser")
        self.open_external_btn.setToolTip("Open this page in your normal web browser")
        self.open_external_btn.clicked.connect(self.open_in_browser)
        row.addWidget(self.address, 1)
        row.addWidget(self.open_external_btn)
        lay.addLayout(row)
        self.web_holder = QVBoxLayout()
        self.web_message = QLabel(
            "Search papers, drug and protein databases, videos and the web from the box above.\n"
            "Nothing is loaded until you search or click a web link.")
        self.web_message.setAlignment(Qt.AlignCenter)
        self.web_message.setWordWrap(True)
        self.web_holder.addWidget(self.web_message)
        lay.addLayout(self.web_holder, 1)
        return panel

    def _act(self, text, slot, shortcut=None, icon=None, tip=""):
        a = QAction(text, self)
        if icon is not None:
            a.setIcon(self.style().standardIcon(icon))
        if shortcut is not None:
            a.setShortcut(shortcut if isinstance(shortcut, QKeySequence) else QKeySequence(shortcut))
        if tip:
            a.setStatusTip(tip)
            a.setToolTip(f"{text}: {tip}")
        a.triggered.connect(lambda _c=False, s=slot: s())
        return a

    def _build_actions(self) -> None:
        S = QStyle.StandardPixmap
        self.act_back = self._act("Back", self.go_back, "Alt+Left", S.SP_ArrowBack, "Previous page")
        self.act_forward = self._act("Forward", self.go_forward, "Alt+Right", S.SP_ArrowForward, "Next page")
        self.act_home = self._act("Contents", self.go_home, "Alt+Home", S.SP_DirHomeIcon, "The manual's contents page")
        self.act_reload = self._act("Reload", self.reload_web, "F5", S.SP_BrowserReload, "Reload the web page")
        self.act_stop = self._act("Stop", self.stop_web, "Esc", S.SP_BrowserStop, "Stop loading the web page")
        self.act_focus_search = self._act("Search", self.focus_search, "Ctrl+L", None, "Jump to the search box")
        self.act_zoom_in = self._act("Larger text", lambda: self.zoom(+1), QKeySequence.ZoomIn)
        self.act_zoom_out = self._act("Smaller text", lambda: self.zoom(-1), QKeySequence.ZoomOut)
        self.act_zoom_reset = self._act("Normal size", lambda: self.zoom(0), "Ctrl+0")
        self.act_close = self._act("Close", self.close, "Ctrl+W")
        self.act_clear = self._act("Clear web data...", self.clear_web_data, None, None,
                                   "Delete cookies and cached pages of the Web tab")
        for a in (self.act_focus_search, self.act_zoom_in, self.act_zoom_out, self.act_zoom_reset, self.act_close):
            self.addAction(a)  # make the shortcuts work even though they are not in a toolbar

        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        m = menubar.addMenu("&Go")
        for a in (self.act_back, self.act_forward, self.act_home):
            m.addAction(a)
        m.addSeparator()
        m.addAction(self.act_close)
        v = menubar.addMenu("&View")
        for a in (self.act_zoom_in, self.act_zoom_out, self.act_zoom_reset):
            v.addAction(a)
        w = menubar.addMenu("&Web")
        for a in (self.act_reload, self.act_stop):
            w.addAction(a)
        w.addSeparator()
        w.addAction(self.act_clear)

    def _build_toolbar(self) -> None:
        bar = QToolBar("Help toolbar", self)
        bar.setMovable(False)
        for a in (self.act_back, self.act_forward, self.act_home, self.act_reload, self.act_stop):
            bar.addAction(a)
        bar.addSeparator()
        self.search_box = QLineEdit()
        self.search_box.setPlaceholderText("Search... (a topic, a DOI, PMID 12345678, a web address)")
        self.search_box.setClearButtonEnabled(True)
        self.search_box.returnPressed.connect(self.run_search)
        self.source_box = QComboBox()
        self.source_box.addItem("This manual", MANUAL_KEY)
        for s in SOURCES:
            self.source_box.addItem(s.name, s.key)
            self.source_box.setItemData(self.source_box.count() - 1, s.description, Qt.ToolTipRole)
        self.source_box.setToolTip("Where to search")
        go = QPushButton("Search")
        go.clicked.connect(lambda _checked=False: self.run_search())
        bar.addWidget(self.search_box)
        bar.addWidget(self.source_box)
        bar.addWidget(go)
        self.addToolBar(bar)
        self.search_box.setMinimumWidth(320)

    # ---- settings ---------------------------------------------------------

    @staticmethod
    def _settings() -> QSettings:
        return QSettings("ligentry", "ligentry")

    def _restore_settings(self) -> None:
        try:
            key = str(self._settings().value("help/source", MANUAL_KEY))
            i = self.source_box.findData(key)
            if i >= 0:
                self.source_box.setCurrentIndex(i)
            self._zoom_steps = int(self._settings().value("help/zoom", 0))
        except (TypeError, ValueError):
            self._zoom_steps = 0
        self._apply_zoom()

    def dispose_web(self) -> None:
        """Destroy the web view and its page now, before the application shuts down.

        QtWebEngine crashes at exit if a page outlives its profile, so this runs
        from help_menu.shutdown (QApplication.aboutToQuit).
        """
        view, self._view = self._view, None
        if view is None:
            return
        try:
            from shiboken6 import delete
            page = view.page()
            view.stop()
            view.setPage(None)
            delete(page)
            delete(view)
        except Exception:
            try:
                view.deleteLater()
            except Exception:
                pass
        self._web_failed = True  # this window will not create another view

    def closeEvent(self, event) -> None:
        try:
            s = self._settings()
            s.setValue("help/source", self.source_box.currentData())
            s.setValue("help/zoom", self._zoom_steps)
        except Exception:
            pass
        if self._view is not None:
            self._view.stop()
        super().closeEvent(event)

    # ---- manual -----------------------------------------------------------

    def show_page(self, page_id: str, record: bool = True) -> None:
        if record:
            self.history.visit(page_id)
        self._render_current()

    def _render_current(self) -> None:
        pid = self.history.current
        if pid.startswith(pages.SEARCH_PREFIX):
            self._hits = search.search(pid[len(pages.SEARCH_PREFIX):])
        self.manual.setHtml(pages.render(pid, self._hits))
        self.manual.verticalScrollBar().setValue(0)
        self.tabs.setCurrentIndex(MANUAL_TAB)
        self._update_nav()

    def open_topic(self, topic_id: str) -> None:
        """Show a manual page (used by F1 and the guide). Unknown ids show the contents."""
        self.show_page(topic_id if topic_id in pages.all_topic_ids() else pages.HOME_ID)

    def go_home(self) -> None:
        self.show_page(pages.HOME_ID)

    def _on_anchor(self, url: QUrl) -> None:
        link = pages.classify_link(url.toString(QUrl.FullyEncoded))
        if link.kind == "topic":
            self.show_page(link.value)
        elif link.kind == "web":
            self.open_web(link.value, link.label)
        elif link.kind == "blocked":
            self.statusBar().showMessage("That address cannot be opened here (only public web pages).", 6000)
        else:
            self.statusBar().showMessage("That link does not lead anywhere.", 4000)

    # ---- search -----------------------------------------------------------

    def focus_search(self) -> None:
        self.search_box.setFocus()
        self.search_box.selectAll()

    def run_search(self, text: str | None = None, source: str | None = None) -> None:
        text = self.search_box.text() if text is None else text
        if source is not None:
            i = self.source_box.findData(source)
            if i >= 0:
                self.source_box.setCurrentIndex(i)
        key = self.source_box.currentData()
        if text is not None and text != self.search_box.text():
            self.search_box.setText(text)
        try:
            direct = resolve_input(text, DEFAULT_SOURCE)  # DOI, PMID, address... go straight to the web
            if not direct.is_search:
                self.open_web(direct.url, direct.label)
                return
            if key == MANUAL_KEY:
                self.show_page(pages.SEARCH_PREFIX + " ".join(text.split()))
                return
            target = resolve_input(text, key)
        except ValueError as exc:
            self.statusBar().showMessage(str(exc), 6000)
            return
        self.open_web(target.url, target.label)

    # ---- web --------------------------------------------------------------

    def _ensure_view(self):
        """The embedded web view, created on first use; None when it is not available."""
        if self._view is not None:
            return self._view
        if self._web_failed:
            return None
        try:
            view = _SafeView()
            view.setPage(_SafePage(QWebEngineProfile.defaultProfile(), view))
            settings = view.settings()
            settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, False)
            settings.setAttribute(QWebEngineSettings.WebAttribute.JavascriptCanOpenWindows, False)
            try:
                settings.setAttribute(QWebEngineSettings.WebAttribute.PdfViewerEnabled, True)
            except AttributeError:
                pass
            QWebEngineProfile.defaultProfile().downloadRequested.connect(self._on_download)
            view.urlChanged.connect(self._on_url_changed)
            view.loadStarted.connect(lambda: self.statusBar().showMessage("Loading..."))
            view.loadFinished.connect(self._on_load_finished)
            view.titleChanged.connect(lambda t: self.statusBar().showMessage(t, 4000) if t else None)
        except Exception as exc:  # no GL, no sandbox, broken install...
            self._web_failed = True
            self.statusBar().showMessage(f"The embedded browser could not start ({exc}).", 8000)
            return None
        self._view = view
        self.web_message.hide()
        self.web_holder.addWidget(view, 1)
        return view

    def open_web(self, url: str, label: str = "") -> None:
        if not is_safe_url(url):
            self.statusBar().showMessage("That address cannot be opened here (only public web pages).", 6000)
            return
        self._last_url = url
        self.address.setText(url)
        view = self._ensure_view()
        if view is None:
            QDesktopServices.openUrl(QUrl(url))
            self.web_message.setText(
                "The embedded browser is not available here, so pages open in your normal browser.\n"
                "The manual works as usual. See Help > Troubleshooting for how to enable it.\n\n"
                f"Last opened: {label or url}")
            self.web_message.show()
            self.tabs.setCurrentIndex(WEB_TAB)
            self.statusBar().showMessage("Opened in your web browser.", 5000)
            self._update_nav()
            return
        self.tabs.setCurrentIndex(WEB_TAB)
        view.load(QUrl(url))
        self.statusBar().showMessage(f"Opening {label or url} ...")
        self._update_nav()

    def _go_address(self) -> None:
        text = self.address.text()
        try:
            target = resolve_input(text, self.source_box.currentData() if self.source_box.currentData() != MANUAL_KEY
                                   else DEFAULT_SOURCE)
        except ValueError as exc:
            self.statusBar().showMessage(str(exc), 6000)
            return
        self.open_web(target.url, target.label)

    def open_in_browser(self) -> None:
        url = self._view.url().toString() if self._view is not None and not self._view.url().isEmpty() else self._last_url
        if not url or not is_safe_url(url):
            self.statusBar().showMessage("There is no web page to open yet.", 4000)
            return
        QDesktopServices.openUrl(QUrl(url))

    def _on_url_changed(self, url: QUrl) -> None:
        self.address.setText(url.toString())
        self._update_nav()

    def _on_load_finished(self, ok: bool) -> None:
        self.statusBar().showMessage("Done." if ok else "The page did not load (offline, or the site blocks embedded browsers). "
                                     "Try Open in browser.", 6000)
        self._update_nav()

    def _on_download(self, download) -> None:
        name = download.downloadFileName() or "download"
        path, _ = QFileDialog.getSaveFileName(self, "Save file", name)
        if not path:
            download.cancel()
            return
        from pathlib import Path

        p = Path(path)
        download.setDownloadDirectory(str(p.parent))
        download.setDownloadFileName(p.name)
        download.accept()
        self.statusBar().showMessage(f"Saving {p.name} ...", 5000)

    def reload_web(self) -> None:
        if self._view is not None and self.tabs.currentIndex() == WEB_TAB:
            self._view.reload()

    def stop_web(self) -> None:
        if self._view is not None:
            self._view.stop()

    def clear_web_data(self) -> None:
        if not WEB_ENGINE_IMPORTED:
            QMessageBox.information(self, "ligentry", "The embedded browser is not available, so there is nothing to clear.")
            return
        ans = QMessageBox.question(self, "ligentry", "Delete the cookies and cached pages of the Web tab?")
        if ans != QMessageBox.Yes:
            return
        try:
            profile = QWebEngineProfile.defaultProfile()
            profile.clearHttpCache()
            profile.cookieStore().deleteAllCookies()
            self.statusBar().showMessage("Web data cleared.", 4000)
        except Exception as exc:
            QMessageBox.warning(self, "ligentry", f"Could not clear the web data: {exc}")

    # ---- navigation and zoom ----------------------------------------------

    def _on_web_tab(self) -> bool:
        return self.tabs.currentIndex() == WEB_TAB and self._view is not None

    def go_back(self) -> None:
        if self._on_web_tab():
            self._view.back()
        else:
            self.history.back()
            self._render_current()

    def go_forward(self) -> None:
        if self._on_web_tab():
            self._view.forward()
        else:
            self.history.forward()
            self._render_current()

    def _update_nav(self) -> None:
        if self._on_web_tab():
            hist = self._view.history()
            self.act_back.setEnabled(hist.canGoBack())
            self.act_forward.setEnabled(hist.canGoForward())
        else:
            self.act_back.setEnabled(self.history.can_back)
            self.act_forward.setEnabled(self.history.can_forward)
        web = self._on_web_tab()
        self.act_reload.setEnabled(web)
        self.act_stop.setEnabled(web)

    _zoom_steps = 0

    def zoom(self, direction: int) -> None:
        self._zoom_steps = 0 if direction == 0 else max(-4, min(8, self._zoom_steps + direction))
        self._apply_zoom()

    def _apply_zoom(self) -> None:
        font = self.manual.font()
        base = 11.0
        font.setPointSizeF(max(6.0, base + self._zoom_steps))
        self.manual.document().setDefaultFont(font)
        if self._view is not None:
            self._view.setZoomFactor(1.0 + 0.1 * self._zoom_steps)
