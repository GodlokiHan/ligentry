"""The welcome guide dialog (shown on first launch and from Help > Welcome guide)."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QHBoxLayout, QLabel, QPushButton, QTextBrowser, QVBoxLayout,
)

from ligentry.help.guide import GUIDE_PAGES


class GuideDialog(QDialog):
    """A few pages with Back / Next. `learn_more(topic_id)` asks for a manual page."""

    learn_more = Signal(str)

    def __init__(self, parent=None, first_run: bool = False) -> None:
        super().__init__(parent)
        self.setWindowTitle("Welcome to ligentry")
        self.resize(640, 460)
        self.index = 0

        self.step_label = QLabel()
        self.step_label.setAlignment(Qt.AlignRight)
        self.title_label = QLabel()
        font = self.title_label.font()
        font.setPointSizeF(font.pointSizeF() * 1.5)
        font.setBold(True)
        self.title_label.setFont(font)
        self.body = QTextBrowser()
        self.body.setOpenLinks(False)
        self.again = QCheckBox("Show this guide again next time")
        self.again.setChecked(False)

        self.back_btn = QPushButton("Back")
        self.next_btn = QPushButton("Next")
        self.next_btn.setDefault(True)
        self.more_btn = QPushButton("Learn more in the manual")
        self.skip_btn = QPushButton("Close")
        self.back_btn.clicked.connect(lambda: self.go(self.index - 1))
        self.next_btn.clicked.connect(self._next)
        self.more_btn.clicked.connect(self._more)
        self.skip_btn.clicked.connect(self.accept)

        top = QHBoxLayout()
        top.addWidget(self.title_label, 1)
        top.addWidget(self.step_label)
        buttons = QHBoxLayout()
        buttons.addWidget(self.skip_btn)
        buttons.addWidget(self.more_btn)
        buttons.addStretch(1)
        buttons.addWidget(self.back_btn)
        buttons.addWidget(self.next_btn)
        lay = QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.body, 1)
        lay.addWidget(self.again)
        lay.addLayout(buttons)
        self.go(0)

    @property
    def show_again(self) -> bool:
        return self.again.isChecked()

    def go(self, index: int) -> None:
        self.index = max(0, min(len(GUIDE_PAGES) - 1, index))
        page = GUIDE_PAGES[self.index]
        self.title_label.setText(page.title)
        self.body.setHtml(page.body)
        self.step_label.setText(f"{self.index + 1} of {len(GUIDE_PAGES)}")
        self.back_btn.setEnabled(self.index > 0)
        last = self.index == len(GUIDE_PAGES) - 1
        self.next_btn.setText("Finish" if last else "Next")

    def _next(self) -> None:
        if self.index >= len(GUIDE_PAGES) - 1:
            self.accept()
        else:
            self.go(self.index + 1)

    def _more(self) -> None:
        self.learn_more.emit(GUIDE_PAGES[self.index].topic)
