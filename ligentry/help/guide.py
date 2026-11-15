"""The welcome guide: its pages and the first-launch rule (no Qt).

The guide shows once, the first time ligentry starts. Closing it records that it was
seen, unless the user ticks "Show this guide again next time". Setting the
environment variable LIGENTRY_NO_GUIDE (to anything) suppresses it, which the test
scripts use so that dialogs never block them.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Protocol

SETTINGS_KEY = "guide/seen_version"
GUIDE_VERSION = 1  # raise this when the guide changes enough to show it again
ENV_SUPPRESS = "LIGENTRY_NO_GUIDE"


@dataclass(frozen=True)
class GuidePage:
    title: str
    body: str  # simple HTML, same subset as the manual
    topic: str  # manual page opened by "Learn more"


GUIDE_PAGES: tuple[GuidePage, ...] = (
    GuidePage(
        "Welcome to ligentry",
        "<p>ligentry helps you explore computational drug design. This short guide shows what is "
        "inside and where to find help. It takes about a minute.</p>"
        "<ul><li><b>Docking</b>: predict how molecules fit a protein pocket.</li>"
        "<li><b>PK/PD simulator</b>: simulate a dose and its effect, and fit your own data.</li>"
        "<li><b>Help</b>: a manual that also searches papers, drug and protein databases and videos.</li></ul>"
        "<p>Results are for learning and exploration. They are not medical advice.</p>",
        "welcome",
    ),
    GuidePage(
        "Docking",
        "<p>In the main window give a <b>PDB ID</b> (a protein structure with a bound ligand) and either a "
        "<b>ChEMBL target</b> or your own <b>SMILES</b> list, then press <b>Run docking</b>.</p>"
        "<ul><li>Ligands are ranked by Vina score: <b>more negative is better</b>.</li>"
        "<li>Click a row to see the pose in 3D. Add your own notes to each ligand.</li>"
        "<li>Runs are saved as folders you can reopen later without docking again.</li></ul>"
        "<p>The first run needs internet to download the structure and ligands.</p>",
        "docking-overview",
    ),
    GuidePage(
        "PK/PD simulator",
        "<p>Open it from <b>Tools &gt; PK/PD simulator</b>.</p>"
        "<ul><li>Choose a dose, route and model; the concentration curve updates as you type.</li>"
        "<li>Add an <b>effect</b> curve, optionally delayed behind the concentration.</li>"
        "<li><b>Load your measured data</b> and fit the parameters; residual plots show how good the fit is.</li></ul>"
        "<p>A worked example with sample data is in the manual.</p>",
        "pkpd-overview",
    ),
    GuidePage(
        "Help and research",
        "<p>Press <b>F1</b> anywhere, or use the <b>Help</b> menu, to open the manual at the right page.</p>"
        "<ul><li>The search box searches the manual, or PubMed, Europe PMC, Google Scholar, ClinicalTrials.gov, "
        "DailyMed, PubChem, ChEMBL, UniProt, RCSB PDB, Wikipedia, YouTube and the web.</li>"
        "<li>You can type a DOI, a PubMed ID or a web address and go straight there.</li>"
        "<li>Nothing goes online unless you ask it to.</li></ul>"
        "<p>You can reopen this guide at any time from <b>Help &gt; Welcome guide</b>.</p>",
        "web-lookup",
    ),
)


class SettingsLike(Protocol):
    def value(self, key: str, default=None): ...
    def setValue(self, key: str, value) -> None: ...


def _seen_version(settings: SettingsLike) -> int:
    try:
        return int(settings.value(SETTINGS_KEY, 0))
    except (TypeError, ValueError):
        return 0


def should_show_guide(settings: SettingsLike, environ=None) -> bool:
    """True on the first launch (or after the guide was updated), unless suppressed."""
    env = os.environ if environ is None else environ
    if env.get(ENV_SUPPRESS):
        return False
    return _seen_version(settings) < GUIDE_VERSION


def mark_guide_seen(settings: SettingsLike, show_again: bool = False) -> None:
    """Record that the guide was shown; with show_again it will appear on the next start too."""
    settings.setValue(SETTINGS_KEY, 0 if show_again else GUIDE_VERSION)
