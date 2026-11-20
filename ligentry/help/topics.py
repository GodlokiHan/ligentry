"""The ligentry manual: every help topic, as data (no Qt).

Each topic body is simple HTML limited to what Qt's rich-text widgets show:
p, ul, ol, li, b, i, code, kbd, pre, h3, h4, table (tr, th, td), br and links.
No colours or styles, so the manual follows the light or dark theme.

Links between topics are written as <a href="topic:<id>">. Links to the web are
ordinary https addresses and open in the Help window's Web tab.

scripts/test_help.py checks that every link and cross-reference resolves, that only
the allowed tags are used, and that the facts quoted here (keyboard shortcuts,
default values, CSV column names) still match the program.
"""
from __future__ import annotations

from dataclasses import dataclass

SECTIONS = (
    "Getting started",
    "Docking",
    "PK/PD simulator",
    "Reference",
    "Research and the web",
)


@dataclass(frozen=True)
class Topic:
    id: str
    title: str
    section: str
    summary: str  # one line, shown in search results and in the contents page
    body: str  # simple HTML, see the module docstring
    keywords: tuple[str, ...] = ()  # extra words that should find this topic
    related: tuple[str, ...] = ()  # ids of topics shown under "See also"
    web: tuple[str, ...] = ()  # searches offered under "Look it up on the web"


_TABLE = 'border="1" cellspacing="0" cellpadding="4"'

TOPICS: tuple[Topic, ...] = (
    # ------------------------------------------------------------------ Getting started
    Topic(
        id="welcome",
        title="Welcome to ligentry",
        section="Getting started",
        summary="What ligentry is, what is in it and where to start.",
        keywords=("overview", "about", "introduction", "what is"),
        related=("first-steps", "safety", "docking-overview", "pkpd-overview", "web-lookup"),
        body="""
<p>ligentry is an open-source desktop app for computational drug design. It has two main tools,
and a built-in manual (this window) that can also take you to papers, databases and videos on the web.</p>
<ul>
<li><b>Docking</b> predicts how a small molecule might sit in the binding pocket of a protein and gives
each molecule a score, then shows the pose in 3D. See <a href="topic:docking-overview">Docking overview</a>.</li>
<li><b>PK/PD simulator</b> simulates how a dose moves through the body (pharmacokinetics, PK) and how the
drug concentration drives an effect (pharmacodynamics, PD). It can also fit its parameters to your
measured data. See <a href="topic:pkpd-overview">PK/PD overview</a>.</li>
<li><b>Validation</b>: the docking engine comes with a redocking benchmark so you can see how well it
works, not just trust it. See <a href="topic:benchmark">Validation benchmark</a>.</li>
</ul>
<p><b>Where to start:</b> <a href="topic:first-steps">First steps</a> is a five-minute tour. If you are
new to the field, read <a href="topic:safety">What ligentry can and cannot tell you</a> first.</p>
<p>Everything in this manual works offline. Searching the web is something you choose to do:
see <a href="topic:web-lookup">Searching the web from Help</a>.</p>
""",
    ),
    Topic(
        id="first-steps",
        title="First steps: a five-minute tour",
        section="Getting started",
        summary="Dock a few ligands, then simulate and fit a dose, using the built-in examples.",
        keywords=("tutorial", "quick start", "getting started", "example", "try"),
        related=("docking-run", "pkpd-example", "shortcuts"),
        body="""
<h3>Try docking</h3>
<ol>
<li>In the main window leave the defaults: PDB ID <code>1HSG</code> (HIV protease) and ChEMBL target
<code>CHEMBL243</code>, with <b>Max ligands</b> at 5.</li>
<li>Press <b>Run docking</b> (<kbd>Ctrl+R</kbd>). The progress bar, elapsed time and log show what is happening.
It needs an internet connection the first time, to fetch the structure and the ligands.</li>
<li>When it finishes, click a row in the table to see that ligand's pose in the 3D view.</li>
<li>Type a remark in the <b>Notes</b> column (double-click it or press <kbd>F2</kbd>) and press
<kbd>Ctrl+S</kbd> to save.</li>
</ol>
<h3>Try the PK/PD simulator</h3>
<ol>
<li>Choose <b>Tools &gt; PK/PD simulator</b> (<kbd>Ctrl+Shift+P</kbd>).</li>
<li>The curve is already drawn for a single 100 mg oral dose. Change the <b>Dose</b> or the
<b>Clearance CL</b> and watch the curve update.</li>
<li>Tick <b>Effect</b> to add an effect curve, then tick <b>Delayed effect</b> and lower <b>ke0</b> to see
the effect lag behind the concentration.</li>
<li>For fitting practice follow <a href="topic:pkpd-example">Worked example: fitting data</a>.</li>
</ol>
<p>Press <kbd>F1</kbd> in any window to open the manual at the page for what you are doing.</p>
""",
    ),
    Topic(
        id="install",
        title="Starting ligentry",
        section="Getting started",
        summary="How to start the app, the launcher, and the welcome guide.",
        keywords=("launch", "start", "open", "launcher", "application menu", "install", "guide", "welcome"),
        related=("files-settings", "cli", "troubleshooting"),
        body="""
<p>ligentry runs on Linux. Windows is not supported yet.</p>
<ul>
<li>From a terminal, with the <code>ligentry</code> environment active: <code>ligentry-gui</code></li>
<li>From the application menu: run <code>python scripts/install_launcher.py</code> once, then open
<b>ligentry</b> from the menu. <code>python scripts/install_launcher.py --uninstall</code> removes it again.</li>
<li>The PK/PD simulator alone: <code>python -m ligentry.gui.pkpd_window</code></li>
</ul>
<h3>The welcome guide</h3>
<p>The first time you start ligentry, a short guide explains the main features. You can close it at any time.
Reopen it from <b>Help &gt; Welcome guide...</b>. If you tick <i>Show this guide again next time</i>
it will appear on the next start as well.</p>
<p>Installation steps are in the README on the project page.</p>
""",
    ),
    Topic(
        id="safety",
        title="What ligentry can and cannot tell you",
        section="Getting started",
        summary="Honest limits: scores are not proof, simulations are not dosing advice.",
        keywords=("disclaimer", "medical advice", "limits", "warning", "trust", "accuracy", "validity"),
        related=("docking-limits", "pkpd-limits", "web-reading"),
        body="""
<p>ligentry is a teaching and exploration tool for computational drug design. Please keep these points in mind.</p>
<ul>
<li>A docking score is a <b>hypothesis</b>, not evidence that a compound works. Docking is one early step
of drug design; real binding, safety and effect have to be tested in the laboratory.</li>
<li>The PK/PD simulator is <b>only as good as the parameters and data you give it</b>. Its output is not a
prediction for any real person and is <b>never dosing advice</b>.</li>
<li>Nothing in this manual or in ligentry is medical advice. For decisions about health or medicines,
speak to a qualified professional.</li>
<li>Pages you open from the Web tab come from other sites and are not reviewed by ligentry.
See <a href="topic:web-reading">Reading medical and scientific sources critically</a>.</li>
</ul>
""",
    ),
    # ------------------------------------------------------------------ Docking
    Topic(
        id="docking-overview",
        title="Docking overview",
        section="Docking",
        summary="What docking does, the steps ligentry takes, and what a score means.",
        keywords=("vina", "autodock", "unidock", "pose", "binding", "score", "affinity", "pocket", "protein", "ligand"),
        related=("docking-inputs", "docking-run", "docking-results", "docking-limits"),
        web=("molecular docking", "AutoDock Vina scoring function"),
        body="""
<p><b>Molecular docking</b> tries many positions and shapes (poses) of a small molecule (the ligand) inside the
binding pocket of a protein (the receptor) and estimates how well each fits. ligentry uses Uni-Dock (GPU) or AutoDock Vina (CPU)
and re-ranks poses with a neural pose rescorer.</p>
<h3>What ligentry does for each run</h3>
<ol>
<li>Downloads the protein structure you named from the RCSB PDB.</li>
<li>Finds the binding pocket from the ligand that is already bound in that crystal structure.</li>
<li>Prepares the receptor (Open Babel, protonated for pH 7.4) and each ligand (RDKit and Meeko).</li>
<li>Docks every ligand and re-ranks candidate poses with the AI neural rescorer.</li>
<li>Ranks the ligands by score and shows the table and a 3D view.</li>
</ol>
<h3>What the score means</h3>
<p>The <b>Vina score</b> is in kcal/mol. <b>More negative means a better predicted fit.</b> It is a screening
heuristic, not a measured binding affinity, and it is most useful for comparing ligands against the same
protein in the same run. See <a href="topic:docking-limits">Docking limitations</a>.</p>
""",
    ),
    Topic(
        id="docking-inputs",
        title="Choosing a target and ligands",
        section="Docking",
        summary="PDB IDs, ChEMBL targets, library files and the SMILES format.",
        keywords=("pdb id", "chembl", "smiles", "input", "ligand list", "paste", "1hsg", "chembl243", "library"),
        related=("docking-run", "docking-stereo", "glossary"),
        web=("RCSB PDB search", "ChEMBL target search", "SMILES notation tutorial"),
        body="""
<h3>PDB ID (with bound ligand)</h3>
<p>A four-character code of a structure in the Protein Data Bank, for example <code>1HSG</code>. The structure must
contain a ligand bound in the pocket you care about, because ligentry builds the search box around it.
You can look up codes at <a href="https://www.rcsb.org/">rcsb.org</a>.</p>
<h3>ChEMBL target</h3>
<p>A ChEMBL target ID such as <code>CHEMBL243</code>. ligentry fetches known active molecules for that target,
up to <b>Max ligands</b> (1 to 100,000).</p>
<h3>Your own ligands (Library file or SMILES)</h3>
<p>Load a <code>.smi</code> or <code>.csv</code> library file (up to 100,000 ligands), or leave the ChEMBL and Library fields empty and paste ligands into the <b>Ligand SMILES</b> box, one per line, as
<code>SMILES [name]</code>:</p>
<pre>CC(=O)Oc1ccccc1C(=O)O aspirin
CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine</pre>
<ul>
<li>If you give no name, the ligand is called <code>lig</code> followed by its line number.</li>
<li>Blank lines and lines starting with <code>#</code> are ignored.</li>
<li>The ChEMBL field wins over Library file and pasted SMILES.</li>
</ul>
""",
    ),
    Topic(
        id="docking-run",
        title="Running a docking job",
        section="Docking",
        summary="Every setting in the docking form and what happens when you press Run.",
        keywords=("run", "start", "exhaustiveness", "max ligands", "output folder", "progress", "form", "settings"),
        related=("docking-inputs", "docking-stereo", "docking-stop", "docking-results"),
        body="""
<table @T@>
<tr><th>Setting</th><th>Default</th><th>What it does</th></tr>
<tr><td>PDB ID (with bound ligand)</td><td>1HSG</td><td>Protein structure, and the pocket from its bound ligand.</td></tr>
<tr><td>ChEMBL target</td><td>CHEMBL243</td><td>Fetch known actives for this target. Empty = use Library file or SMILES box.</td></tr>
<tr><td>Max ligands</td><td>5</td><td>How many ligands to fetch or dock from a library (1 to 100,000).</td></tr>
<tr><td>Ligand SMILES</td><td>empty</td><td>Your own ligands, used when ChEMBL and Library file are empty.</td></tr>
<tr><td>Exhaustiveness</td><td>8</td><td>How hard the docking engine searches (1 to 64). Higher is slower and a little more thorough.</td></tr>
<tr><td>Stereoisomers</td><td>on</td><td>See <a href="topic:docking-stereo">Stereoisomers</a>.</td></tr>
<tr><td>Max isomers per ligand</td><td>8</td><td>Cap on isomers docked per ligand.</td></tr>
<tr><td>Output folder</td><td>runs/gui</td><td>Where the run is written. Reusing a folder replaces its results.</td></tr>
</table>
<p>Press <b>Run docking</b> (<kbd>Ctrl+R</kbd>). The window shows a progress bar (ligands done), elapsed
and estimated remaining time, a CPU meter and a log. The estimate appears after the first ligand finishes.</p>
<p>Docking runs in its own process at a lower priority, so the window stays responsive even while all CPU cores
are busy. You can <a href="topic:docking-stop">stop</a> a run at any time and keep what is finished.</p>
<p>Before it starts, ligentry checks that the PDB ID has four characters and that you gave either a ChEMBL target,
a library file, or some SMILES. If you have unsaved changes to earlier results it asks whether to save them first.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="docking-stereo",
        title="Stereoisomers",
        section="Docking",
        summary="How unspecified stereocentres are expanded into isomers and the best one kept.",
        keywords=("isomer", "enantiomer", "chirality", "stereocentre", "stereocenter", "cap", "enumerate"),
        related=("docking-run", "docking-results", "glossary"),
        web=("stereoisomers drug activity", "chirality in drug design"),
        body="""
<p>A SMILES string often leaves stereocentres or double bonds unspecified, but different isomers can bind very
differently. With <b>Enumerate unspecified stereoisomers</b> ticked (the default), ligentry docks every isomer
of such a ligand and keeps the best-scoring one as the ligand's result.</p>
<ul>
<li>Only the <i>unspecified</i> centres are expanded. Anything your SMILES already fixes is kept.</li>
<li><b>Max isomers per ligand</b> (default 8) caps the work. If a ligand has more isomers, a reproducible
random sample of that size is docked.</li>
<li>Salts and counter-ions are removed: only the largest fragment is docked.</li>
<li>The <b>Isomers</b> column shows how many isomers were docked, and the <b>Status</b> column notes when only
some were, for example <code>ok (3/4 isomers docked)</code>.</li>
<li>Copy SMILES (<kbd>Ctrl+Shift+C</kbd>) copies the SMILES of the best isomer.</li>
</ul>
<p>Docking time grows with the number of isomers, so untick the box (or lower the cap) for a quick look.</p>
""",
    ),
    Topic(
        id="docking-results",
        title="Reading the results table",
        section="Docking",
        summary="Columns, ranking, status messages and notes.",
        keywords=("rank", "score", "status", "notes", "table", "failed", "kcal/mol", "results", "best pose"),
        related=("docking-viewer", "docking-runs", "docking-limits", "troubleshooting"),
        body="""
<table @T@>
<tr><th>Column</th><th>Meaning</th></tr>
<tr><td>Rank</td><td>Position by AI pose score (or Vina score when rescoring is off), best first.</td></tr>
<tr><td>Name</td><td>The ligand's name from ChEMBL, library file, or from your SMILES list.</td></tr>
<tr><td>AI Score</td><td>Neural pose rescorer probability (0 to 1) that the top pose is native-like.</td></tr>
<tr><td>Pred pKd</td><td>Predicted binding affinity (-log10 Kd) from the neural rescorer.</td></tr>
<tr><td>Vina (kcal/mol)</td><td>Best Vina score. More negative is better. <code>n/a</code> means docking failed.</td></tr>
<tr><td>Isomers</td><td>How many stereoisomers were docked.</td></tr>
<tr><td>Status</td><td><code>ok</code>, or a note such as <code>ok (stopped: 2/4 isomers docked)</code> or the reason for a failure.</td></tr>
<tr><td>Notes</td><td>Your own text. Double-click or press <kbd>F2</kbd> to edit. Saved with the run.</td></tr>
</table>
<p>Scores of a few kcal/mol apart can be within the noise of the method, so do not over-read small differences.
A good score, a plausible pose and agreement with what is known about the target together are much stronger
than a score alone.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="docking-viewer",
        title="The 3D viewer",
        section="Docking",
        summary="What the 3D picture shows and how to move it.",
        keywords=("3d", "view", "pose", "rotate", "zoom", "pocket", "residues", "colours", "colors", "reset view"),
        related=("docking-results", "docking-runs"),
        body="""
<p>Select a row in the results table to see its pose.</p>
<ul>
<li>The <b>pale blue tube</b> is the protein backbone (the C-alpha trace).</li>
<li>The <b>grey sticks</b> are the pocket: residues with an atom within 5 &Aring; of the ligand.</li>
<li>The <b>ligand</b> is drawn as spheres coloured by element: carbon green, nitrogen blue, oxygen red,
sulfur yellow, phosphorus orange, fluorine and chlorine green, bromine dark red. Hydrogens are hidden.</li>
</ul>
<p>Drag with the mouse to rotate and use the wheel to zoom. If you get lost, use <b>View &gt; Reset 3D view</b>,
which re-centres on the selected ligand.</p>
<p>The viewer needs the receptor file <code>&lt;PDB&gt;_protein.pdb</code> inside the run folder.</p>
""",
    ),
    Topic(
        id="docking-runs",
        title="Saving, opening and exporting runs",
        section="Docking",
        summary="Reload a finished run without docking again; save, save as, export a CSV.",
        keywords=("open run", "save", "save as", "export", "csv", "recent", "run folder", "results.csv", "remove ligand", "notes"),
        related=("docking-results", "files-settings", "shortcuts"),
        body="""
<p>Every run is a folder. It holds <code>results.csv</code>, a <code>poses</code> folder with one
<code>.pdbqt</code> file per ligand, the receptor <code>&lt;PDB&gt;_protein.pdb</code> and
<code>job.json</code> (the settings used).</p>
<table @T@>
<tr><th>Action</th><th>Shortcut</th><th>What it does</th></tr>
<tr><td>New run</td><td>Ctrl+N</td><td>Clears the table so you can start fresh.</td></tr>
<tr><td>Open run...</td><td>Ctrl+O</td><td>Loads a finished run folder without docking again. <b>Open recent</b> lists the last eight.</td></tr>
<tr><td>Save</td><td>Ctrl+S</td><td>Writes the table, with your notes, back to the run folder.</td></tr>
<tr><td>Save as...</td><td>Ctrl+Shift+S</td><td>Copies the whole run folder to a new folder and continues there.</td></tr>
<tr><td>Export CSV...</td><td>Ctrl+E</td><td>Writes the table to any CSV file, leaving the run folder alone.</td></tr>
<tr><td>Remove ligand</td><td>Ctrl+Backspace</td><td>Removes the selected row. Its pose file stays on disk; Save rewrites results.csv without it.</td></tr>
</table>
<p>An asterisk after the folder name in the title bar means there are unsaved changes. ligentry asks before you
close or replace unsaved results.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="docking-stop",
        title="Stopping a run",
        section="Docking",
        summary="Stop keeps everything finished so far.",
        keywords=("stop", "cancel", "abort", "interrupt", "partial results", "ctrl+."),
        related=("docking-run", "troubleshooting"),
        body="""
<p>Press <b>Stop</b> (or <kbd>Ctrl+.</kbd>) to end a run early. ligentry finishes the docking step it is in, keeps
every ligand that is done, writes <code>results.csv</code> and shows the partial results. A ligand that was stopped
part-way is marked, for example <code>ok (stopped: 2/4 isomers docked)</code>.</p>
<p>The window cannot be closed while docking is running; stop the run first.</p>
""",
    ),
    Topic(
        id="docking-limits",
        title="Docking limitations",
        section="Docking",
        summary="What the docking model leaves out.",
        keywords=("limits", "limitations", "rigid receptor", "water", "protonation", "scoring", "heuristic"),
        related=("safety", "benchmark", "docking-overview"),
        web=("limitations of molecular docking", "docking scoring function accuracy"),
        body="""
<ul>
<li>Vina scores are a <b>screening heuristic</b>, not binding affinities. A good score is a hypothesis.</li>
<li>Ligands are docked in their <b>neutral form</b>; there is no protonation-state modelling for ligands.</li>
<li>Crystal <b>waters and cofactors are removed</b> from the receptor, and the receptor is <b>rigid</b>.</li>
<li>The pocket comes from the ligand bound in the PDB structure you pick. A different structure of the same
protein can give a different answer.</li>
<li>Stereoisomer enumeration multiplies the docking time for ligands with many unspecified centres.</li>
</ul>
<p>The <a href="topic:benchmark">validation benchmark</a> shows how often the engine reproduces known crystal poses
in this setup, including where it fails.</p>
""",
    ),
    Topic(
        id="benchmark",
        title="Validation benchmark",
        section="Docking",
        summary="How well redocking reproduces crystal poses, and how to run it yourself.",
        keywords=("redock", "rmsd", "validation", "baseline", "accuracy", "seeds", "benchmark", "1hsg", "4dfr"),
        related=("docking-limits", "cli", "glossary"),
        web=("redocking RMSD validation docking", "docking benchmark success rate"),
        body="""
<p>The benchmark <b>redocks</b> the crystal ligand of 8 diverse complexes with 5 random seeds each. A run passes when
the top-ranked pose is within 2 &Aring; heavy-atom RMSD of the crystal pose; a complex passes when most seeds do.</p>
<table @T@>
<tr><th>PDB</th><th>Target / ligand</th><th>Seeds passed</th></tr>
<tr><td>3PTB</td><td>Trypsin / benzamidine</td><td>5/5</td></tr>
<tr><td>1STP</td><td>Streptavidin / biotin</td><td>5/5</td></tr>
<tr><td>1M17</td><td>EGFR kinase / erlotinib</td><td>5/5</td></tr>
<tr><td>3ERT</td><td>Estrogen receptor alpha / 4-hydroxytamoxifen</td><td>5/5</td></tr>
<tr><td>1KV1</td><td>p38 MAP kinase / pyrazole urea fragment</td><td>5/5</td></tr>
<tr><td>1HSG</td><td>HIV-1 protease / indinavir</td><td>2/5</td></tr>
<tr><td>1EVE</td><td>Acetylcholinesterase / donepezil</td><td>0/5</td></tr>
<tr><td>4DFR</td><td>DHFR / methotrexate</td><td>0/5</td></tr>
</table>
<p>Overall 5 of 8 complexes pass (27 of 40 runs). For 7 of 8 complexes a near-native pose exists among the returned
poses, so most failures are <i>ranking</i> failures of the scoring function, not sampling failures.</p>
<p>Run it with <code>ligentry benchmark --seeds 42 7 123 2024 99 --out runs/benchmark</code> (about 11 minutes).
Results are repeatable for a given environment and seed, but not across RDKit versions, so compare pass rates over several
seeds rather than single runs.</p>
""".replace("@T@", _TABLE),
    ),
    # ------------------------------------------------------------------ PK/PD
    Topic(
        id="pkpd-overview",
        title="PK/PD simulator overview",
        section="PK/PD simulator",
        summary="The window, the workflow and what each part is for.",
        keywords=("pharmacokinetics", "pharmacodynamics", "pk", "pd", "simulator", "window", "tools", "layout"),
        related=("pkpd-dosing", "pkpd-model", "pkpd-effect", "pkpd-fit", "pkpd-effect-fit", "pkpd-joint-fit", "pkpd-example", "pkpd-limits"),
        web=("pharmacokinetics basics", "pharmacodynamics Emax model", "PK/PD modelling introduction"),
        body="""
<p>Open it with <b>Tools &gt; PK/PD simulator</b> (<kbd>Ctrl+Shift+P</kbd>) or run
<code>python -m ligentry.gui.pkpd_window</code>.</p>
<p><b>Pharmacokinetics (PK)</b> is what the body does to the drug: absorption, distribution, elimination. It gives
the drug concentration over time. <b>Pharmacodynamics (PD)</b> is what the drug does to the body: the effect that
goes with a concentration.</p>
<h3>The window</h3>
<ul>
<li><b>Left panel</b> (scrolls): <a href="topic:pkpd-dosing">Dosing</a>, <a href="topic:pkpd-model">Pharmacokinetics</a>,
<a href="topic:pkpd-effect">Effect</a>, <a href="topic:pkpd-readout">Time window and readout</a>,
<a href="topic:pkpd-fit">Measured data and fitting</a>, <a href="topic:pkpd-effect-fit">Effect data and fitting</a>,
<a href="topic:pkpd-joint-fit">Joint PK + PD fitting</a> and the export button.</li>
<li><b>Right</b>: the plots. Concentration on top; residual and effect panels appear when you load data or tick
Effect. The toolbar above them zooms, pans and saves a picture.</li>
</ul>
<h3>A typical session</h3>
<ol>
<li>Set the regimen and PK parameters until the curve looks right, or</li>
<li>load measured concentrations and <b>Fit parameters</b> to estimate them,</li>
<li>then add an effect and, if you have effect data, <b>Fit effect</b> (or load both datasets and use <b>Fit PK + PD jointly</b>).</li>
</ol>
<p>The curve updates a moment after every change. A red message under the form (&quot;Cannot simulate&quot;)
means a value is not allowed; the text says which.</p>
""",
    ),
    Topic(
        id="pkpd-dosing",
        title="Dosing: route, dose and schedule",
        section="PK/PD simulator",
        summary="Oral, IV bolus, and IV infusion doses, lag time Tlag, repeated dosing and interval.",
        keywords=("route", "oral", "iv", "intravenous", "bolus", "infusion", "duration", "tlag", "lag time", "dose", "interval", "number of doses", "regimen", "schedule"),
        related=("pkpd-model", "pkpd-readout", "glossary"),
        web=("oral versus intravenous bolus pharmacokinetics", "zero order intravenous infusion pharmacokinetics", "multiple dosing accumulation steady state"),
        body="""
<table @T@>
<tr><th>Setting</th><th>Default</th><th>Meaning</th></tr>
<tr><td>Route</td><td>Oral</td><td><b>Oral</b>: each dose enters the gut (times bioavailability F) after optional lag time Tlag and absorbs with rate ka.
<b>IV bolus</b>: the whole dose is injected straight into the blood at once.
<b>IV infusion</b>: each dose is infused at constant rate dose / duration over the infusion duration.</td></tr>
<tr><td>Dose</td><td>100 mg</td><td>Amount in each dose.</td></tr>
<tr><td>Number of doses</td><td>1</td><td>How many doses are given.</td></tr>
<tr><td>Interval between doses</td><td>12 h</td><td>Time between doses (used when there is more than one).</td></tr>
<tr><td>Infusion duration</td><td>1 h</td><td>Duration of each IV infusion (used only for IV infusion).</td></tr>
<tr><td>Oral lag time Tlag</td><td>0 h</td><td>Delay before oral absorption begins (used only for oral dosing).</td></tr>
</table>
<p>With an oral dose the concentration stays zero until <b>Tlag</b> and then rises smoothly. With an IV bolus it <b>jumps</b>
at the moment of the dose and then falls. With an IV infusion it climbs over the infusion duration and peaks when the
infusion stops. Overlapping repeated infusions add their infusion rates automatically.</p>
<p>The first dose starts at time 0.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="pkpd-model",
        title="The pharmacokinetic model",
        section="PK/PD simulator",
        summary="One or two compartments, absorption, clearance, Michaelis-Menten, bioavailability.",
        keywords=("compartment", "ka", "volume", "clearance", "v2", "q", "intercompartmental", "peripheral", "vmax", "km",
                  "michaelis", "bioavailability", "elimination", "absorption", "half-life"),
        related=("pkpd-dosing", "pkpd-readout", "pkpd-limits", "units", "glossary"),
        web=("one compartment model pharmacokinetics", "two compartment model pharmacokinetics",
             "Michaelis-Menten elimination saturable"),
        body="""
<p>The body is modelled as one or two well-mixed compartments. The concentration you see is C = amount in the
central compartment divided by its volume.</p>
<table @T@>
<tr><th>Parameter</th><th>Default</th><th>Meaning</th></tr>
<tr><td>Compartments</td><td>One</td><td><b>One</b>: the drug distributes instantly. <b>Two</b>: it also exchanges with a
peripheral compartment (tissues), which gives a fast early fall followed by a slower terminal phase.</td></tr>
<tr><td>Absorption rate ka</td><td>1 1/h</td><td>How fast an oral dose enters the blood. Greyed out for IV routes.</td></tr>
<tr><td>Volume V (central)</td><td>40 L</td><td>Volume of the central compartment.</td></tr>
<tr><td>Peripheral volume V2</td><td>60 L</td><td>Two compartments only.</td></tr>
<tr><td>Intercompartmental clearance Q</td><td>8 L/h</td><td>Two compartments only: how fast drug moves between the two.</td></tr>
<tr><td>Elimination</td><td>Linear</td><td><b>Linear</b>: removal proportional to concentration, set by clearance CL.
<b>Michaelis-Menten</b>: removal saturates at high concentration (Vmax in mg/h, Km in mg/L).</td></tr>
<tr><td>Clearance CL</td><td>5 L/h</td><td>Volume of blood cleared of drug per hour (linear elimination).</td></tr>
<tr><td>Bioavailability F</td><td>1.0</td><td>Fraction of an oral dose that reaches the gut compartment. Greyed out for IV routes.</td></tr>
</table>
<p>Boxes that do not apply to your choices are greyed out. Half-lives are shown in the readout; for one compartment
with linear elimination the half-life is ln 2 &times; V / CL. With Michaelis-Menten elimination the half-life shown describes the
low-concentration tail.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="pkpd-effect",
        title="The effect model (Emax)",
        section="PK/PD simulator",
        summary="Baseline, maximum effect, EC50 and the Hill coefficient.",
        keywords=("effect", "emax", "ec50", "hill", "e0", "baseline", "inhibition", "sigmoid", "dose response", "concentration-effect"),
        related=("pkpd-delayed", "pkpd-effect-fit", "glossary"),
        web=("Emax model pharmacodynamics", "EC50 Hill coefficient concentration-effect"),
        body="""
<p>Tick the <b>Effect</b> box to add an effect curve. The effect follows the sigmoid Emax relationship:</p>
<p><code>E = E0 + Emax &times; C<sup>h</sup> / (EC50<sup>h</sup> + C<sup>h</sup>)</code></p>
<table @T@>
<tr><th>Parameter</th><th>Default</th><th>Meaning</th></tr>
<tr><td>Baseline E0</td><td>0</td><td>Effect with no drug.</td></tr>
<tr><td>Emax</td><td>100</td><td>Largest change from baseline. <b>Negative for an inhibitory effect.</b></td></tr>
<tr><td>EC50</td><td>1 mg/L</td><td>Concentration giving half of Emax.</td></tr>
<tr><td>Hill coefficient h</td><td>1</td><td>Steepness of the curve. Larger is steeper.</td></tr>
</table>
<p>The effect unit is whatever you measure (a percentage, a blood pressure change, ...). ligentry only uses it as a number.</p>
<p>By default the effect follows the plasma concentration with no delay. If the real effect lags behind the concentration,
use <a href="topic:pkpd-delayed">Delayed effect</a>.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="pkpd-delayed",
        title="Delayed effect (effect compartment)",
        section="PK/PD simulator",
        summary="ke0, the effect site and hysteresis.",
        keywords=("ke0", "delayed", "effect compartment", "effect site", "hysteresis", "lag", "equilibration", "delay"),
        related=("pkpd-effect", "pkpd-effect-fit", "pkpd-residuals", "pkpd-limits"),
        web=("effect compartment model ke0", "hysteresis loop pharmacodynamics", "biophase distribution delay"),
        body="""
<p>Some drugs act at a site that the blood concentration reaches with a delay, so the effect peaks <i>after</i> the
plasma concentration. Tick <b>Delayed effect (effect compartment)</b> in the Effect box and set <b>ke0</b>, the rate
(1/h) at which an imaginary effect-site concentration Ce chases the plasma concentration C:</p>
<p><code>dCe/dt = ke0 &times; (C &minus; Ce)</code></p>
<p>The effect is then computed from Ce instead of C, with the same Emax formula.</p>
<ul>
<li>The equilibration half-life is ln 2 / ke0 and is shown in the readout. Small ke0 means a long delay.</li>
<li>A very large ke0 gives back the direct effect.</li>
<li>The concentration plot shows Ce as a <b>dashed green &quot;effect site&quot; curve</b>.</li>
<li>The effect compartment does not change the PK: it holds a negligible amount of drug.</li>
<li>Plotting effect against concentration for a delayed effect makes a loop (<b>hysteresis</b>): the same
concentration gives different effects on the way up and on the way down.</li>
</ul>
<p>To estimate ke0 from data, tick <b>ke0</b> beside the other effect parameters in
<a href="topic:pkpd-effect-fit">Effect data and fitting</a>.</p>
""",
    ),
    Topic(
        id="pkpd-readout",
        title="Time window and readout",
        section="PK/PD simulator",
        summary="Cmax, Tmax, AUC, half-lives and how the plotted time range is chosen.",
        keywords=("cmax", "tmax", "auc", "half-life", "time window", "end time", "readout", "metrics", "rmse"),
        related=("pkpd-model", "glossary", "units"),
        web=("AUC Cmax Tmax pharmacokinetic parameters", "elimination half-life calculation"),
        body="""
<p>The text under the form summarises the simulated curve.</p>
<ul>
<li><b>Cmax</b> and <b>Tmax</b>: the highest concentration and when it occurs (over the simulated points).</li>
<li><b>AUC(0-t)</b>: area under the concentration curve over the plotted window (trapezoid rule), in mg&middot;h/L.
It measures total exposure.</li>
<li><b>Half-lives</b>: terminal half-life, plus the distribution half-life for two compartments.</li>
<li><b>Effect-site equilibration half-life</b> when the delayed effect is on.</li>
<li><b>RMSE</b> against your data when concentrations or effects are loaded.</li>
</ul>
<h3>Time window</h3>
<p>By default the window is chosen automatically: long enough to show the last dose being absorbed or infused and about five
terminal half-lives of elimination, and extended to cover your data. Untick <b>Choose the time window automatically</b>
to type your own <b>End time</b> (up to 100 000 h).</p>
""",
    ),
    Topic(
        id="pkpd-data-format",
        title="Data file format",
        section="PK/PD simulator",
        summary="How to lay out concentration and effect CSV files.",
        keywords=("csv", "file", "format", "load data", "header", "columns", "delimiter", "decimal comma", "tsv", "import"),
        related=("pkpd-fit", "pkpd-effect-fit", "pkpd-example", "troubleshooting"),
        body="""
<p>Both loaders read plain text with one observation per line: <b>time in hours</b> and then either the
<b>concentration in mg/L</b> or the <b>effect</b> in your own unit.</p>
<pre>time,conc
0.5,1.8
1,3.1</pre>
<ul>
<li>A <b>header line</b> is optional. When there is one, columns are found by name, in any order, and extra columns are ignored.
Time columns may be called time, t, h, hr, hours. Concentration columns: conc..., c, cp. Effect columns:
effect..., response..., e, pd, y. A header with two unrecognised names is read as time, then value.</li>
<li>Without a header the first column is time and the second is the value. An effect file without a header must have exactly two columns.</li>
<li>Separators: comma, semicolon, tab or spaces. With semicolons or tabs a comma is a decimal separator (1,5 means 1.5).</li>
<li>Blank lines and lines starting with <code>#</code> are skipped. Rows are sorted by time.</li>
<li>Times and concentrations must not be negative. Effects may be negative.</li>
</ul>
<p>If a line is wrong, the error message gives its line number.</p>
""",
    ),
    Topic(
        id="pkpd-fit",
        title="Fitting PK parameters to concentrations",
        section="PK/PD simulator",
        summary="Load data, choose parameters, fit, read the table and notes, use the result.",
        keywords=("fit", "fitting", "estimate", "cv", "weighting", "proportional", "rmse", "optimiser", "load data", "use fitted values",
                  "uncertainty", "least squares", "tlag"),
        related=("pkpd-data-format", "pkpd-effect-fit", "pkpd-joint-fit", "pkpd-residuals", "pkpd-example", "pkpd-limits"),
        web=("nonlinear least squares pharmacokinetic parameter estimation", "population versus individual PK fitting"),
        body="""
<ol>
<li>Set the <b>Dosing</b> and the model structure (route, compartments, elimination) to match how your data were
collected, and enter starting values that are roughly right.</li>
<li>Press <b>Load data...</b> and choose your concentration file (<a href="topic:pkpd-data-format">format</a>).
The points appear on the plot.</li>
<li>Tick the parameters to fit (including <b>Tlag</b> for oral data if there is a delay before absorption starts).
Boxes that do not apply are greyed out (for example ka and Tlag for IV bolus or IV infusion).</li>
<li>Choose the <b>Weighting</b>: <i>Proportional</i> (relative error, the usual choice, so small concentrations
count as much as large ones) or <i>None</i> (absolute error).</li>
<li>Press <b>Fit parameters</b>. The table shows each estimate, its unit and its <b>CV %</b>; below it are the RMSE and notes.</li>
<li>Press <b>Use fitted values</b> to copy the estimates into the form.</li>
</ol>
<h3>Reading the result</h3>
<ul>
<li><b>CV %</b> is the standard error as a percentage of the estimate. Above 50 % the parameter is poorly determined.</li>
<li>The notes warn about common traps: <b>flip-flop</b> between ka and the elimination rate (oral data can swap them),
apparent values (<b>V/F</b>, <b>CL/F</b>) after oral dosing, two-compartment fits that need early samples, and phases
that are too close together to tell apart.</li>
<li>Always look at the <a href="topic:pkpd-residuals">residual plot</a> as well as the numbers.</li>
</ul>
<p>The fit is a local search, so it can land in a wrong answer if the starting values are far off. Try different starts if a fit looks odd.</p>
""",
    ),
    Topic(
        id="pkpd-effect-fit",
        title="Fitting the effect model",
        section="PK/PD simulator",
        summary="Load effect data, choose parameters (including ke0), fit and read the result.",
        keywords=("effect fit", "fit effect", "ec50", "emax", "e0", "hill", "ke0", "start from the values above", "standard error", "se"),
        related=("pkpd-effect", "pkpd-delayed", "pkpd-joint-fit", "pkpd-residuals", "pkpd-data-format", "pkpd-example"),
        web=("fitting Emax model EC50 estimation", "sparse data pharmacodynamic parameter identifiability"),
        body="""
<p>The effect depends on the concentration, and ligentry takes the concentration at each effect time from the
<b>PK model currently in the form</b>. So <b>fit or enter the PK parameters first</b>, then fit the effect.</p>
<ol>
<li>Press <b>Load effect data...</b> and choose a file with time and effect (<a href="topic:pkpd-data-format">format</a>).</li>
<li>Tick the effect parameters to fit: <b>E0, Emax, EC50, Hill</b>, and <b>ke0</b> if the effect may lag behind the concentration.
Parameters you leave unticked keep the values in the Effect box. (Hill is off by default because it needs data across
the whole curve.)</li>
<li>Press <b>Fit effect</b>. The table shows each estimate with its standard error (&plusmn; SE).</li>
<li>Press <b>Use fitted effect values</b> to copy the estimates into the form. If ke0 was fitted this also switches the delayed effect on.</li>
</ol>
<p>By default the fitted parameters start from values guessed from your data. Tick <b>Start from the values above</b> to start
from the form instead. For ke0, the default is to try a range of values first and begin from the best.</p>
<h3>Warnings you may see</h3>
<ul>
<li><b>EC50 above the highest concentration</b>: the data never approach saturation, so Emax and EC50 trade off against each other and are weakly determined.</li>
<li><b>A parameter on its limit</b>: the data cannot determine it. For ke0 on its upper limit this means <b>no detectable delay</b>: a direct effect fits as well.</li>
<li><b>All concentrations are zero</b>: your effect times are before the drug is present. Check the dosing and the times.</li>
<li>With a delayed effect, ke0 and EC50 are only as reliable as the PK model: a wrong PK model can look like a delay.</li>
</ul>
""",
    ),
    Topic(
        id="pkpd-joint-fit",
        title="Joint PK + PD fitting",
        section="PK/PD simulator",
        summary="Simultaneously fit PK and PD parameters against concentration and effect data.",
        keywords=("joint", "simultaneous", "pk+pd", "joint fit", "fit pk + pd jointly", "use joint fitted values", "uncertainty propagation"),
        related=("pkpd-fit", "pkpd-effect-fit", "pkpd-delayed", "pkpd-residuals", "pkpd-example"),
        web=("simultaneous PK PD modelling", "sequential versus simultaneous PK/PD estimation"),
        body="""
<p>Sequential fitting estimates PK parameters first and then treats those predicted concentrations as fixed, noiseless
inputs when fitting the effect. <b>Joint PK + PD fitting</b> optimizes all ticked PK and PD parameters simultaneously
against both the measured concentrations and the measured effects.</p>
<h3>When to use joint fitting</h3>
<ul>
<li>When effect measurements are sampled densely around the peak or hysteresis loop and can help pin down absorption
or elimination rates alongside plasma samples.</li>
<li>When you want parameter uncertainties (CV % for PK, &plusmn; SE for PD) computed from the full joint covariance matrix
so PK uncertainty propagates into EC50 and ke0.</li>
</ul>
<h3>How to run a joint fit</h3>
<ol>
<li>Load both concentration data (<b>Load data...</b>) and effect data (<b>Load effect data...</b>). The
<b>Fit PK + PD jointly</b> button becomes active once both datasets are loaded.</li>
<li>Tick the PK parameters to fit in the PK fitting box and the PD parameters (including <b>ke0</b> for a delayed effect)
in the Effect fitting box.</li>
<li>Press <b>Fit PK + PD jointly</b>. ligentry balances the concentration and effect residual scales automatically,
runs a fast warm-start, and polishes all free PK and PD parameters together.</li>
<li>Press <b>Use joint fitted values</b> to copy both the PK and PD estimates back into the form.</li>
</ol>
""",
    ),
    Topic(
        id="pkpd-residuals",
        title="Reading the residual plots",
        section="PK/PD simulator",
        summary="What good and bad residuals look like and what to do about them.",
        keywords=("residual", "residuals", "pattern", "goodness of fit", "diagnostic", "systematic", "scatter"),
        related=("pkpd-fit", "pkpd-effect-fit", "pkpd-delayed"),
        web=("residual plot interpretation model fit", "goodness of fit plots pharmacometrics"),
        body="""
<p>A <b>residual</b> is observed minus predicted. The small panels under the concentration and effect plots show one
stem per data point.</p>
<ul>
<li><b>Good:</b> points scatter randomly around zero, with no pattern.</li>
<li><b>Bad:</b> a smooth curve, runs of the same sign, or residuals that grow with concentration.</li>
</ul>
<table @T@>
<tr><th>Pattern</th><th>Possible cause</th></tr>
<tr><td>Concentration residuals positive early, negative later (or the reverse)</td><td>Wrong number of compartments, or wrong route or dose entered.</td></tr>
<tr><td>Effect residuals positive on the way up and negative on the way down</td><td>The effect lags the concentration: try the <a href="topic:pkpd-delayed">delayed effect</a>.</td></tr>
<tr><td>Spread grows with the value</td><td>Use proportional weighting.</td></tr>
</table>
<p>A low RMSE with a clear pattern still means the model is wrong.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="pkpd-example",
        title="Worked example: fitting data",
        section="PK/PD simulator",
        summary="Made-up concentration and effect data you can save and fit yourself.",
        keywords=("example", "tutorial", "sample data", "practice", "walkthrough", "try it"),
        related=("pkpd-fit", "pkpd-effect-fit", "pkpd-data-format", "first-steps"),
        body="""
<p>These made-up data come from a one-compartment oral model (dose 100 mg, ka 1.2 1/h, V 30 L, CL 5 L/h), rounded to
three digits. Save each block as a text file (for example <code>conc.csv</code> and <code>effect.csv</code>).</p>
<h3>Concentrations</h3>
<pre>time,conc
0.25,0.845
0.5,1.44
1,2.11
2,2.42
3,2.24
4,1.96
6,1.42
8,1.02
12,0.524
24,0.0709</pre>
<h3>Effects (at the same times)</h3>
<pre>time,effect
0.25,45.8
0.5,59.0
1,67.8
2,70.8
3,69.1
4,66.2
6,58.7
8,50.5
12,34.4
24,6.6</pre>
<h3>Steps</h3>
<ol>
<li>Open the PK/PD simulator. Keep <b>Oral</b>, <b>One compartment</b> and <b>Linear</b>. Set <b>Dose</b> to 100 mg.</li>
<li>Press <b>Load data...</b> and choose <code>conc.csv</code>. Keep ka, V and CL ticked and press <b>Fit parameters</b>.
You should get values close to ka 1.2, V 30 and CL 5 (they will differ slightly because the data are rounded).
Press <b>Use fitted values</b>.</li>
<li>Tick <b>Effect</b>. Press <b>Load effect data...</b>, choose <code>effect.csv</code>, keep E0, Emax and EC50 ticked and press
<b>Fit effect</b>. Expect E0 near 0, Emax near 100 and EC50 near 1.</li>
<li>Now tick <b>ke0</b> as well and fit again. These data have no delay, so ke0 finishes on its upper limit and the notes say
there is no detectable delay. That is the correct answer.</li>
</ol>
<p>Then try changing a data value by hand to see how the residuals and CV % react.</p>
""",
    ),
    Topic(
        id="pkpd-export",
        title="Exporting a simulation",
        section="PK/PD simulator",
        summary="The CSV written by Export simulation CSV and what its columns mean.",
        keywords=("export", "csv", "save", "columns", "simulation output"),
        related=("pkpd-readout", "units"),
        body="""
<p>Press <b>Export simulation CSV...</b> to save the plotted curves. The columns, in order:</p>
<table @T@>
<tr><th>Column</th><th>Meaning</th></tr>
<tr><td>time_h</td><td>Time in hours.</td></tr>
<tr><td>conc_mg_per_L</td><td>Plasma concentration.</td></tr>
<tr><td>gut_mg</td><td>Amount left in the gut (zero for IV routes).</td></tr>
<tr><td>central_mg</td><td>Amount in the central compartment.</td></tr>
<tr><td>peripheral_mg</td><td>Amount in the peripheral compartment (only with two compartments).</td></tr>
<tr><td>effect_site_mg_per_L</td><td>Effect-site concentration (only with a delayed effect).</td></tr>
<tr><td>effect</td><td>The effect (only when the Effect box is ticked).</td></tr>
</table>
<p>To save a picture of the plots, use the disk icon in the plot toolbar.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="pkpd-limits",
        title="PK/PD limitations",
        section="PK/PD simulator",
        summary="What the simulator does not model and what fits cannot tell apart.",
        keywords=("limits", "limitations", "assumptions", "flip-flop", "identifiability", "apparent volume", "not modelled"),
        related=("safety", "pkpd-fit", "pkpd-effect-fit"),
        web=("pharmacokinetic model assumptions limitations", "flip-flop kinetics"),
        body="""
<ul>
<li>It is a modelling tool, not a prediction of how a drug behaves in people, and nothing here is dosing advice.</li>
<li>At most <b>two compartments</b>; <b>oral, IV bolus, or IV infusion</b> dosing (no other extravascular routes).</li>
<li>The effect is <b>direct or one effect compartment</b>. There are no indirect-response, tolerance or disease-progression models.</li>
<li>From <b>oral data alone, F and V cannot be separated</b>: fitted volumes are V/F and V2/F, and clearances are CL/F and Q/F.</li>
<li><b>Flip-flop</b>: with oral data, absorption and elimination can swap and fit equally well.</li>
<li>A two-compartment fit needs samples from the early, rapidly falling phase. Without them V2 and Q are poorly determined.</li>
<li>EC50 and Emax trade off when the data never approach EC50.</li>
<li>Fits describe <b>one</b> data set. There is no population (many-subject) modelling.</li>
</ul>
""",
    ),
    # ------------------------------------------------------------------ Reference
    Topic(
        id="shortcuts",
        title="Keyboard shortcuts",
        section="Reference",
        summary="All shortcuts in the main window, the PK/PD window and the Help window.",
        keywords=("shortcut", "keys", "keyboard", "hotkey", "ctrl", "f1"),
        related=("docking-runs", "docking-run"),
        body="""
<h3>Main window</h3>
<table @T@>
<tr><th>Keys</th><th>Action</th></tr>
<tr><td>Ctrl+N</td><td>New run</td></tr>
<tr><td>Ctrl+O</td><td>Open run...</td></tr>
<tr><td>Ctrl+S</td><td>Save</td></tr>
<tr><td>Ctrl+Shift+S</td><td>Save as...</td></tr>
<tr><td>Ctrl+E</td><td>Export CSV...</td></tr>
<tr><td>Ctrl+Q</td><td>Quit</td></tr>
<tr><td>F2</td><td>Edit note of the selected ligand</td></tr>
<tr><td>Ctrl+Shift+C</td><td>Copy SMILES</td></tr>
<tr><td>Ctrl+Backspace</td><td>Remove ligand</td></tr>
<tr><td>Ctrl+R</td><td>Run docking</td></tr>
<tr><td>Ctrl+.</td><td>Stop</td></tr>
<tr><td>Ctrl+Shift+P</td><td>PK/PD simulator</td></tr>
<tr><td>F1</td><td>Open this manual</td></tr>
<tr><td>Ctrl+Shift+F</td><td>Search papers and the web</td></tr>
</table>
<h3>Help window</h3>
<table @T@>
<tr><th>Keys</th><th>Action</th></tr>
<tr><td>Ctrl+L</td><td>Jump to the search box</td></tr>
<tr><td>Alt+Left, Alt+Right</td><td>Back, Forward</td></tr>
<tr><td>Alt+Home</td><td>Manual home page</td></tr>
<tr><td>Ctrl++, Ctrl+-, Ctrl+0</td><td>Larger text, smaller text, normal size</td></tr>
<tr><td>Ctrl+W</td><td>Close the Help window</td></tr>
</table>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="glossary",
        title="Glossary",
        section="Reference",
        summary="Short definitions of the terms used in ligentry.",
        keywords=("definition", "terms", "meaning", "what is", "abbreviation", "auc", "cmax", "rmsd", "smiles", "chembl", "pdb", "tlag", "infusion"),
        related=("units", "pkpd-model", "docking-overview"),
        body="""
<table @T@>
<tr><th>Term</th><th>Meaning</th></tr>
<tr><td>AUC</td><td>Area under the concentration-time curve; total drug exposure (mg&middot;h/L).</td></tr>
<tr><td>Bioavailability (F)</td><td>Fraction of an oral dose that reaches the blood.</td></tr>
<tr><td>ChEMBL</td><td>Public database of bioactive molecules and their targets.</td></tr>
<tr><td>Clearance (CL)</td><td>Volume of blood cleared of drug per hour (L/h).</td></tr>
<tr><td>Cmax, Tmax</td><td>Highest concentration and the time it occurs.</td></tr>
<tr><td>CV %</td><td>Standard error of a fitted value as a percentage of that value.</td></tr>
<tr><td>EC50</td><td>Concentration giving half of the maximum effect.</td></tr>
<tr><td>Emax</td><td>Maximum change in effect from baseline.</td></tr>
<tr><td>Exhaustiveness</td><td>How much effort Vina puts into searching for poses.</td></tr>
<tr><td>Flip-flop kinetics</td><td>Absorption slower than elimination, so the two rates can be confused when fitting oral data.</td></tr>
<tr><td>Hill coefficient</td><td>Steepness of the concentration-effect curve.</td></tr>
<tr><td>Hysteresis</td><td>Effect that differs at the same concentration on the way up and on the way down, because of a delay.</td></tr>
<tr><td>IV infusion</td><td>Constant zero-order delivery (dose / duration) into the central compartment over a set duration.</td></tr>
<tr><td>ke0</td><td>Rate constant of equilibration between plasma and the effect site (1/h).</td></tr>
<tr><td>Ligand / receptor</td><td>The small molecule / the protein it binds.</td></tr>
<tr><td>Michaelis-Menten elimination</td><td>Saturable removal: Vmax (mg/h) and Km (mg/L).</td></tr>
<tr><td>PDB</td><td>Protein Data Bank, the archive of 3D structures.</td></tr>
<tr><td>Pose</td><td>One position and shape of the ligand in the pocket.</td></tr>
<tr><td>RMSD</td><td>Root-mean-square distance between two sets of atoms (&Aring;). Below 2 &Aring; counts as a match in the benchmark.</td></tr>
<tr><td>RMSE</td><td>Root-mean-square error between data and model, in the data's unit.</td></tr>
<tr><td>SE</td><td>Standard error: the uncertainty of a fitted value.</td></tr>
<tr><td>SMILES</td><td>A line of text that describes a molecule's structure.</td></tr>
<tr><td>Stereoisomers</td><td>Molecules with the same bonds but different 3D arrangement.</td></tr>
<tr><td>Tlag</td><td>Oral absorption lag time (h) between administration and the start of absorption.</td></tr>
<tr><td>Vina score</td><td>Docking score in kcal/mol; more negative means a better predicted fit.</td></tr>
</table>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="units",
        title="Units",
        section="Reference",
        summary="The units used throughout the PK/PD simulator.",
        keywords=("unit", "mg", "litre", "hour", "mg/l"),
        related=("glossary", "pkpd-export"),
        body="""
<table @T@>
<tr><th>Quantity</th><th>Unit</th></tr>
<tr><td>Dose, amounts</td><td>mg</td></tr>
<tr><td>Volume (V, V2)</td><td>L</td></tr>
<tr><td>Clearance (CL, Q)</td><td>L/h</td></tr>
<tr><td>Rate constants (ka, ke0)</td><td>1/h</td></tr>
<tr><td>Vmax</td><td>mg/h</td></tr>
<tr><td>Concentration, Km, EC50</td><td>mg/L</td></tr>
<tr><td>Time, infusion duration, Tlag</td><td>h</td></tr>
<tr><td>AUC</td><td>mg&middot;h/L</td></tr>
<tr><td>Vina score</td><td>kcal/mol</td></tr>
</table>
<p>The units are fixed. If your data are in other units (ng/mL, minutes, micromolar), convert them before loading.
Note that 1 ng/mL equals 0.001 mg/L.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="cli",
        title="Command line",
        section="Reference",
        summary="The ligentry command and its subcommands.",
        keywords=("terminal", "command", "cli", "script", "batch", "redock", "benchmark", "plugins"),
        related=("benchmark", "install"),
        body="""
<pre>ligentry plugins
ligentry dock --pdb 1HSG --chembl CHEMBL243 --n 20 --out runs/hiv_protease
ligentry redock --pdb 1HSG --seed 42 --out runs/redock_1hsg
ligentry benchmark --seeds 42 7 123 2024 99 --out runs/benchmark</pre>
<ul>
<li><code>plugins</code> lists the installed tool wrappers (receptor prep, ligand prep, Vina, Uni-Dock, pose rescorer).</li>
<li><code>dock</code> runs a docking job like the main window does.</li>
<li><code>redock</code> redocks one crystal ligand and reports its RMSD.</li>
<li><code>benchmark</code> runs the full validation set (about 11 minutes on CPU, ~6 minutes on GPU).</li>
</ul>
<p>A run made on the command line can be opened later with <b>File &gt; Open run...</b>.</p>
""",
    ),
    Topic(
        id="files-settings",
        title="Files, settings and privacy",
        section="Reference",
        summary="What ligentry stores, where, and what it sends over the network.",
        keywords=("privacy", "data", "settings", "stored", "network", "internet", "cookies", "telemetry", "offline", "clear"),
        related=("docking-runs", "web-lookup"),
        body="""
<h3>What is stored</h3>
<ul>
<li><b>Run folders</b> are written where you choose (default <code>runs/gui</code> in the folder you started from).</li>
<li><b>Settings</b> (the recent-runs list and whether you have seen the welcome guide) are in
<code>~/.config/ligentry/ligentry.conf</code>. Delete the file to reset them.</li>
<li><b>Web tab data</b> (cookies, cache) is kept in your user data folder, as in any browser. Use
<b>Web &gt; Clear web data</b> in the Help window to remove it.</li>
</ul>
<h3>What goes over the network</h3>
<ul>
<li>When you run docking, ligentry downloads the structure from the RCSB PDB and, if you chose a ChEMBL target, ligands from ChEMBL.</li>
<li>The Help window's Web tab loads whichever page you ask for. Those sites can see your address and may set cookies.</li>
<li>Nothing else is sent anywhere. There is no usage tracking.</li>
</ul>
<p>The PK/PD simulator and this manual work with no internet connection.</p>
""",
    ),
    Topic(
        id="troubleshooting",
        title="Troubleshooting",
        section="Reference",
        summary="Common messages and what to do about them.",
        keywords=("error", "problem", "fix", "not working", "failed", "message", "blank", "crash", "cannot", "no results", "not available"),
        related=("docking-results", "pkpd-data-format", "files-settings"),
        body="""
<h3>Docking</h3>
<table @T@>
<tr><th>Message or symptom</th><th>What to do</th></tr>
<tr><td>Enter a 4-character PDB ID.</td><td>Type a four-character code such as 1HSG.</td></tr>
<tr><td>Enter a ChEMBL target or paste SMILES.</td><td>Fill in one of the two ligand sources.</td></tr>
<tr><td>No ligands to dock.</td><td>The ChEMBL target returned nothing, or the SMILES box had no usable lines.</td></tr>
<tr><td>Score shows n/a, status begins with failed</td><td>That ligand could not be prepared or docked; the status says why, for example a SMILES that RDKit cannot parse.</td></tr>
<tr><td>Docking is still running...</td><td>Press Stop, wait for it to finish, then close.</td></tr>
<tr><td>No results.csv in the folder</td><td>Choose the run folder itself (the one containing results.csv), not a parent.</td></tr>
<tr><td>The 3D viewer is unavailable</td><td>The run folder has no <code>&lt;PDB&gt;_protein.pdb</code> file. The table still works.</td></tr>
</table>
<h3>PK/PD</h3>
<table @T@>
<tr><th>Message or symptom</th><th>What to do</th></tr>
<tr><td>Cannot simulate: ...</td><td>A value is outside the allowed range. The message names it.</td></tr>
<tr><td>Line N: ... is not a number</td><td>Fix that line of your data file.</td></tr>
<tr><td>The optimiser did not converge</td><td>Try other starting values, fewer free parameters, or check the dosing.</td></tr>
<tr><td>Too few observations for uncertainty</td><td>You need more data points than fitted parameters.</td></tr>
<tr><td>All concentrations at the measurement times are zero</td><td>The effect times are before the dose, or the PK values are wrong.</td></tr>
<tr><td>A fitted value shows differently in the box</td><td>The box has a limited range or number of decimals; a note says which values changed.</td></tr>
</table>
<h3>Help window</h3>
<ul>
<li><b>The Web tab opens pages in your normal browser instead.</b> The embedded web view is not available in this installation.
The manual still works. If PySide6 was installed without its web engine component, install <code>PySide6-Addons</code> into the same
environment with <code>python -m pip install PySide6-Addons</code> (same version as PySide6: check with <code>python -m pip show PySide6</code>).</li>
<li><b>The Web tab is blank or closes.</b> Start ligentry with <code>QTWEBENGINE_CHROMIUM_FLAGS=&quot;--disable-gpu&quot; ligentry-gui</code>, or with
<code>QTWEBENGINE_DISABLE_SANDBOX=1</code> as a last resort.</li>
<li><b>A site will not load.</b> You may be offline, or the site may be refusing your network. Some sites (Google, for example) answer &quot;unusual traffic&quot; on certain networks. Try <b>Open in browser</b>, or use another source such as PubMed or Europe PMC.</li>
</ul>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="credits",
        title="About, license and credits",
        section="Reference",
        summary="License and the open-source projects ligentry builds on.",
        keywords=("license", "licence", "apache", "credits", "about", "thanks", "authors", "open source", "citation"),
        related=("welcome",),
        body="""
<p>ligentry is released under the Apache-2.0 license (see the <code>LICENSE</code> file). It builds on AutoDock Vina, Uni-Dock, RDKit,
Meeko, Open Babel, PySide6, PyVista, VTK, NumPy, SciPy, PyTorch, ONNX Runtime and Matplotlib, and on the RCSB PDB and ChEMBL databases. Each has its
own license, so check them before redistributing a bundled copy.</p>
<p>If you use ligentry in published work, please also cite the tools and databases it uses, in particular AutoDock Vina for docking.</p>
""",
    ),
    # ------------------------------------------------------------------ Research and the web
    Topic(
        id="web-lookup",
        title="Searching the web from Help",
        section="Research and the web",
        summary="Use the search box to search the manual, papers, medical sites, databases and videos.",
        keywords=("search", "browser", "web", "internet", "papers", "journal", "pubmed", "videos", "youtube", "doi", "pmid",
                  "research", "look up", "online"),
        related=("web-sources", "web-reading", "files-settings", "troubleshooting"),
        body="""
<p>The Help window works like a small browser. It has two tabs: <b>Manual</b> (this text, offline) and <b>Web</b>.</p>
<h3>The search box</h3>
<p>Type in the box at the top, choose <b>where to search</b> in the list beside it, and press Enter:</p>
<ul>
<li><b>This manual</b> searches every page here and shows the best matches.</li>
<li><b>PubMed, Europe PMC, PubMed Central</b> search research papers.</li>
<li><b>ClinicalTrials.gov, DailyMed, MedlinePlus</b> cover trials, drug labels and health information.</li>
<li><b>PubChem, ChEMBL, UniProt, RCSB PDB</b> cover molecules, targets and structures.</li>
<li><b>Wikipedia, YouTube, Web search</b> give overviews, explanatory videos and general results.</li>
</ul>
<h3>Shortcuts you can type</h3>
<table @T@>
<tr><th>You type</th><th>You get</th></tr>
<tr><td><code>PMID 12345678</code></td><td>That PubMed record</td></tr>
<tr><td><code>PMC1234567</code></td><td>That free full-text article</td></tr>
<tr><td><code>10.1000/xyz123</code> or <code>doi:10.1000/xyz123</code></td><td>The paper with that DOI</td></tr>
<tr><td><code>NCT01234567</code></td><td>That clinical study</td></tr>
<tr><td><code>PDB 1HSG</code></td><td>That structure at the RCSB PDB</td></tr>
<tr><td>a web address such as <code>who.int</code></td><td>That page</td></tr>
</table>
<h3>The Web tab</h3>
<p>It has Back, Forward, Reload, Stop, an address bar and <b>Open in browser</b> (to continue in your normal browser). PDF papers open
in the page, and downloads ask where to save the file. Links on any manual page that point to the web open here. Only ordinary
public http and https pages open; local files and addresses on your own network are refused.</p>
<p>Web searching is something you start. ligentry sends nothing to these sites on its own, and the manual always works offline.</p>
""".replace("@T@", _TABLE),
    ),
    Topic(
        id="web-sources",
        title="Where to look: reliable sources",
        section="Research and the web",
        summary="A short guide to the main sites for papers, drugs, molecules, proteins and learning.",
        keywords=("sources", "databases", "journals", "papers", "drug labels", "learning", "courses", "videos", "preprint", "reference"),
        related=("web-lookup", "web-reading"),
        web=("pharmacokinetics lecture", "molecular docking tutorial", "pharmacodynamics explained"),
        body="""
<h3>Research papers</h3>
<ul>
<li><a href="https://pubmed.ncbi.nlm.nih.gov/">PubMed</a>: citations and abstracts of biomedical papers. Free to search.</li>
<li><a href="https://europepmc.org/">Europe PMC</a>: papers with free full text where available, and preprints.</li>
<li><a href="https://pmc.ncbi.nlm.nih.gov/">PubMed Central</a>: free full-text archive.</li>
</ul>
<h3>Drugs, trials and health information</h3>
<ul>
<li><a href="https://dailymed.nlm.nih.gov/">DailyMed</a>: current US drug labels.</li>
<li><a href="https://clinicaltrials.gov/">ClinicalTrials.gov</a>: registry of clinical studies.</li>
<li><a href="https://medlineplus.gov/">MedlinePlus</a>: plain-language health information.</li>
</ul>
<h3>Molecules, targets and structures</h3>
<ul>
<li><a href="https://pubchem.ncbi.nlm.nih.gov/">PubChem</a> and <a href="https://www.ebi.ac.uk/chembl/">ChEMBL</a>: chemistry and bioactivity.</li>
<li><a href="https://www.uniprot.org/">UniProt</a>: proteins.</li>
<li><a href="https://www.rcsb.org/">RCSB PDB</a>: 3D structures (the source of your docking receptors).</li>
</ul>
<h3>Learning</h3>
<ul>
<li>Search YouTube or the web for lecture series on pharmacokinetics, pharmacodynamics or molecular docking, and prefer
courses from universities or established educators.</li>
<li>Each manual page ends with ready-made searches for its topic.</li>
</ul>
""",
    ),
    Topic(
        id="web-reading",
        title="Reading medical and scientific sources critically",
        section="Research and the web",
        summary="Peer review, preprints, study types and why a web page is not medical advice.",
        keywords=("evidence", "peer review", "preprint", "reliability", "misinformation", "critical", "quality", "bias", "citation"),
        related=("safety", "web-lookup", "web-sources"),
        web=("how to read a scientific paper", "levels of evidence medicine"),
        body="""
<ul>
<li><b>Check where it comes from.</b> A journal paper, a regulator's page and a personal blog are different kinds of source.</li>
<li><b>Peer review is a filter, not a guarantee.</b> A <b>preprint</b> (for example on bioRxiv or medRxiv) has <i>not</i> been peer reviewed.</li>
<li><b>Study type matters.</b> A test-tube or animal result does not show an effect in people. Within human studies,
randomised controlled trials and systematic reviews carry more weight than single case reports.</li>
<li><b>Look at the size and the numbers.</b> Small studies and results that nobody has repeated deserve caution.</li>
<li><b>Look for conflicts of interest and the date.</b> Medical knowledge changes.</li>
<li><b>Videos and summaries</b> are shortcuts to the idea; go to the paper or label for the details.</li>
</ul>
<p>Nothing you find here replaces advice from a doctor or pharmacist about your own health or medicines.</p>
""",
    ),
)

_BY_ID = {t.id: t for t in TOPICS}


def get_topic(topic_id: str) -> Topic:
    try:
        return _BY_ID[topic_id]
    except KeyError:
        raise KeyError(f"no help topic '{topic_id}'") from None


def has_topic(topic_id: str) -> bool:
    return topic_id in _BY_ID


def topics_in_section(section: str) -> tuple[Topic, ...]:
    return tuple(t for t in TOPICS if t.section == section)
