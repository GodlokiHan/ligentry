# Ligentry

Open-source desktop app for computational drug design: GPU-accelerated molecular docking with an AI neural pose rescorer, an interactive 3D viewer, and a built-in validation benchmark, plus a pharmacokinetic/pharmacodynamic (PK/PD) simulator with sequential and simultaneous (joint) fitting, a built-in manual with a first-launch guide, and a small browser for looking up papers and databases. Protein-structure work is planned.

**Status: early (v0.1.0).** Linux only for now. Windows is not supported yet.

## What it does

### Docking & AI Pose Rescoring

- Dock small-molecule libraries into a protein target taken from the PDB. The binding pocket is defined by the ligand already bound in the crystal structure.
- **Dual docking engines:** runs **Uni-Dock** (NVIDIA GPU batch docking) when installed or **AutoDock Vina** (CPU fallback), paired with multi-core parallel 3D conformer preparation (RDKit + Meeko) and Open Babel receptor protonation at pH 7.4.
- **Two-tier AI neural pose rescorer (ONNX Runtime):**
  - **CPU default model (`cpu_v1`):** bundled out of the box (~30 KB) and runs on any machine without a GPU.
  - **Optional GPU model pack (`gpu_v1`):** deeper neural interaction network installable via **Tools > AI pose rescorer models...** (or from a local `.tar.gz` release archive) for NVIDIA GPUs.
  - Re-ranks candidate poses by predicted native-pose probability (`ai_pose_score`) and estimates binding affinity (`pred_pkd`) alongside Vina scores (`vina_score_kcal_mol`).
- **Flexible ligand sources:** fetch known actives from ChEMBL, load large `.smi` / `.csv` ligand libraries directly from disk (up to **100,000 ligands**), or paste SMILES strings.
- Ranked results table, an interactive 3D view of the pose in its pocket, per-ligand notes, and save/open/export of whole runs.
- Unspecified stereocentres are expanded into explicit stereoisomers (capped, reproducible), and the best-scoring isomer represents each ligand.
- Docking runs in a separate process, so the window stays responsive and **Stop** keeps everything finished so far.
- A redocking benchmark re-docks crystal ligands from known complexes and measures heavy-atom RMSD against the crystal pose.

### PK/PD simulator

- One- or two-compartment model (central plus an optional peripheral compartment, set by its volume V2 and the intercompartmental clearance Q).
- Oral (with optional absorption lag time `Tlag`), intravenous (IV) bolus, or zero-order IV infusion dosing, single or repeated doses.
- Linear (clearance) or Michaelis-Menten elimination, with either number of compartments.
- Half-lives reported for the current parameters: terminal and distribution half-life for two compartments, terminal only for one.
- Sigmoid Emax effect model (a negative Emax models inhibition), either **direct** (the effect follows the plasma concentration at once) or **delayed** through an effect compartment: the effect follows a hypothetical effect-site concentration Ce that lags plasma, dCe/dt = ke0 (C - Ce). The delay is described by ke0 (1/h) and its equilibration half-life ln2/ke0, which is shown with the results.
- **Sequential or Joint PK/PD fitting:** fit PK parameters (`ka, V, CL, Vmax, Km, V2, Q, Tlag`) to plasma concentrations and PD parameters (`E0, Emax, EC50, Hill, ke0`) to measured effects either sequentially or **simultaneously (jointly)** via **Fit PK + PD jointly**, balancing concentration and effect residual scales and propagating PK uncertainty into PD estimates.
- Export the simulation as CSV (time, concentration, gut, central and, for two compartments, peripheral amounts, plus the effect-site concentration for a delayed effect and the effect when it is switched on).

Units throughout: dose in mg, volume in L, clearance in L/h, rate constants in 1/h, time in h, concentration in mg/L.

### Help, guide and research lookup

- **Welcome guide** on the first launch (reopen it from **Help > Welcome guide**).
- **Manual** (Help > ligentry manual, or **F1**, which opens the page for the window you are in): 38 offline pages covering every feature, a worked example with sample data, a glossary, units, keyboard shortcuts and troubleshooting. It has full-text search, back/forward and a contents page.
- **Web tab** in the Help window, an embedded browser for papers, drug and protein databases and videos. The search box can search the manual, or PubMed, Europe PMC, PubMed Central, ClinicalTrials.gov, DailyMed, MedlinePlus, PubChem, ChEMBL, UniProt, RCSB PDB, Wikipedia, YouTube or the web. Typing a DOI, `PMID 12345678`, `PMC1234567`, `NCT01234567`, `PDB 1HSG` or a web address goes straight to that page.
- Nothing is loaded from the web until you search or click a web link, only public http/https pages open, and the manual works offline. If the embedded browser is unavailable, pages open in your normal browser. Some sites refuse requests from some networks (Google Scholar did in testing, so it is not offered).

## Install (Linux)

You need [micromamba](https://mamba.readthedocs.io/) (or conda). A GPU is not required for CPU Vina docking and the bundled `cpu_v1` AI rescorer.

```fish
git clone https://github.com/GodlokiHan/ligentry.git
cd ligentry
micromamba env create -f envs/environment.yml --override-channels -c conda-forge
micromamba activate ligentry
pip install -e . --no-deps
python scripts/install_launcher.py   # optional: adds "ligentry" to your application menu
```

The environment file pins the versions the benchmark baseline was measured with. Install the Qt and VTK packages from the pip entries in that file, not from conda: the conda builds could not reach the system OpenGL on the development machine. Use `--no-deps` with `pip install -e .` so pip does not replace the pinned packages. The pip `PySide6` package includes the web engine used by the Help window's Web tab; if it is missing, pages open in your normal browser instead.

### Optional: NVIDIA GPU Acceleration (Uni-Dock + GPU AI Rescorer)

1. **Uni-Dock GPU engine:** install `unidock` in a dedicated micromamba environment (`micromamba create -n unidock -c conda-forge unidock`) and symlink the `unidock` binary onto your `PATH`. Ligentry automatically detects `unidock` and configures its CUDA library path at runtime.
2. **GPU AI rescorer model pack:** open `ligentry-gui`, go to **Tools > AI pose rescorer models...**, and install `ligentry-rescore-gpu-v1.tar.gz` from the GitHub Releases page (or from a local archive).

## Use

### Docking (GUI)

Run `ligentry-gui`, or open **ligentry** from your application menu.

1. Enter a PDB ID with a bound ligand (for example `1HSG`).
2. Enter a ChEMBL target ID (for example `CHEMBL243`), click **Load .smi / .csv...** to select a ligand library file, or paste SMILES, one per line as `SMILES [name]`.
3. Choose the docking engine (**Auto**, **Uni-Dock GPU**, or **AutoDock Vina CPU**) and whether to re-rank poses with the AI neural rescorer.
4. Press **Run docking** (Ctrl+R). Select a row to inspect its 3D pose, AI score, predicted pKd, and Vina energy.

Use **File > Open run...** to reload an earlier run folder without docking again. Press **F1** at any time for the manual.

### PK/PD simulator

Open it from **Tools > PK/PD simulator** (or the toolbar button) in the docking window, or run it on its own:

```fish
python -m ligentry.gui.pkpd_window
```

1. Choose the route (oral, IV bolus, or IV infusion), then set the dose, regimen, infusion duration or oral lag time (`Tlag`), and PK parameters. Choose **Two compartments** to enable V2 and Q. The concentration curve updates as you change them.
2. Optionally load measured concentrations (a CSV with time in h and concentration in mg/L) and press **Fit parameters** to estimate the PK parameters. With IV bolus or IV infusion dosing, absorption rate `ka` and `Tlag` are disabled.
3. Tick **Effect** to add the effect curve. To model a lag between concentration and effect, tick **Delayed effect (effect compartment)** in the same box and set the effect-site rate ke0; the effect-site concentration is then drawn as a dashed curve on the concentration plot. Optionally load measured effects (a CSV with time in h and effect), choose which parameters to fit (tick **ke0** to fit a delayed effect; it is off by default) and press **Fit effect**.
4. When both concentration and effect datasets are loaded, **Fit PK + PD jointly** estimates all ticked PK and PD parameters simultaneously and propagates PK uncertainty into PD estimates. **Use fitted values**, **Use fitted effect values**, and **Use joint fitted values** copy the estimates back into the form.

Tips for fitting a two-compartment model: IV bolus data with early, rapidly falling samples identify V, CL, V2 and Q best. Oral data alone make the fit harder and give apparent values (V/F, CL/F, V2/F, Q/F).

Tips for fitting a delayed effect: measure the effect densely around and after the concentration peak, where the lag is visible. Compare the residuals of the direct and delayed models before trusting a fitted ke0.

The manual has a worked example with sample data you can save and fit (Help > ligentry manual > PK/PD simulator > Worked example: fitting data).

### Command line

```fish
ligentry plugins
ligentry dock --pdb 1HSG --chembl CHEMBL243 --n 20 --out runs/hiv_protease
ligentry redock --pdb 1HSG --seed 42 --out runs/redock_1hsg
ligentry benchmark --seeds 42 7 123 2024 99 --out runs/benchmark
```

## Validation

### Docking

The benchmark redocks the crystal ligand of 8 diverse complexes. A run passes when the top-ranked pose is within 2 Å heavy-atom RMSD of the crystal pose.

| PDB | Target / ligand | Vina CPU baseline (5 seeds) | Uni-Dock GPU + AI Rescorer (seed 42) |
|---|---|---|---|
| 3PTB | Trypsin / benzamidine | 5/5 (PASS) | PASS (0.35 Å) |
| 1STP | Streptavidin / biotin | 5/5 (PASS) | PASS (0.96 Å) |
| 1M17 | EGFR kinase / erlotinib | 5/5 (PASS) | PASS (1.66 Å) |
| 3ERT | Estrogen receptor α / 4-hydroxytamoxifen | 5/5 (PASS) | PASS (1.32 Å) |
| 1KV1 | p38 MAP kinase / pyrazole urea fragment | 5/5 (PASS) | PASS (0.69 Å) |
| 1HSG | HIV-1 protease / indinavir | 2/5 (FAIL) | FAIL (4.34 Å top; 1.20 Å in top 4) |
| 1EVE | Acetylcholinesterase / donepezil | 0/5 (FAIL) | PASS (0.85 Å) |
| 4DFR | DHFR / methotrexate | 0/5 (FAIL) | FAIL (8.31 Å) |

- **CPU Vina baseline:** 5 of 8 complexes pass (27 of 40 runs, 67.5%).
- **Uni-Dock GPU + CrossDocked2020-trained AI Rescorer:** **6 of 8 complexes pass (75.0% top-1 recovery, 1.14 Å median top-pose RMSD)**, rescuing `1EVE` (donepezil) to **0.85 Å** top-pose accuracy while completing all 8 targets in ~6 minutes on an NVIDIA RTX 2050.
- The measured baseline is in [`benchmarks/baseline_garuda.json`](benchmarks/baseline_garuda.json). 4DFR is documented in [`benchmarks/README.md`](benchmarks/README.md) as a scoring limitation.

### PK/PD

The simulator is checked against analytic solutions rather than a benchmark data set:

- IV bolus, IV infusion, and oral two-compartment curves match the closed-form solutions (including overlapping infusions and shifted oral lag time `Tlag`).
- Dose = central + peripheral + eliminated amount (mass balance), and AUC equals dose/CL (times F for oral dosing).
- Limiting cases: Q → 0 recovers the one-compartment model, and very fast distribution behaves as one compartment of volume V1 + V2, and very short IV infusions recover the IV bolus limit.
- Repeated doses obey superposition, and saturated Michaelis-Menten elimination removes Vmax mg/h from the whole body.
- The effect-site concentration matches its closed form after an IV bolus (including the limit where ke0 equals the elimination rate), its AUC equals the plasma AUC, a very large ke0 recovers the direct effect, and a delayed effect shows hysteresis: the same plasma concentration gives different effects on the way up and down.
- Sequential and joint fits recover known parameters from simulated data (to 0.1–1% without noise, within 25% with noise), including `Tlag` and `ke0`. Data from a direct effect push ke0 to its upper limit, which the fit notes report as no detectable delay.

These tests show the maths is implemented correctly. They do not show that any parameter set describes a real drug.

## Limitations

**Docking**

- Vina scores and AI rescorer predictions are screening heuristics, not measured binding affinities. A good docking score is a hypothesis, not evidence that a compound works.
- Ligands are docked in their neutral form. There is no protonation-state modelling.
- Crystal waters and cofactors are removed from the receptor, and the receptor is rigid.
- The pocket comes from the ligand bound in the PDB structure you choose.
- Stereoisomer enumeration multiplies docking time for ligands with many unspecified centres (default cap: 8 isomers per ligand).

**PK/PD**

- The output is only as good as the parameters and data you supply. It is a modelling tool, not a prediction of how a drug behaves in people, and nothing here is dosing advice.
- At most two compartments, oral, IV bolus, or IV infusion dosing only (no other extravascular routes), and the effect is either direct or a single effect compartment. Other delay models (for example indirect-response models or tolerance) are not available.
- From oral data alone, bioavailability F and volume V cannot be separated: the fitted volumes are V/F and V2/F, and clearances are CL/F and Q/F.
- Absorption and elimination can be swapped (flip-flop kinetics) when only oral data are fitted.
- A two-compartment fit needs data that show both phases. When the distribution and terminal phases are close together, or early samples are missing, V2 and Q are poorly determined and the fit notes will say so.
- EC50 and Emax trade off against each other when the measured concentrations never approach EC50, so the fit notes may warn about wide uncertainty.
- ke0 and EC50 are only as reliable as the PK model that supplies the plasma concentrations: a wrong PK model can look like a delay. With sparse effect data ke0 is poorly determined, and a ke0 that finishes on its limits (0.001 or 100 per hour) is a sign the data cannot determine it.

**Help and web**

- Pages opened in the Web tab come from other sites and are not reviewed by ligentry. Nothing in ligentry or its manual is medical advice.
- The search addresses of external sites may change over time; the plain web search always works.

## Project layout

```
ligentry/
  core/        plugin interface, registry, and AI model tier manager
  data/        PDB and ChEMBL access
  models/      bundled CPU default ONNX pose rescorer (cpu_v1)
  plugins/     receptor prep (Open Babel), ligand prep (RDKit + Meeko),
               docking (AutoDock Vina & Uni-Dock GPU), AI neural pose rescorer
  workflows/   docking, redocking, benchmark
  gui/         docking window, 3D viewer, background docking job, AI models dialog,
               PK/PD window, launcher, Help window, welcome guide dialog
  pkpd/        PK/PD models, sequential & joint fitting, data loaders, residuals (no Qt dependency)
  help/        manual pages, search, page rendering, web sources, guide logic (no Qt dependency)
  cli.py       command-line interface
scripts/       tests, launcher installer, ML dataset builder & rescorer trainer, 4DFR investigation
benchmarks/    benchmark documentation and baseline
```

Each wrapped docking tool is a plugin behind a common interface, so engines can be swapped without touching the workflows. Scientific logic lives in Qt-free modules with their own test scripts; the windows stay thin.

## Tests

```fish
for t in scripts/test_*.py; echo $t; python $t | tail -2; end
```

Each script ends with `ALL PASSED` or a count such as `1333 passed, 0 failed`. The PK/PD scripts are `test_pkpd`, `test_twocomp`, `test_effect`, `test_effect_compartment`, `test_observations`, `test_residuals`, `test_effect_residuals`, `test_pkpd_workflow`, `test_joint_fit`, `test_infusion_tlag`, `test_pkpd_window` and `test_pkpd_window_delayed`. They need only NumPy and SciPy, apart from the two window tests, which also need Qt.

The help scripts are `test_help` (manual, search, web helpers and guide logic; no Qt; it also checks that the manual still matches the program's shortcuts, defaults and CSV columns) and `test_help_gui` (Help window, guide dialog and Help menus, offscreen). To try the embedded web view in the second one too: `env LIGENTRY_TEST_WEB=1 python scripts/test_help_gui.py`.

## Roadmap

- PK/PD: extravascular routes beyond oral, transit-compartment absorption, SBML import.
- Docking: native PDBQT domain fine-tuning for the neural pose rescorer (aimed at `1HSG` and `4DFR`), OpenMM pose refinement, a protonation step if it improves the full benchmark, command-line flags for stereoisomer options.
- Protein-structure workflows.
- Replace the Open Babel subprocess with a pure-Python route, and port to Windows.

## Privacy

ligentry has no usage tracking. It contacts the network only when you run docking (to download the structure from the RCSB PDB and ligands from ChEMBL), when you explicitly click **Download & install GPU pack** in the AI Models dialog, and when you search or open a page in the Help window. Settings (recent runs, whether the guide was shown) are stored in `~/.config/ligentry/ligentry.conf`, and optional downloaded AI model packs are stored in `~/.local/share/ligentry/models/`. The Web tab keeps cookies and cache like any browser; **Web > Clear web data** in the Help window removes them.

## License and credits

Ligentry is released under the Apache-2.0 license (see `LICENSE`). It builds on AutoDock Vina, Uni-Dock, RDKit, Meeko, Open Babel, PySide6, PyVista, VTK, NumPy, SciPy, PyTorch, ONNX Runtime, Matplotlib, and the RCSB PDB and ChEMBL databases. Each has its own license, so check them before redistributing a bundled copy.