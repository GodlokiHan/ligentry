# Ligentry

Open-source desktop app for computational drug design: molecular docking with a 3D viewer and a built-in validation benchmark, plus a pharmacokinetic/pharmacodynamic (PK/PD) simulator. Protein-structure work is planned.

**Status: early (v0.1.0).** Linux only for now. Windows is not supported yet.

## What it does

### Docking

- Dock a set of ligands into a protein target taken from the PDB. The binding pocket is defined by the ligand already bound in the crystal structure.
- Ligands come from ChEMBL (known actives for a target) or pasted SMILES.
- Ranked results table, an interactive 3D view of the pose in its pocket, per-ligand notes, and save/open/export of whole runs.
- Unspecified stereocentres are expanded into explicit stereoisomers (capped, reproducible), and the best-scoring isomer represents each ligand.
- Docking runs in a separate process, so the window stays responsive and **Stop** keeps everything finished so far.
- A redocking benchmark re-docks crystal ligands from known complexes and measures RMSD against the crystal pose.

Docking uses [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina). Ligands are prepared with RDKit and Meeko, receptors with Open Babel.

### PK/PD simulator

- One-compartment oral model with linear or Michaelis-Menten elimination, single or repeated doses.
- Direct sigmoid Emax effect model (a negative Emax models inhibition).
- Load measured concentration data and fit the PK parameters, with standard errors, residuals and RMSE shown beside the plot.
- Load measured effect data and fit E0, Emax, EC50 and (optionally) the Hill coefficient, using the concentrations predicted by the PK model. Effect residuals are plotted too: a smooth up-then-down pattern suggests a delayed effect that this model cannot capture.
- Export the simulation as CSV.

Units throughout: dose in mg, volume in L, clearance in L/h, rate constants in 1/h, time in h, concentration in mg/L.

## Install (Linux)

You need [micromamba](https://mamba.readthedocs.io/) (or conda). A GPU is not required.

```fish
git clone https://github.com/GodlokiHan/ligentry.git
cd ligentry
micromamba env create -f envs/environment.yml --override-channels -c conda-forge
micromamba activate ligentry
pip install -e . --no-deps
python scripts/install_launcher.py   # optional: adds "ligentry" to your application menu
```

The environment file pins the versions the benchmark baseline was measured with. Install the Qt and VTK packages from the pip entries in that file, not from conda: the conda builds could not reach the system OpenGL on the development machine. Use `--no-deps` with `pip install -e .` so pip does not replace the pinned packages.

## Use

### Docking (GUI)

Run `ligentry-gui`, or open **ligentry** from your application menu.

1. Enter a PDB ID with a bound ligand (for example `1HSG`).
2. Enter a ChEMBL target ID (for example `CHEMBL243`) or paste SMILES, one per line as `SMILES [name]`.
3. Press **Run docking** (Ctrl+R). Select a row to see its pose.

Use **File > Open run...** to reload an earlier run folder without docking again.

### PK/PD simulator

Open it from **Tools > PK/PD simulator** (or the toolbar button) in the docking window, or run it on its own:

```fish
python -m ligentry.gui.pkpd_window
```

1. Set the dose, regimen and PK parameters. The concentration curve updates as you change them.
2. Optionally load measured concentrations (a CSV with time in h and concentration in mg/L) and press **Fit** to estimate the PK parameters.
3. Tick **Effect** to add the effect curve. Optionally load measured effects (a CSV with time in h and effect), choose which parameters to fit and press **Fit effect**. Fit the PK first: the effect fit takes concentrations from the PK model currently in the form.
4. **Use fitted effect values** copies the estimates back into the form.

### Command line

```fish
ligentry plugins
ligentry dock --pdb 1HSG --chembl CHEMBL243 --n 20 --out runs/hiv_protease
ligentry redock --pdb 1HSG --seed 42 --out runs/redock_1hsg
ligentry benchmark --seeds 42 7 123 2024 99 --out runs/benchmark   # about 11 minutes
```

## Validation (docking)

The benchmark redocks the crystal ligand of 8 diverse complexes with 5 seeds each. A run passes when the top-ranked pose is within 2 Å heavy-atom RMSD of the crystal pose, and a complex passes when most seeds do.

| PDB | Target / ligand | Seeds passed |
|---|---|---|
| 3PTB | Trypsin / benzamidine | 5/5 |
| 1STP | Streptavidin / biotin | 5/5 |
| 1M17 | EGFR kinase / erlotinib | 5/5 |
| 3ERT | Estrogen receptor α / 4-hydroxytamoxifen | 5/5 |
| 1KV1 | p38 MAP kinase / pyrazole urea fragment | 5/5 |
| 1HSG | HIV-1 protease / indinavir | 2/5 |
| 1EVE | Acetylcholinesterase / donepezil | 0/5 |
| 4DFR | DHFR / methotrexate | 0/5 |

Overall: 5 of 8 complexes pass, and 27 of 40 individual runs (67.5%). For 7 of 8 complexes, a near-native pose exists among the nine returned poses, so most failures are ranking failures of the scoring function rather than sampling failures. 4DFR is documented in [`benchmarks/README.md`](benchmarks/README.md) as a Vina scoring limitation.

The measured baseline is in [`benchmarks/baseline_garuda.json`](benchmarks/baseline_garuda.json). Results are deterministic for a given environment and seed but not across RDKit versions, so compare multi-seed pass rates rather than single runs.

## Limitations

**Docking**

- Vina scores are a screening heuristic, not binding affinities. A good docking score is a hypothesis, not evidence that a compound works.
- Ligands are docked in their neutral form. There is no protonation-state modelling.
- Crystal waters and cofactors are removed from the receptor, and the receptor is rigid.
- The pocket comes from the ligand bound in the PDB structure you choose.
- Stereoisomer enumeration multiplies docking time for ligands with many unspecified centres (default cap: 8 isomers per ligand).

**PK/PD**

- The output is only as good as the parameters and data you supply. It is a modelling tool, not a prediction of how a drug behaves in people.
- One compartment and direct effect only. A delayed effect shows up as a pattern in the effect residuals but is not modelled.
- From oral data alone, bioavailability F and volume V cannot be separated: the fitted volume is V/F and the fitted clearance is CL/F.
- Absorption and elimination can be swapped (flip-flop kinetics) when only oral data are fitted.
- EC50 and Emax trade off against each other when the measured concentrations never approach EC50, so the fit notes may warn about wide uncertainty.

## Project layout

```
ligentry/
  core/        plugin interface and registry
  data/        PDB and ChEMBL access
  plugins/     receptor prep (Open Babel), ligand prep (RDKit + Meeko), docking (Vina)
  workflows/   docking, redocking, benchmark
  gui/         docking window, 3D viewer, background docking job, PK/PD window, launcher
  pkpd/        PK/PD models, fitting, data loaders, residuals (no Qt dependency)
  cli.py       command-line interface
scripts/       tests, launcher installer, 4DFR investigation
benchmarks/    benchmark documentation and baseline
```

Each wrapped docking tool is a plugin behind a common interface, so engines can be swapped without touching the workflows. Scientific logic lives in Qt-free modules with their own test scripts; the windows stay thin.

## Tests

```fish
for t in scripts/test_*.py; echo $t; python $t | tail -2; end
```

Each script ends with `ALL PASSED` (or a count of failures). The window code is not covered by automated tests yet.

## Roadmap

- PK/PD: two-compartment model, delayed (effect-compartment) response, joint PK and PD fitting, SBML import.
- Docking: OpenMM pose refinement and rescoring (aimed at the 1HSG and 1EVE ranking failures), a protonation step if it improves the full benchmark, command-line flags for stereoisomer options.
- Protein-structure workflows.
- Replace the Open Babel subprocess with a pure-Python route, and port to Windows.

## License and credits

Ligentry is released under the Apache-2.0 license (see `LICENSE`). It builds on AutoDock Vina, RDKit, Meeko, Open Babel, PySide6, PyVista, VTK, NumPy, SciPy, Matplotlib, and the RCSB PDB and ChEMBL databases. Each has its own license, so check them before redistributing a bundled copy.
