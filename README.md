# biosim

Open-source molecular docking with a desktop GUI, a 3D viewer and a built-in validation benchmark. It is the first stage of a planned platform for computational drug design: docking now, with drug-effectiveness modelling and protein-structure work planned.

**Status: early (v0.1.0).** Linux only for now. Windows is not supported yet.

## What it does

- Dock a set of ligands into a protein target taken from the PDB. The binding pocket is defined by the ligand already bound in the crystal structure.
- Ligands come from ChEMBL (known actives for a target) or pasted SMILES.
- Ranked results table, an interactive 3D view of the pose in its pocket, per-ligand notes, and save/open/export of whole runs.
- Unspecified stereocentres are expanded into explicit stereoisomers (capped, reproducible), and the best-scoring isomer represents each ligand.
- Docking runs in a separate process, so the window stays responsive and **Stop** keeps everything finished so far.
- A redocking benchmark that re-docks crystal ligands from known complexes and measures RMSD against the crystal pose.

Docking uses [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina). Ligands are prepared with RDKit and Meeko, receptors with Open Babel.

## Install (Linux)

You need [micromamba](https://mamba.readthedocs.io/) (or conda). A GPU is not required.

```fish
git clone <repo-url> biosim
cd biosim
micromamba env create -f envs/environment.yml --override-channels -c conda-forge
micromamba activate biosim
pip install -e .
python scripts/install_launcher.py   # optional: adds "biosim" to your application menu
```

The environment file pins the versions the benchmark baseline was measured with. Install the Qt and VTK packages from the pip entries in that file, not from conda: the conda builds could not reach the system OpenGL on the development machine.

## Use

**GUI:** run `biosim-gui`, or open **biosim** from your application menu.

1. Enter a PDB ID with a bound ligand (for example `1HSG`).
2. Enter a ChEMBL target ID (for example `CHEMBL243`) or paste SMILES, one per line as `SMILES [name]`.
3. Press **Run docking** (Ctrl+R). Select a row to see its pose.

**Command line:**

```fish
biosim plugins
biosim dock --pdb 1HSG --chembl CHEMBL243 --n 20 --out runs/hiv_protease
biosim redock --pdb 1HSG --seed 42 --out runs/redock_1hsg
biosim benchmark --seeds 42 7 123 2024 99 --out runs/benchmark   # about 11 minutes
```

## Validation

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

- Vina scores are a screening heuristic, not binding affinities. Docking results are hypotheses, not evidence of activity.
- Ligands are docked in their neutral form. There is no protonation-state modelling.
- Crystal waters and cofactors are removed from the receptor, and the receptor is rigid.
- The pocket comes from the ligand bound in the PDB structure you choose.
- Stereoisomer enumeration multiplies docking time for ligands with many unspecified centres (default cap: 8 isomers per ligand).

## Project layout

```
biosim/
  core/        plugin interface and registry
  data/        PDB and ChEMBL access
  plugins/     receptor prep (Open Babel), ligand prep (RDKit + Meeko), docking (Vina)
  workflows/   docking, redocking, benchmark
  gui/         Qt window, 3D viewer, background docking job
  cli.py       command-line interface
scripts/       tests and the 4DFR investigation
benchmarks/    benchmark documentation and baseline
```

Each wrapped tool is a plugin behind a common interface, so engines can be swapped without touching the workflows.

## Tests

```fish
python scripts/test_stereo.py
python scripts/test_runio.py
python scripts/test_job.py
```

## Roadmap

Drug-effectiveness modelling (pharmacokinetic/pharmacodynamic models) and protein-structure workflows are planned but not started. Replacing the Open Babel subprocess with a pure-Python route and a Windows port are also open.

## License and credits

biosim is released under the Apache-2.0 license (see `LICENSE`). It builds on AutoDock Vina, RDKit, Meeko, Open Babel, PySide6, PyVista, VTK, and the RCSB PDB and ChEMBL databases. Each has its own license, so check them before redistributing a bundled copy.
