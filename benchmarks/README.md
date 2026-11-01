# ligentry redocking benchmark

Redocking of 8 crystal ligands across diverse targets. Each ligand is rebuilt
from its RCSB chemical component SMILES (random 3D conformer, stereochemistry
taken from the crystal), docked with AutoDock Vina into a box built around the
native ligand, and compared to the crystal pose by symmetry-aware heavy-atom
RMSD (no re-alignment). A run passes when the top-ranked pose is within 2.0 A.
Each complex is run with several seeds (the seed drives both the RDKit
conformer embedding and Vina); a complex passes when the top pose passes in a
majority of seeds.

Run it:

```
ligentry benchmark --seeds 42 7 123 2024 99 --out runs/benchmark
```

Subset: `ligentry benchmark --only 3PTB 1STP --seeds 42 7 --out runs/bench_subset`

## Baseline (Garuda Linux, `baseline_garuda.json`)

Settings: exhaustiveness 16, 9 poses, seeds 42 7 123 2024 99, CPU only,
RDKit 2025.09.6 (pinned in `envs/environment.yml`).

| PDB  | Ligand | Target                                | Seeds passed | Median top RMSD (A) | Best any pose (A) | Result |
|------|--------|---------------------------------------|--------------|---------------------|-------------------|--------|
| 1HSG | MK1    | HIV-1 protease / indinavir            | 2/5          | 4.357               | 0.819             | FAIL   |
| 3PTB | BEN    | Trypsin / benzamidine                 | 5/5          | 0.390               | 0.327             | PASS   |
| 1STP | BTN    | Streptavidin / biotin                 | 5/5          | 0.940               | 0.744             | PASS   |
| 1M17 | AQ4    | EGFR kinase / erlotinib               | 5/5          | 1.795               | 1.359             | PASS   |
| 3ERT | OHT    | Estrogen receptor alpha / 4-OHT       | 5/5          | 1.377               | 1.230             | PASS   |
| 4DFR | MTX    | DHFR / methotrexate                   | 0/5          | 5.761               | 4.942             | FAIL   |
| 1KV1 | BMU    | p38 MAP kinase / pyrazole urea (BMU)  | 5/5          | 0.642               | 0.532             | PASS   |
| 1EVE | E20    | Acetylcholinesterase / donepezil      | 0/5          | 3.474               | 0.867             | FAIL   |

- Complexes passing by majority of seeds: 5/8
- Per-run top-pose pass rate: 27/40 (0.675)
- Complexes with a near-native pose (<= 2.0 A) among any pose of any seed: 7/8
- Median of per-complex median top RMSD: 1.586 A

"Best any pose" is the minimum over all poses of all seeds.

Results are deterministic for a given environment and seed, but not across
RDKit versions: the same single-seed (42) benchmark gave 6/8 on RDKit
2026.03.1 and 5/8 on RDKit 2025.09.6, with 1HSG, 1EVE and 1M17 changing
between pass and fail. Pin the environment and compare multi-seed pass rates,
not single runs.

## Known failures

### 1HSG and 1EVE: ranking, not sampling

Vina finds near-native poses (best 0.82 A and 0.87 A) but often does not rank
them first. Treated as a known scoring limitation.

### 4DFR (DHFR / methotrexate): scoring limitation

4DFR fails in every seed (0/5, median top-pose RMSD 5.76 A, best pose of any
rank 4.94 A). We investigated it and attribute the failure to Vina's scoring,
not to a bug in ligand preparation or in the search. It stays in the suite as a
documented hard case, and the pipeline is deliberately not tuned to pass it.

**Diagnostics.** Same redocking, five seeds (42, 7, 123, 2024, 99),
exhaustiveness 16 unless noted. Reproduce with `python scripts/diag_4dfr.py`.

| Variant | Ligand state                      | Seeds passed | Median top RMSD (A) | Best any pose (A) |
|---------|-----------------------------------|--------------|---------------------|-------------------|
| A       | neutral (baseline)                | 0/5          | 5.761               | 4.942             |
| B       | carboxylates deprotonated         | 0/5          | 5.952               | 2.690             |
| C       | B + pteridine N1-H+ (bound state) | 0/5          | 8.634               | 5.719             |
| D       | neutral, exhaustiveness 32        | 0/5          | 5.838               | 4.942             |

- Search depth is not the cause. Raising exhaustiveness to 32 (D) finds the
  same best pose as the baseline in two seeds and the same wrong top pose
  repeatedly. An earlier single-seed test at exhaustiveness 64 gave the same
  picture (top pose 5.858 A, best of 9 poses 5.246 A).
- Ligand protonation helps sampling but not ranking. Charging the carboxylates
  (B) brings the closest pose to 2.7-2.85 A in four of five seeds, but the
  top-ranked pose is still about 6 A away. Protonating N1 (C) makes it worse.
- Receptor preparation looks sound. Asp27 is anionic in the receptor PDBQT (no
  hydrogen on either carboxylate oxygen). The 4DFR file contains no cofactor.
  Only `ATOM` records are kept, so the 428 crystallographic waters and the ions
  are not part of the docking receptor.

**Where the error is.** `python scripts/analyze_4dfr_poses.py` splits the RMSD
by fragment. In most baseline top poses the pteridine head is about 8 A from
the crystal position (median 7.96 A) while the benzoyl linker and glutamate
tail are about 3.7 A off. In variant B's closest poses the benzoyl is 0.4-0.9 A
and the glutamate about 2.3-2.6 A from the crystal, but the pteridine is still
about 4 A away. The failure is mainly the deep pteridine pocket, not only a
floppy tail.

**Scoring, not sampling.** `python scripts/rescore_4dfr.py` scores the crystal
pose and the docked poses with the same Vina call. The crystal pose scores
worse than the docked top pose in all 10 runs (variants A and B, five seeds
each):

| Variant | Gap, native minus docked top (kcal/mol) | Median |
|---------|-----------------------------------------|--------|
| A       | +0.90 to +1.32                          | +1.24  |
| B       | +0.55 to +1.75                          | +0.96  |

The intermolecular term alone also favours the wrong pose (crystal about -12.3
kcal/mol against -13.2 to -15.1 for the docked top poses), so this is not an
intramolecular-strain artefact. The crystal pose is a stable Vina minimum:
local minimisation moves it only about 0.3 A. Vina therefore ranks wrong poses
above the experimental one, and more sampling cannot fix that.

**Caveats and untested factors.**

- The gap (about 1 kcal/mol) is small compared with typical scoring-function
  error.
- Conserved active-site waters are removed from the receptor and may matter for
  this complex. Keeping waters near the ligand has not been tried.
- Vina has no explicit electrostatics, so a strongly charged ligand such as
  methotrexate is a weak case for it.
- The crystal-pose hydrogens were placed geometrically, not force-field
  relaxed.

**Decision.**

- 4DFR stays in `BENCHMARK_SET` as a known hard case.
- Ligand preparation keeps the neutral default. `run_redock(..., protonation=...)`
  still accepts `"anion"` and `"bound"` for experiments.
- Any future protonation step (for example Dimorphite-DL) should be adopted
  only if it improves the full multi-seed benchmark, not 4DFR alone.

## Notes

- 1KV1 uses ligand code `BMU`, a small pyrazole urea fragment, not BIRB-796.
- Crystal waters are not modelled (see the 4DFR caveats above).
- To extend the suite, add a row to `BENCHMARK_SET` in
  `ligentry/workflows/benchmark.py`.
