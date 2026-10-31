# biosim redocking benchmark

Redocking of 8 crystal ligands across diverse targets. Each ligand is rebuilt
from its RCSB chemical component SMILES (random 3D conformer, stereochemistry
taken from the crystal), docked with AutoDock Vina into a box built around the
native ligand, and compared to the crystal pose by symmetry-aware heavy-atom
RMSD (no re-alignment). Success is a top-ranked pose within 2.0 A.

Run it:

```
biosim benchmark --out runs/benchmark
```

Subset: `biosim benchmark --only 3PTB 1STP --out runs/bench_subset`

## Baseline (Garuda Linux, `baseline_garuda.json`)

Settings: exhaustiveness 16, 9 poses, CPU only.

| PDB  | Ligand | Target                                | Top RMSD (A) | Best of 9 (A) | Result |
|------|--------|---------------------------------------|--------------|---------------|--------|
| 1HSG | MK1    | HIV-1 protease / indinavir            | 0.639        | 0.639         | PASS   |
| 3PTB | BEN    | Trypsin / benzamidine                 | 0.386        | 0.357         | PASS   |
| 1STP | BTN    | Streptavidin / biotin                 | 0.691        | 0.691         | PASS   |
| 1M17 | AQ4    | EGFR kinase / erlotinib               | 5.853        | 1.417         | FAIL   |
| 3ERT | OHT    | Estrogen receptor alpha / 4-OHT       | 1.432        | 1.205         | PASS   |
| 4DFR | MTX    | DHFR / methotrexate                   | 5.858        | 5.246         | FAIL   |
| 1KV1 | BMU    | p38 MAP kinase / pyrazole urea (BMU)  | 0.623        | 0.623         | PASS   |
| 1EVE | E20    | Acetylcholinesterase / donepezil      | 0.848        | 0.848         | PASS   |

- Top-pose success: 6/8
- Success with any of the 9 poses: 7/8
- Median top-pose RMSD: 0.769 A

Small RMSD differences are expected between machines and runs (Vina is
stochastic). A different pass/fail pattern is not.

## Known failures

### 1M17 (EGFR / erlotinib): scoring tie, not a sampling failure

Vina finds near-native poses but does not rank them first:

| Pose | Vina score (kcal/mol) | RMSD (A) |
|------|-----------------------|----------|
| 1    | -7.269                | 5.853    |
| 2    | -7.221                | 1.833    |
| 3    | -7.118                | 1.417    |

Poses 2 and 3 are within 0.15 kcal/mol of pose 1, far inside Vina's scoring
error. Treated as a known scoring limitation. Crystal waters are not modelled.

### 4DFR (DHFR / methotrexate): ligand preparation

Raising exhaustiveness from 16 to 64 did not help: top pose scored -8.912
kcal/mol at 5.858 A RMSD, and the best of 9 poses stayed at 5.246 A. Sampling
is therefore not the limiting factor.

The pipeline neutralizes ligands, while crystallographic methotrexate binds
with an ionized glutamate and (likely) a protonated pteridine ring that form
salt bridges in the pocket. This is the suspected cause and has not been
confirmed. Planned fix: protonation-aware ligand preparation.

## Notes

- 1KV1 uses ligand code `BMU`, a small pyrazole urea fragment, not BIRB-796.
- To extend the suite, add a row to `BENCHMARK_SET` in
  `biosim/workflows/benchmark.py`.
