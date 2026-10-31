# biosim redocking benchmark

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
biosim benchmark --seeds 42 7 123 2024 99 --out runs/benchmark
```

Subset: `biosim benchmark --only 3PTB 1STP --seeds 42 7 --out runs/bench_subset`

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

### 4DFR (DHFR / methotrexate): ligand preparation

Fails in every seed. Raising exhaustiveness from 16 to 64 did not help (top
pose 5.858 A, best of 9 poses 5.246 A), so sampling is not the limiting factor.
The pipeline neutralizes ligands, while crystallographic methotrexate binds
with an ionized glutamate and (likely) a protonated pteridine ring that form
salt bridges in the pocket. This is the suspected cause and has not been
confirmed. Planned fix: protonation-aware ligand preparation.

## Notes

- 1KV1 uses ligand code `BMU`, a small pyrazole urea fragment, not BIRB-796.
- Crystal waters are not modelled.
- To extend the suite, add a row to `BENCHMARK_SET` in
  `biosim/workflows/benchmark.py`.
