# External GFN1-xTB/GFN2-xTB comparison source

Complete periodic GFN1-xTB and GFN2-xTB inputs, outputs, and analysis belong to
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
They are not duplicated in this repository.

Part I imports only the compact values needed for like-for-like comparisons:

| Benchmark | Protocol | GFN1-xTB MAE | GFN2-xTB MAE |
|---|---|---:|---:|
| DMC-ICE13 | fixed 3x3x3 | 8.005 kJ mol-1 H2O-1 | 3.463 kJ mol-1 H2O-1 |
| DMC-ICE13 | phase-wise adaptive, at most 4x4x4 | 8.006 kJ mol-1 H2O-1 | 3.461 kJ mol-1 H2O-1 |
| LC10 lattice constants | common ten-solid set | 0.1451 A | 0.0624 A |
| LC10 cohesive energies | common ten-solid set | 1.5439 eV atom-1 | 1.2993 eV atom-1 |

The g-xTB-specific structures, inputs, results, and validation evidence remain
in this repository. Shared non-GFN reference data are retained where required
to make a Part-I table independently interpretable.
