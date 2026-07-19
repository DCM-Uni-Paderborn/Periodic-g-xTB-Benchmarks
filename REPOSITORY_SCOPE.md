# Repository scope

This repository is split by paper:

- `main`: Part I reference implementation, scientific benchmarks, analytical
  derivatives, symmetry validation, PBC tests, and their reproducibility data;
- `part-II`: numerical-method performance and scalability material for Part II.

The `main` tree contains DMC-ICE13, LC10, molecular-limit, force/stress,
symmetry, K290/SPGLIB, lower-dimensional PBC, and Bloch/supercell validation
material. It does not duplicate complete GFN1-xTB or GFN2-xTB datasets; those
remain in `DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`.

The repository originated from the periodic-GFN benchmark collection so that
shared geometries and reference values retain provenance. The current tree is
curated independently for periodic g-xTB Part I. Removed exploratory material
remains recoverable from Git history but is not presented as current paper
data.
