# Repository and branch scope

The public repository is split by paper responsibility:

- `main` is the Part-I branch. It owns the periodic g-xTB application
  benchmarks, integration validation, paper-facing tables/figures, and the
  retained technical molecular-crystal development archive.
- `part-II` is the acceleration branch. Its tip contains only reference-oracle
  comparisons, force/stress validation, timing/memory data, scaling tests, and
  exact provenance for accelerated Brillouin-zone-coupled nonlocal exchange.

Part-II evidence is self-contained. Inputs mentioned through historical paths
inside immutable manifests are also stored in the corresponding campaign or
validation archive; no Part-I benchmark directory is required to reproduce a
Part-II check.

Complete GFN1-xTB and GFN2-xTB method-owned benchmark inputs and results are
canonical in
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
They are not duplicated on `part-II`.

Removing Part-I datasets from this branch is lossless: their current versions
remain on `main`, and all earlier branch snapshots remain reachable through
Git history.
