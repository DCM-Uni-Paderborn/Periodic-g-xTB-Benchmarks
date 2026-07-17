# Repository scope and provenance

This repository was initialized from revision
`8aa717992a0b8cccc9db3c418df4235b00fedeb0` of
`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks` so that benchmark geometries and
reference data retain their git provenance.

The repository boundary is now:

- `Periodic-GFN2-Benchmarks`: canonical periodic GFN1/GFN2 paper material,
  including complete method-owned inputs and result datasets;
- `Periodic-g-xTB-Benchmarks`: all g-xTB Part-I/Part-II inputs, workflows,
  validation gates, and provisional/final results. Only compact GFN values
  required by explicit cross-method comparisons remain here.

The verified source revision and byte-identical removal inventory are recorded
in `GFN_BASELINE_SOURCE.md`. Shared structures and non-GFN reference data stay
in both repositories where they are needed to make each benchmark usable.

The g-xTB migration was made before any g-xTB commit or branch was pushed to
the public GFN2 repository. That repository therefore requires no history
rewrite; its remote `main` remains the unchanged source revision above.

Raw HPC working directories are not repositories and may temporarily retain
their original directory names. Publication artifacts record content hashes
and build identities rather than treating those transient paths as canonical
provenance.
