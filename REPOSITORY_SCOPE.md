# Repository scope and provenance

This repository was initialized from revision
`8aa717992a0b8cccc9db3c418df4235b00fedeb0` of
`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks` so that benchmark geometries,
reference data, and the frozen GFN1/GFN2 comparison workflow retain their git
provenance.

The repository boundary is now:

- `Periodic-GFN2-Benchmarks`: canonical GFN1/GFN2 benchmark material only;
- `Periodic-g-xTB-Benchmarks`: all g-xTB inputs, workflows, validation gates,
  provisional and final results, and the frozen GFN1/GFN2 comparison snapshot
  needed for common-subset analyses.

The g-xTB migration was made before any g-xTB commit or branch was pushed to
the public GFN2 repository. That repository therefore requires no history
rewrite; its remote `main` remains the unchanged source revision above.

Raw HPC working directories are not repositories and may temporarily retain
their original directory names. Publication artifacts record content hashes
and build identities rather than treating those transient paths as canonical
provenance.
