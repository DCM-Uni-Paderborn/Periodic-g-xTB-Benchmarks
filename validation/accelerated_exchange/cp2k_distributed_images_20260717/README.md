# Distributed BvK-image exchange qualification

This directory is the immutable qualification bundle for the CP2K
`KGROUP_PARTIAL_DISTRIBUTED_IMAGES` implementation.  It compares the
distributed nonlinear BvK-image kernel with the explicit dense full-mesh
reference for energies, Fock responses, forces, and stresses.

## Revisions

- CP2K branch: `codex/gxtb-partial-distributed-images`
- CP2K commit: `172c141daf` (`Distribute g-XTB exchange over BvK images`)
- save_tblite branch: `codex/distributed-image-kernel`
- save_tblite commit: `56b233288f821c6f30e96dd10e03dd238d9e7e3b`
- CP2K executable SHA-256:
  `28f7b78990ff9675746c016d145b4cdbe1e43fc8e5d4a76e1d7871b843553dce`

## Result

The frozen verifier passed 30 dense/distributed pairs and six deliberate
failure cases.  The matrix covers RKS and UKS, Gamma and shifted meshes,
K290, symmetry reduction, one-, two-, and three-dimensional periodicity,
nondivisor batch sizes, and the `P > N_image` empty-importer case.

- maximum external energy difference: `7.105427357601002e-15 Ha`
- maximum external force difference: `2.1396645149e-15 Ha bohr^-1`
- maximum external stress difference: `3.99999953515362e-06 bar`
- maximum internal forward residual: `2.04281e-14`
- maximum internal reverse residual: `6.66133814775094e-16`
- verified affinity intervals: `66`
- concurrent core conflicts: `0`

`SHA256SUMS` authenticates the complete raw bundle.  Recheck it from this
directory with `shasum -a 256 -c SHA256SUMS`.  `maxima.json`, `summary.tsv`,
and `affinity_concurrency_audit.json` are the concise machine-readable
results; `formal_runs/` contains the complete raw outputs and per-rank
pre-exec affinity proofs.

The single-shot wall times in `summary.tsv` are qualification metadata, not
publication-quality scaling statistics.  Performance claims require repeated
runs with the hardened cross-process CPU-reservation launcher.
