# CP2K bounded image-batch runtime qualification

This immutable campaign archives the completed local CP2K runtime checks for
the caller-driven, bounded image-batch g-xTB exchange transaction.  The final
run ended at `2026-07-16T16:04:44Z`, which defines the campaign identifier.

## Scope

The archive contains 34 completed outputs and all ten distinct CP2K inputs
used by them.  Every archived output contains `PROGRAM ENDED`.  Coverage
includes:

- complete `2x2x2` RKS meshes with legacy evaluation and image batch sizes
  1, 2, 7, 8, and 9;
- the same batch-size-one case with two MPI ranks and with two OpenMP threads;
- SPGLIB-reduced `2x2x2` energy/force runs with legacy evaluation and batch
  sizes 1, 3, and 8;
- SPGLIB and K290 analytical force/stress DEBUG checks with legacy evaluation
  and batch sizes 1 and 8;
- 1D and 2D force/stress DEBUG checks for full and SPGLIB-reduced meshes;
- UKS `2x2x2` runs with legacy evaluation and batch sizes 1 and 8;
- the `6x6x6` (`K=216`) timing curve for legacy evaluation and batch sizes
  8, 27, 72, and 216;
- both the first post-fix Debug pair and the final operation-preserving
  minimal-zero-safe legacy/batch-size-one oracle pair.

`run_manifest.tsv` maps every archived file to its original path and records
the intended mode, batch size, symmetry path, dimensionality, spin treatment,
and build class.  `derived/summary.json` is the full machine-readable record;
`derived/summary.tsv` is its compact per-run projection.

## Final operation-preserving zero-safe regression

The authoritative final Debug pair is
`debug_zero_safe_minimal_legacy.out` and
`debug_zero_safe_minimal_b1_oracle.out`.  Both return normally at the CP2K
level and print the same total energy,
`-40.473758159181223` hartree.  The largest printed stream-oracle residuals in
the batched run are `1.110223e-16` for energy, `5.551115e-17` for the folded
Fock response and covariance response, and `1.396274e-16` for the dual
response.  The raw output remains the authority; these values are parsed again
by `scripts/summarize_runtime.py`.

The provider change writes the harmonic mean in a zero-safe form: a zero
angular-momentum-channel parameter returns the mathematically defined limiting
value zero, while the original arithmetic expression is retained verbatim for
every finite nonzero pair.  The earlier post-fix Debug pair is retained as
historical evidence; the `minimal` pair is the final qualification basis.

## K=216 memory/work trade-off

The `6x6x6` case reports CP2K total times of 7.194 s for the legacy full-image
path and 90.113, 30.898, 15.659, and 9.299 s for batch sizes 8, 27, 72, and
216, respectively.  All printed energies agree within
`7.105427357601002e-15` hartree.  This is a memory/transform-recomputation
trade-off, not an acceleration claim: small batches bound image storage but
re-submit and transform all 216 k blocks many times.  At batch size 216 the
batched path approaches full-image storage and still remains slower than the
legacy path in this small molecular test.

## Provenance boundary

`provenance/source_state.txt` records the CP2K and save_tblite base revisions,
working-tree diff hashes, final post-fix Debug binary/library hashes, mtimes,
and build paths.  Both CMake caches are preserved verbatim.

The original RelWithDebInfo launcher/library were rebuilt after the earlier
qualification series completed.  The overwritten executable cannot be hashed
retroactively, so no current binary hash is falsely assigned to those runs.
Their unedited inputs and outputs, CP2K base revision, source worktree identity,
and configuration are retained.  The final post-fix Debug pair does have exact
binary, linked-library, provider-library, source-state, and configuration
provenance.

The `b1`, `b2`, and similar labels record the batch size selected for the run.
The original ad-hoc launch commands were not separately logged, which is an
explicit limitation of the pre-fix series.  CP2K independently prints MPI
rank count, OpenMP thread count, version, input, completion, timing, energies,
and qualification residuals into each archived output.

## Reproduction of derived files

Run from this directory:

```text
python3 scripts/summarize_runtime.py
```

The script parses only the archived raw data and manifest.  It regenerates
`derived/summary.json`, `derived/summary.tsv`, and `SHA256SUMS`.  The checksum
manifest covers every campaign file except itself and Python bytecode caches.

No paper or source-code file was modified while creating this archive.  No
running calculation was inspected, stopped, copied while open, or restarted.
