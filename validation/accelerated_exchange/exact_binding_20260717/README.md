# Exact MPI rank-binding evidence (2026-07-17)

This directory freezes the external qualification harness and final Terok
smoke tests of the exact rank-binding implementation through repository commit
`f4eae80`.  It keeps publication timings separate from historical schema-1
runs that used a shared `taskset` mask and therefore cannot prove which CPU
executed each MPI rank.

## What was verified

- The selected repository integration suite completed **118/118 tests** on the
  rebased tree.  It covers literal disjoint PE lists, forbidden launcher
  overrides, complete removal of inherited `OMPI_MCA_*`/`PRTE_MCA_*` controls,
  single-PU topology checks, cross-process CPU locks, and two Linux `/proc`
  overlap preflights (pool construction and immediately before launch).
- Lock acquisition is fail-safe under `BaseException`: the currently flocked
  handle is registered before JSON metadata, flush, or `fsync`, and every
  post-acquisition constructor step remains inside one cleanup region. Injected
  failures at `json.dump`, `flush`, `fsync`, and the final `threading.Lock`
  retained their tracebacks while the same CPU was immediately reacquired,
  without garbage collection.
- Rank identity and singleton masks are sampled throughout the process
  lifetime.  Sequential same-rank/same-mask PID generations are accepted, but
  rank migration, a successor mask change, an unranked CP2K process, or two
  concurrently live PIDs for one rank fail closed.  Concurrent rank/PID groups
  and their sample indices are persisted.  Record revalidation reconstructs
  every rank, mask, PID-generation, and gate summary from the child-process
  evidence and reparses the hashed Open MPI binding log.
- The standalone Linux harness completed **12/12 tests** on `terok`; the raw
  transcript is `harness/terok_unit_test_transcript.txt`.
- An injected live process named `cp2k.inject`, with no MPI environment and
  allowed mask 224, was rejected before launch on CPU 224.  The raw negative
  transcript is `harness/live_overlap_negative_transcript.txt`.
- The final two-rank Release CP2K `ENERGY_FORCE` smoke ran from
  `2026-07-17T03:02:46Z` to `03:02:50Z`.  Rank 0 remained on CPU 222 and rank 1
  on CPU 223 for all samples.  CP2K returned zero, both reservation/overlap
  gates and the Open MPI binding report passed, and the record is
  `production_scaling_eligible`.
- The debug/LeakSanitizer build with `LSAN_OPTIONS=detect_leaks=0` likewise
  returned zero on CPUs 220 and 221 and is `production_scaling_eligible`.

The unsuppressed LeakSanitizer diagnostic is deliberately retained as negative
evidence.  It returned CP2K code 23 because of known MPI/PMIx teardown leaks.
During sample 99, rank 1 had concurrently live PIDs 3899691 and 3899739,
although both masks remained exactly CPU 221.  The temporal duplicate-rank gate
therefore failed and the record is `timing_non_scaling`.  This is intentional:
only non-overlapping sequential generations may be combined.

## Archived material

- `harness/` contains the fixed 48-run DENSE/STREAMED/QUALIFY matrix driver,
  verifier, fail-closed unit and injected-negative tests, 0D/1D/2D/3D inputs,
  and the allocation/lifetime audit.  `harness/EXACT_BINDING_RERUN.md` is the
  reproducible production protocol.
- `terok_smoke/common/` contains the exact `benchmark_execution.py` snapshot,
  smoke driver, and input used on the host.
- `terok_smoke/release/` and `terok_smoke/pdbg_lsan_suppressed/` contain the two
  positive observations, complete CP2K outputs, complete Open MPI launcher
  reports, driver stdout, and driver return codes.
- `terok_smoke/lsan_diagnostic_non_scaling/` contains the unsuppressed negative
  diagnostic with its persisted concurrent-generation sample.
- `SHA256SUMS` binds every archived file.  Verify from this directory with
  `sha256sum -c SHA256SUMS`.

The CP2K and MPI binaries are not copied here.  Their absolute paths and
launch-time SHA-256 identities are in each `*.execution-smoke.json` record:

- Release CP2K: `44f7d651f88272272dde73b8fb59c9f0c0b37850256fd3ed165c9e1cbcab8129`
- Debug/LSan CP2K: `8cdecb1e925ef1a2e6391da7298bd720d7bf4f335ab8d405127097db9e53033d`
- MPI launcher: `6ffdc5f4649751c7e8d2f6d72ac0ec2c69c9512eeacbe041ffc4f8042cfa00ee`

These two-rank smokes validate the production binding and provenance path; they
are not a substitute for the complete 48-run performance matrix.  Only fully
revalidated schema-2 runs from that matrix may contribute publication speedups.
