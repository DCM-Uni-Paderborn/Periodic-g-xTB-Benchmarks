# Partial k-to-R provider validation manifest

Archived: 2026-07-16

Provider repository: `DCM-Uni-Paderborn/save_tblite`

Provider branch: `codex/partial-k-to-r-accumulator`

Provider base: `257ba442684c39454175e5192c8a2342b4c6380f`

Qualified provider commit: `e57095e15e88a3cf5cecb70cf1379f68ff05867c`

Benchmark archive base: `055353ed1281b0da14e5439cecc212eb01e09d42`
(`origin/main` when this archive branch was created)

## Scope

This directory is a self-contained qualification record for the MPI-free,
bounded, exact partial k-to-R accumulator provider ABI.  It records the source
snapshot and patch, macOS Debug/Release unit and CTest logs, controlled Linux
Debug/Release logs from terok, the controlled two-rank smoke test, and an MPI
link/source audit.  The provider contract and numerical test matrix are
described in `README.md`.

The controlled terok qualifications used one OpenMP/BLAS thread per process.
Unrestricted jobs already running on terok were left untouched and are not
used as qualification evidence.

## Files

- `README.md`: ABI contract, payload formulae, test matrix, and results.
- `branch-status.txt`: qualified provider branch and full commit.
- `source-worktree.tar.gz`: `git archive` of the qualified provider commit.
- `worktree.diff.gz`: deterministic gzip (`-n`) of the binary-safe diff from
  the provider base to the qualified commit.
- `source-files.sha256`: hashes of the nine modified provider source/test
  files, with paths relative to the provider repository root.
- `local-*-exchange.log`: controlled macOS Debug/Release exchange tests.
- `local-*-ctest-exchange.log`: macOS CTest exchange targets.
- `terok-linux-*-exchange.log`: controlled Linux Debug/Release tests.
- `terok-linux-mpi2-release-smoke.log`: controlled two-rank Release smoke.
- `terok-linux-mpi-freedom.log`: provider MPI symbol/link/source audit.
- `SHA256SUMS`: SHA-256 digest of every other file in this directory.

## Integrity check

From this directory, run:

```sh
shasum -a 256 -c SHA256SUMS
```
