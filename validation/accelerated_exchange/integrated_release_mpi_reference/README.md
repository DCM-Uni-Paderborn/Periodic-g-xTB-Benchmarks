# Integrated CP2K Release/MPI reference qualification

This directory freezes the focused end-to-end qualification of the integrated periodic g-xTB
acceleration tree.  It combines the streamed reverse consumer, owner-generated CP2K k-point-group
bridge, separable regular-mesh transform, and validated cross-mesh restart work on the same CP2K
source tree after merging the then-current upstream master.

The scientific input is a two-atom silicon cell on a shifted `2 x 2 x 2` MacDonald mesh with
`FULL_GRID ON`, one MPI rank per CP2K k-point group, analytical forces, and analytical stress.  It
therefore has `nred = nfull = 8`, while the `P = 1, 2, 4` KGROUP runs use exactly `1, 2, 4` official
CP2K k-point groups.  The dense complete-mesh path is retained as the numerical reference.

## Test design

The warm-up/reference archive contains 15 successful Release runs:

- `P = 1, 2, 4` for DENSE, BATCHED with streamed reverse, KGROUP_OWNER with streamed reverse, and
  the separable transform with dense derivatives;
- one KGROUP_OWNER/QUALIFY run at each MPI size, including the dense forward and reverse oracles.

The timing archive contains 36 additional successful Release runs.  At each MPI size, all four
production modes were repeated three times in alternating orders.  The preceding 15-run matrix is
kept separate as warm-up and is not included in the timing statistics.  Every run used the same
CPU set for a given MPI size, one OpenMP thread per MPI rank, one OpenBLAS thread per MPI rank, and
the same executable, input, environment, and process-tree RSS sampler.

All 51 positive runs returned zero and printed `PROGRAM ENDED`.  The KGROUP qualification oracle
gave at most `3.219647e-15` hartree in energy and `5.551115e-17` in the folded Fock response, with
zero shell-potential difference.  The streamed reverse oracle gave at most `1.110223e-16` in the
overlap adjoint and zero direct-force and direct-stress differences for this cell.

Across all 36 measured production runs, the printed energy was uniquely
`-579.0508385538171` hartree.  Relative to the same-rank dense reference, the largest force-component
difference was `2.40355431e-16` hartree/bohr.  The largest stress-component difference was
`6.63103935e-10` bar and occurred only in an off-diagonal component at roundoff; every diagonal
stress component was identical at printed precision.  See `invariance_summary.tsv`.

`timing_summary.tsv` reports the median, minimum, maximum, and median absolute deviation from the
three measured repetitions.  Whole-process RSS differences are negligible for this small case,
because replicated CP2K state dominates.  BATCHED has visible transaction overhead here, whereas
KGROUP_OWNER and the separable transform are within the observed end-to-end variability of DENSE.
In particular, the relatively broad `P = 4` ranges preclude a speedup claim.  The provider kernel
also remains single-state in this bridge until a public mergeable partial-k-to-R accumulator ABI is
available.  These data establish numerical equivalence and selectable-path robustness; they do not
claim acceleration or distributed provider-kernel scaling.

## Raw data and verification

The two archives under `raw/` contain every input, output, selected diagnostic, return code,
timing/RSS JSON file, per-run SHA-256 manifest, runner script, campaign script, and top-level
manifest.  Their expected hashes are recorded in this directory's `SHA256SUMS`.

From this directory, verify the outer archive hashes with:

```sh
sha256sum -c SHA256SUMS
```

To verify an archive recursively, extract it into an empty directory and run `sha256sum -c` on its
top-level `*_SHA256SUMS` file.  Both archives were verified this way before publication.
