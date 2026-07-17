# Terminal `/proc/environ` race evidence (2026-07-17)

This directory preserves the unmodified files that exposed a Linux procfs
teardown race in the exact-rank qualification harness.  They were copied from

`/home/kuehne88/work/gxtb-exact-star-matrix-fullbuild-20260717/harness`

on `terok` after the campaign had stopped.  No CP2K calculation was restarted
and no file in the original result tree was changed.

## Observation

The `o2_uks_tr_p4_dense` launcher reported the complete binding

```text
rank 0 -> CPU 196
rank 1 -> CPU 197
rank 2 -> CPU 198
rank 3 -> CPU 199
```

and CP2K printed one normal `PROGRAM ENDED` marker.  During the final monitor
sample, PID 3906903 still had state `R` and singleton mask `197`, but its
environment read no longer yielded `OMPI_COMM_WORLD_RANK`.  The old monitor
therefore changed its previously proven rank 1 to `None` and failed closed.
The initial `run.json` consequently remained
`timing_pending_full_revalidation`; the run was never accepted as timing
evidence.

This is not evidence of CPU overlap.  A read-only audit after the failure found
no live CP2K processes and no simultaneously running pair among the 35 created
matrix run directories with intersecting ordered PE lists.  Repeated CPU
numbers in process displays were sequential slot reuse or stale `PSR` values
from unrelated zombie processes.

## Hardened rule

An empty or unreadable environment is now retained only as a provisional
terminal suffix after a rank was explicitly proven.  The PID, Linux procfs
`starttime`, and exact singleton mask must remain unchanged, and there must be
no earlier migration, mask, or duplicate-rank violation.  The provisional rank
still participates in duplicate detection and its PID is monitored directly,
even if the descendant list temporarily omits it.

The suffix becomes valid only when the same `(PID, starttime)` reaches `Z`/`X`
or disappears and the complete Open MPI binding report is present.  A live
child after launcher exit, PID reuse, a changed mask, a nonempty environment
with missing/invalid rank, or any later explicit-rank observation remains a
sticky failure.  Raw observations and terminal resolution are persisted and
independently revalidated.

In addition, every runtime sample performs a hostwide overlap scan.  Exclusions
are immutable `(PID,starttime)` identities rather than raw PIDs and are checked
again after procfs sampling, so PID reuse cannot hide a foreign owner.  Exactly
one PID generation is accepted per rank.  Every post-launch `BaseException`
drains the launcher group and any directly tracked escaped rank before the CPU
reservation can be returned.

## Original identities

| File | Original size / modification time (CEST) | SHA-256 |
|---|---:|---|
| `campaign.log` | 4923 B / 2026-07-17 05:15:06.108 | `cf23a8fc653f0d7bfafed69055ec562394ec46c013478ac180c8a614bbbc8218` |
| `o2_uks_tr_p4_dense.run.json` | 2670 B / 2026-07-17 05:14:58.276 | `a9226af7455af056074992ea47e43797dfd3dbb99ca434e69275a0b284653c58` |
| `o2_uks_tr_p4_dense.launcher.log` | 216 B / 2026-07-17 05:14:58.440 | `b372c621b7ded3a0c198c87c8e14807c67828e5ecd16467d94a60ab2812ad583` |
| `o2_uks_tr_p4_dense.cp2k.out` | 40193 B / 2026-07-17 05:15:01.372 | `80cd40fedb0c92bb96d667477fbc0a136d2dd2b454c0b798601ebdf728ab260c` |

The repository-level `../SHA256SUMS` also binds these copied files and this
explanation.
