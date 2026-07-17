# Canonical GFN1/GFN2 benchmark source

Complete periodic GFN1-xTB and GFN2-xTB inputs, raw/curated results, figures,
and the corresponding baseline execution scripts live in
[`DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks).
The source revision used for the repository split is
[`279413c7d7ad8b0b553a3e423293ad581f3960bf`](https://github.com/DCM-Uni-Paderborn/Periodic-GFN2-Benchmarks/commit/279413c7d7ad8b0b553a3e423293ad581f3960bf).

On 2026-07-17, 512 method-owned GFN1/GFN2 files were removed from the current
tips of this repository. Of these, 509 paths were verified to have the same
Git blob at the same path in the source revision above. The inventory is:

- 205 DMC-ICE13 GFN inputs, derived data/figures, and baseline-only scripts;
- 301 X23b GFN inputs, derived data/figures, and its baseline-only launcher;
- one byte-identical LC snapshot table; and
- two byte-identical top-level X23b baseline utilities.

The remaining three files were a GFN-only result digest, its numeric index,
and a GFN-only X23b runner. Their scientific data already exist in the
canonical repository; the runner's additional method guard was first migrated
to the canonical copy.

The removal is lossless: the data remain in the canonical GFN repository and
all former file versions remain in this repository's earlier Git history.
Shared benchmark geometries,
experimental/DFT reference data, and files containing g-xTB results remain
here. Compact GFN values retained inside explicit g-xTB-versus-GFN comparison
tables are comparison metadata, not a second raw-data source.

Some historical LC12 snapshot files are intentionally retained because their
full contents are not present in the cited GFN source revision. They are not
current g-xTB results and must not be used as the canonical periodic GFN paper
dataset.
