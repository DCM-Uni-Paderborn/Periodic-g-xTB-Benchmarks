# Repeated ACP timing and whole-process RSS probe

This archive records a same-build, alternating DENSE/STREAMED CP2K comparison
for the bounded atom-centered-potential (ACP) mesh contraction.  It is a scoped
small-system probe, not a scaling study.

## Accepted campaign

The accepted campaign used the production CH4 K290 2x2x2 input, one CP2K rank
on singleton CPU 90, and one thread for OpenMP and every recorded BLAS runtime.
It contains one warm-up per mode followed by five measured repetitions per
mode.  The order alternates within each measured pair.

- CP2K executable SHA-256:
  `df6466552a495e94e710174dbf468ec1765f1a4690ef955f9129ef983f35790b`
- provider archive SHA-256:
  `fe210c64a4c4fa6897668a8657dd234046143c55bda2f7c9279d24108f2f152a`
- input SHA-256:
  `7e1d23a4d9e0d66df81caa3744ca005342eb802599750f6888e8de38d2dd7f4f`
- final energy in every run: `-40.473748967057013` Ha
- maximum paired DENSE/STREAMED force difference:
  `6.6786853e-17` Ha/bohr
- maximum paired DENSE/STREAMED stress difference: `3.6928468e-12` bar

| Mode | measured n | median wall (s) | MAD (s) | wall range (s) | median peak RSS (KiB) | MAD (KiB) | RSS range (KiB) |
|---|---:|---:|---:|---:|---:|---:|---:|
| DENSE | 5 | 21.365591 | 0.080420 | 21.215829--21.446011 | 151508 | 648 | 150336--152844 |
| STREAMED | 5 | 21.238421 | 0.068826 | 21.160578--22.647820 | 151888 | 72 | 151816--152244 |

The STREAMED median is 0.127170 s (0.595%) below DENSE, whereas its sampled
whole-process median RSS is 380 KiB (0.251%) higher.  These small differences
do not support an end-to-end speedup or whole-process RSS-reduction claim for
this case.  The result is compatible with the intended optimization boundary:
the exact provider allocation counters prove removal of the full ACP Bloch
tensor and quadratic projector-difference set, but unrelated CP2K storage
dominates this small calculation.

Before every launch the campaign recorded all live process RSS values and
`MemAvailable`.  It subtracted the remaining-growth allowance of the one
unrelated healthy CP2K job and a 100-GiB candidate peak.  The smallest resulting
margin was 308890012 KiB (294.58 GiB), above the required 128-GiB floor.

## Raw provenance and failed attempts

`raw/gxtb-acp-timing-20260722/` is the byte-for-byte accepted remote tree.
Its original `SHA256SUMS` was generated before `campaign_end` was appended to
`provenance/campaign.txt`; consequently that original manifest has exactly one
expected stale metadata entry.  The raw file is preserved rather than silently
rewritten.  The portable archive-level `SHA256SUMS` binds its final bytes and
all other files.

The three earlier controller attempts are retained under `raw/` but excluded
from statistics:

1. `attempt1` used two nonexistent build paths and launched no calculation.
2. `attempt2` requested unavailable `/usr/bin/time` and launched no accepted
   calculation.
3. `attempt3` exposed a `set -o pipefail` race when process sampling ended.  It
   did not complete the prescribed campaign ledger and is not accepted even
   though its two warm-up CP2K outputs ended normally.

Run `python3 verify.py` from any directory to verify the portable manifest,
the expected raw-manifest ordering defect, all acceptance gates, paired
observables, and the recomputed medians/MADs.
