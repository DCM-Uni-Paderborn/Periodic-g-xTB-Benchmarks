# Exact-build automatic-policy regression

This archive closes the release-provenance gate for the default periodic
g-xTB production policy.  The complete focused CP2K regression directory was
run after both source branches had been published and after CP2K had been
rebuilt from the final public revisions.  The executable itself reports the
accepted CP2K revision.

## Accepted result

- regression directory: `xTB/regtest-tblite-gxtb`
- input calculations: 48
- matcher results: 196 correct, 0 wrong, 0 failed
- normally terminated outputs: 48 of 48
- CP2K source: `g-xTB-pbc` at
  `f44008823d3319547f34ef335561256816a1a031`
- save_tblite source: `cp2k-integration` at
  `718629fbc86e0b362491cf70dd4198d0d82082b5`
- CP2K executable SHA-256:
  `8850e8a39c14fbd172ab89a7992cee69e492b3d6ab039451985c312711e3e0aa`
- CP2K shared-library SHA-256:
  `efaa4fcbbb36a00e13e630f05f14a9a4b87c8e04b12b7f6ee9fa9d54630694b0`

The run used two MPI ranks, one OpenMP thread, and one thread for every
recorded BLAS runtime.  It is a correctness regression, not a timing or
scaling campaign.  The separate SPGLIB-only directory was not selected by
this build because SPGLIB was absent from its compile flags; the archived
Part-II SPGLIB qualification remains the evidence for that optional backend.

The default CH4 K290 calculation printed
`mode=AUTO exchange=SYMMETRY_FUSED gradient=STREAMED
transform=MIXED_RADIX_FFT`, selected the streamed ACP forward contraction and
sparse ACP reverse, terminated normally, and reached
`-40.468866070692435` Ha.  The molecular open-shell O2 control reproduced the
save_tblite CLI energy to `5.911715561524e-12` Ha; its maximum gradient
difference was `3.300347844518e-06` Ha/bohr.

## Molecular-limit source history

The actual molecular-limit correction is provider-side: save_tblite commit
`cbec298` fixes the molecular quadrupole large-cell limit, and the final
`718629f` commit refreshes the affected CLI references.  CP2K has no separate
molecular-limit correction in this source interval.  Its final commits add
and qualify exact periodic ACP streaming and make that path the automatic
multi-point default.  The selected source histories are retained in
`source_history.txt`.

## Reproduction and retained files

`command.txt` gives the exact focused-regression invocation.  The two compact
input/output pairs expose the default AUTO markers and the molecular CLI
control directly.  `raw/focused_gxtb_regression_raw.tar.gz` contains the
complete 3.7-MiB working directory for all 48 calculations, including inputs,
outputs, restart files, and the matcher manifest.

Run `python3 verify.py` from any directory to verify every retained SHA-256,
the source and binary identities, all 48 normal terminations, the 196 matcher
definitions, the automatic ACP selectors, and the two reported numerical
results.
