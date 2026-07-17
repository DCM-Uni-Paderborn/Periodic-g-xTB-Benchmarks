# Distributed image-range and separable direct-DFT qualification

This immutable campaign archives two distinct acceleration components for
Brillouin-zone-coupled nonlocal exchange in periodic g-xTB:

1. CP2K distributes disjoint, contiguous Born--von Karman image ranges over
   MPI ranks and adds their partial exchange energy, shell potential, and
   folded Fock response after proving exact-once image coverage.
2. `save_tblite` factorizes the regular-mesh k-to-R and R-to-k transforms into
   direct one-dimensional DFTs along grid pencils.

Both the provider transform and its selectable CP2K/Quickstep consumer are
qualified against the dense transform oracle in this archive.

The second method is a **separable direct DFT, not an FFT**.  Its arithmetic
cost is `O(Nrow * Nk * sum(nmesh))`; no `O(Nk log Nk)` claim is made.  The
dense phase-table transform stays selectable as the numerical oracle and as
the route for unsupported grids.

## Qualification boundary

The distributed image-range forward path is qualified end to end in CP2K.
Nine paired cases cover full, SPGLIB, K290, and time-reversal meshes; RKS and
UKS; shifted and unshifted grids; `P > K` and uneven `K/P`; and true 1D, 2D,
and 3D periodic boundary conditions.  Every one of the 18 production outputs
listed in the summary ended normally.  Across the 3D DEBUG displacement
series, the largest legacy/distributed energy difference is
`7.105427358e-15 Eh`.  Over the full comparison matrix, the largest numeric
force-component difference is `3.330900228e-16 Eh/bohr`, the largest
analytical-stress difference is `1.399999746e-6 bar`, and the largest folded
Fock oracle residual is `3.053113e-16 Eh`.

The provider-level separable transform passes the Debug/fcheck exchange suite
31/31, the g-xTB suite with 40 ordinary passes plus four expected diagnostic
failures, and the Release exchange suite 31/31.  RKS and UKS energy, Fock,
shell-potential, overlap-adjoint, force-response, and stress comparisons use
the dense backend as oracle.  The focused absolute gate is `1e-11`.  The
transform-only timing table preserves both favorable multidimensional grids
and the expected unfavorable long one-dimensional pencils.

The selectable CP2K consumer is qualified by nine tabulated dense/candidate
pairs and 18 normally ended pair entries represented by 17 distinct raw
outputs because the two-rank dense baseline is shared by two comparisons.
The cases cover complete,
SPGLIB, K290, shifted, and time-reversal-reduced meshes; RKS and UKS; true 1D,
2D, and 3D periodicity; analytical energy, force, and stress; a 23-energy
DEBUG trajectory; a two-rank distributed smoke test; and the combined
distributed-forward/separable-derivative selection.  Its largest paired
energy, force-component, and analytical-stress differences are respectively
`7.105427357601002e-15 Ha`, `1.7318587359999999e-16 Ha/bohr`, and
`8.000006346264854e-7 bar`.

A separate release CP2K `6x6x6` CH4 probe contains three sequential dense and
three sequential separable runs, all with `PROGRAM ENDED` and 11
`build_tblite_ks_matrix` calls.  The median exclusive kernel timer decreases
from `0.185` to `0.154 s`, corresponding to `1.20x`.  The median inclusive
timer changes only from `1.131` to `1.103 s` (`1.03x`), while median total
wall times are `7.19` and `7.18 s`.  At the resolution and scale of this
small-system probe, the total wall time is therefore indistinguishable: these
data show reduced transform-kernel work but **no measurable end-to-end
speedup**.

Two additional complete-mesh identity-diagnostic outputs are preserved but
are explicitly **not** production passes.  Both lack `PROGRAM ENDED`: the
dense diagnostic reaches the qualification-only mode-selection ABORT, while
the separable diagnostic reports its dense full-mesh Fock agreement and then
reaches the qualification-only identity-duality ABORT.  Neither appears in a
tabulated production pair.

This archive does **not** qualify a distributed or batch-bounded derivative
algorithm.  The analytical force and stress results agree between the
forward modes, but the reverse/gradient calculation still evaluates the
complete mesh.  It therefore supports a derivative correctness statement,
not derivative memory scaling or MPI speedup.  Likewise, the separable
backend's CP2K consumer is now qualified for correctness, but this does not
turn transform-only timings into whole-SCF timing claims.  Serial image
batching remains a bounded-memory versus recomputation trade-off and is not
presented as a serial speedup.

## Evidence layout

- `raw/distributed_cp2k/` is a byte-for-byte copy of all available CP2K
  inputs, outputs, the original narrative, table, and inner checksum manifest.
  Seven additional outputs without `PROGRAM ENDED` are retained as diagnostic
  or incomplete evidence and are excluded from the qualification rows.
- `raw/provider/` preserves the provider test transcripts, source hashes,
  separable-transform oracle description, benchmark program, and timing TSV.
- `raw/cp2k_separable_consumer/` is a byte-for-byte copy of the CP2K consumer
  inputs, 32 raw outputs, original summaries, and inner checksum manifest.
  Its 30 normally ended outputs include all 17 distinct outputs used in the
  nine tabulated correctness pairs and all six release-timing runs; the two
  remaining identity diagnostics are retained as non-passing qualification
  evidence.  `release_bench/timing_summary.tsv` is the machine-readable
  timing authority.
- `derived/qualification_summary.json` is the authoritative machine-readable
  index.  It records method identities, boundaries, raw case rows, aggregate
  residuals, timing rows, and the completion inventory.
- `provenance/source_state.json` records repository revisions and the exact
  limits of source/binary attribution.
- `SHA256SUMS` covers every campaign file except itself and Python bytecode.

No manuscript, supporting-information source, PDF, or corresponding hash is
stored in this campaign.

## Rebuild and verify the index

From this directory run:

```text
python3 scripts/build_campaign_summary.py
```

The command first verifies both copied CP2K evidence sets against their
original inner manifests, requires every production output referenced by
either table to contain `PROGRAM ENDED`, requires the only non-ended consumer
outputs to be the two named identity diagnostics, regenerates the aggregate
JSON, checks all six release-timing outputs and the exact 3+3 backend matrix,
regenerates the campaign manifest, and verifies that manifest immediately.
