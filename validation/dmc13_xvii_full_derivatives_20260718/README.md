# Complete ice-XVII derivative qualification

This package records the force/stress diagnosis and final validation of the
native CP2K periodic g-xTB response.  The production implementation uses two
algebraically exact ACP reverse representations:

- time-reversal-invariant meshes are contracted in CP2K's real-space image
  convention, with zero density outside the overlap support;
- genuinely complex meshes use the direct Bloch-space ACP reverse sweep.

The split avoids a storage-convention ambiguity at self-inverse k points
without changing the forward energy.  No approximate truncation or numerical
screening is introduced by this choice.

The ice-XVII `2 x 2 x 2` full-grid and symmetry-reduced calculations have
identical printed energies, forces, and stresses.  The reduced native result
also agrees with the explicit 144-atom save_tblite CLI supercell to the values
reported in `summary.json`.  Independent central differences validate atom-1
`x` force and homogeneous `xx` strain.  The complex `3 x 1 x 1` H2 test and
the self-inverse `2 x 2 x 2` CH4 tests exercise both ACP reverse paths.

`cli_supercell_validation/` retains the independent direct-CLI gradient,
collective-displacement finite differences, strained structures, and absolute
energies.  `raw_diagnostics.tar.gz` preserves the native full/reduced outputs,
native strain differences, regression probes, intermediate ablations, and the
representation trials used to isolate the response error.  The historical
`SHA256SUMS` entries for compacted outputs are resolved to their archive
members by `verify_full_derivatives.py`; the verifier then reconstructs every
number reported in the Part-I manuscript and Supporting Information.
`CLI_SHA256SUMS` independently covers every retained direct-CLI input and raw
output.

`run_current_binary_requalification.sh` repeats the complete ice-XVII
full-grid/reduced-grid and central-difference sequence only after the capped
DMC-ICE13 endpoint queue has finished.  All six runs use the same qualified
CP2K executable.  A case that was completed early under the same hash and
singleton-affinity gates is validated and reused rather than run twice.  After
archiving that campaign, run
`verify_current_binary_requalification.py CAMPAIGN_ROOT`; the verifier checks
the executable and input hashes of every run, reconstructs the Cartesian
displacement and homogeneous strain from the inputs, and independently gates
energy, all 54 force components, all nine stress components, and both finite
differences.  The optional `--inputs-only` mode validates the manifest and
perturbation geometry before any production calculation is available.

The corresponding CP2K implementation is on `g-xTB-pbc`; exact source
revisions are intentionally recorded here in the repository rather than in
the manuscript.
