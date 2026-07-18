# Native periodic g-xTB derivative hardening

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

`raw_diagnostics.tar.gz` preserves all intermediate ablations, unsuccessful
representation trials, displaced structures, inputs, and outputs used to
isolate the issue.  `SHA256SUMS` covers every directly stored artifact and the
diagnostic archive.

The corresponding CP2K implementation is on `g-xTB-pbc`; exact source
revisions are intentionally recorded here in the repository rather than in
the manuscript.

`current_build_gate` repeats the reduced ice-XVII `2 x 2 x 2` energy, force,
and stress calculation after a fresh CMake reconfiguration of the final
source revision.  It reproduces the archived production output exactly in
every printed energy, force, and stress component.  The executable identity,
raw output, and machine-readable zero-difference summary are retained there.
