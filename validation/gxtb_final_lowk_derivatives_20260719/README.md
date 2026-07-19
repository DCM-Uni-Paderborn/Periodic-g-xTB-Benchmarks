# Final-build low-k and derivative qualification

This package repeats the inexpensive periodic g-xTB implementation gates with
the final clean CP2K/save_tblite build used for the DMC-ICE13 same-build
qualification.  It contains 23 independent calculations:

- 11 three-dimensional derivative tests covering implicit and explicit
  Gamma, a full and time-reversal-reduced `3 x 1 x 1` mesh, a commensurate
  Gamma supercell, and full/K290/SPGLIB `2 x 2 x 2` routes;
- 12 one- and two-dimensional tests covering shifted and Gamma-centered
  meshes, commensurate supercells, full meshes, K290, and SPGLIB.

Every calculation was run as a single process on an exactly pinned logical
CPU with all OpenMP and BLAS thread counts fixed to one.  The launcher checked
the shared reservation registry and live-process affinities before execution.
All 23 calculations returned zero, contain a normal CP2K termination marker,
and reproduce the recorded executable and input hashes.

`verification.json` is the authoritative machine-readable result.  The
finite-difference maxima are

- `1e-8` hartree/bohr for the three-dimensional force checks and `2e-8`
  hartree/bohr for the partial-periodicity force checks;
- `3.454e-8` hartree and `1.25643e-7` hartree for the corresponding summed
  virial differences.

At printed precision, full, K290, and SPGLIB energies are identical for the
one-dimensional, two-dimensional, and cubic `2 x 2 x 2` tests.  The full and
time-reversal-reduced `3 x 1 x 1` energies are also identical.  The H2 Gamma
supercell differs from the native full mesh by only
`5.473e-11` hartree per primitive cell.

Most importantly, this final build closes the previously open partial-PBC
Born--von Karman supercell residual.  The one-dimensional supercell and
Gamma-centered native result are identical at printed precision.  The
two-dimensional residual is `4.070e-9` hartree, or
`1.069e-5` kJ/mol per primitive cell.  Relative to the archived earlier
qualification, the two-dimensional absolute residual is smaller by a factor
of about `1.16e4`; the one-dimensional residual decreases to zero at printed
precision.  The archived earlier manifest is retained as
`legacy_partial_pbc_manifest.json` and its hash is embedded in the new
verification summary.

The first controller attempt stopped before launching any scientific process
because its result-parent directories had not yet been created.  This
non-scientific setup failure is preserved in
`controller.initial_setup_failure.log`; the corrected persistent controller
and all scientific outputs are retained separately.

Run the independent verifier with

```bash
python3 verify_final_lowk_derivatives.py . \
  --legacy-partial-manifest legacy_partial_pbc_manifest.json \
  --output verification.reproduced.json
```

`SHA256SUMS` covers every retained input, output, affinity proof, controller
record, and verification artifact present when the package was frozen.
