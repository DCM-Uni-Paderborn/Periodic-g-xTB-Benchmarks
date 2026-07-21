# Current-build low-k and derivative rerun

This package repeats all 23 inexpensive periodic g-XTB low-k and derivative
qualification cases with the exact CP2K executable used for the accepted
DMC-ICE13 series.  Its SHA-256 is
`b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`.

The campaign contains 11 three-dimensional derivative calculations and 12
one- and two-dimensional periodicity diagnostics.  Every calculation ran
sequentially on the disjoint singleton CPU 136, with OpenMP and all recorded
BLAS thread counts fixed to one.  The launcher checked the reservation
registry, live affinities, available memory, executable identity, and input
identity before execution.  All 23 outputs returned zero, terminate normally,
and pass the independent verifier.

The rerun reproduces the archived qualified numerical summary exactly:

- maximum three-dimensional finite-difference residuals are `1e-8`
  hartree/bohr for forces, `3.454e-8` hartree for the virial sum, and
  `6.189449e-7` GPa for the converted stress sum;
- the corresponding partial-periodicity maxima are `2e-8` hartree/bohr,
  `1.25643e-7` hartree, and `1.099905e-6` GPa;
- full, K290, and SPGLIB energies agree at printed precision for the tested
  one-dimensional, two-dimensional, and cubic `2 x 2 x 2` routes;
- the H2 Gamma-supercell result differs from the native full `3 x 1 x 1`
  mesh by `-5.472667e-11` hartree per primitive cell;
- implicit and explicit Gamma H2O differ by `3.666401e-12` hartree;
- the one-dimensional Born--von Karman residual is zero at printed precision,
  while the two-dimensional diagnostic residual is `4.069989e-9` hartree.

`verification.json` is the authoritative machine-readable result.  The
successful controller ended by writing `STATUS` with `status=PASS`, after
logging completion of the 23rd case.  The controller script did not itself
emit the verifier's legacy
`controller.exit_status` sidecar, so the archived sidecar records zero as
deduced from that terminal PASS record; all individual scientific exit codes
remain independently preserved and checked.

Reproduce the decision with:

```bash
python3 verify_final_lowk_derivatives.py . \
  --legacy-partial-manifest legacy_partial_pbc_manifest.json \
  --output verification.reproduced.json
```
