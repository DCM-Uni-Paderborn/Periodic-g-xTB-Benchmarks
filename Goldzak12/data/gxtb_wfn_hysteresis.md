# LC12 g-xTB wavefunction-continuation diagnostic

## Scope and protocol

This diagnostic tests whether the anomalous LiH and MgO equation-of-state curves arise from multiple self-consistent roots. It is deliberately excluded from the production LC12 data. Each target calculation uses the converged CP2K Bloch density/orbital restart (`*-RESTART.kp`) from a specified source volume. All Hamiltonian and numerical settings otherwise remain those of the frozen Terok campaign: `METHOD GXTB`, `ACCURACY 0.05`, native save_tblite Fock DIIS, shifted MACDONALD k444 mesh, SPGLIB reduction, 300 K electronic temperature, and `EPS_SCF 1e-9`.

Only the CP2K Bloch density/orbital guess is transferred. The internal save_tblite charge, dipole, quadrupole, and FDIIS-history state is not restartable in the frozen build and is reconstructed or reinitialized. The frozen CP2K build prints an explicit `WFN_RESTART` line only when `SCF/PRINT/RESTART LOG_PRINT_KEY` is enabled. The unchanged inputs did not enable that logging option. Every successful target therefore records the independent evidence tuple

- `Density guess: RESTART` reported by CP2K,
- normal `PROGRAM ENDED` completion,
- no missing-file/atomic-guess fallback,
- no atom-count-mismatch fallback, and
- no restart-read error.

The campaign identity is CP2K revision `18d37c946413dba1b848f57563c46d16b866ce20`, CP2K binary SHA256 `c6b51be7e356170dcb39a597d0e389bd701586e6131365ba317da3968c36eea7`, save_tblite revision `1449febde312874cd0fac4227919f5ba4e4b69b8`, and build-manifest SHA256 `b6ba8a9aeeebc5feca2f42c212bb79fcd19c95162d6f11963ca57711fa8d663f`. The diagnostic driver SHA256 is `cd64afb89c7f8c80d12d1ac0801bea6a887ff4dfc8caec02fc2bf75c5f54bf21`.

## Result

Eleven of thirteen directed continuations completed. The four regular MgO paths between scales 0.80/0.82 and 0.94/0.96/0.98 reproduce the independent target roots within `2.3e-13` hartree. Both paths that require a calculation at MgO scale 0.85 fail the complete-mesh Fock symmetry invariant: the independent 0.85 seed fails at FDIIS step 27 with residual `2.0771e-4`, while continuation from 0.82 fails at target step 21 with residual `1.1204e-4`. Restarting from the nearby regular root therefore does not cure the MgO instability.

LiH exhibits decisive root hysteresis. Continuation from scale 1.375 to 1.300 yields `-32.77037359901141` hartree, which is `-0.31231396744967` hartree below the independently initialized 1.300 result. Conversely, continuation from 1.000 to 1.375 yields `-32.43710154175437` hartree, which is `+0.37951192868372` hartree above the independently initialized 1.375 result. The other three LiH paths reproduce the corresponding independent targets to within `5e-12` hartree.

The distinct LiH roots also differ qualitatively. At scale 1.300 the independently initialized high root has Fermi energy `-0.30720778` hartree and CP2K Mulliken net charges of about `-0.218/+0.218` on Li/H, whereas the root continued from 1.375 has Fermi energy `-0.64936174` hartree and charges of about `-1.086/+1.086`. At scale 1.375 the root continued from 1.000 retains weak charges of about `-0.185/+0.185`, while the independently initialized low root is strongly reorganized. This is a self-consistent branch/charge-state problem rather than equation-of-state fit noise. The energetically lower branch must not be accepted automatically: its strongly charge-inverted Mulliken state requires model-level validation.

Consequently, neither LiH nor MgO is eligible for a final g-xTB LC12 lattice constant or cohesive energy in the frozen campaign.  For the remaining common ten solids, the provisional lattice-constant statistics recomputed from `eos_fits.csv` (SHA256 `e2ca4d372052d3bed4938dc09961223a4b5a338ae2ef27465101c18590bcc868`) are ME `-0.16331768679` angstrom, MAE `0.16384382993` angstrom, RMSE `0.18858340070` angstrom, and MaxAE `0.28257364610` angstrom.  On the identical ten-system subset the GFN1-xTB and GFN2-xTB MAEs are `0.14511781378` and `0.06240968509` angstrom, respectively.  These g-xTB values remain a pre-#5582 diagnostic snapshot, not production results.  No g-xTB LC12 cohesive-energy MAE is reported.
