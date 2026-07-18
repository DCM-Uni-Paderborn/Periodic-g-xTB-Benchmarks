# Native/CLI derivative validation for ice XVII at 2 x 2 x 2

The native calculation uses the primitive ice-XVII cell, a Gamma-centred
`2 x 2 x 2` mesh, and SPGLIB symmetry reduction.  The direct `save_tblite`
calculation uses the explicit 144-atom Born--von-Karman supercell at Gamma.

`native_vs_cli_summary.json` is generated with
`tools/compare_derivatives.py`.  It divides the CLI energy by eight, averages
the eight translational-image gradients for each primitive atom, changes the
gradient sign to obtain forces, and converts the CLI virial to stress using
the explicit POSCAR volume.

The independent central differences displace all eight replicas of the first
primitive atom together for the force check.  The stress check applies a
homogeneous `xx` strain to the complete supercell.  All displaced structures
and raw JSON energies are retained below the two `finite_difference_*`
directories.
