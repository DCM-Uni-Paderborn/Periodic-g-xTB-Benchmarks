# Direct-CLI versus CP2K-native geometry equivalence

This gate compares the primitive cells and atom positions underlying all
thirteen DMC-ICE13 direct-CLI and CP2K-native calculations.  CP2K scaled
coordinates are converted to Cartesian Angstrom and compared with the
canonical absolute-coordinate POSCARs from which the explicit BvK
supercells are generated.

Run `python3 verify_geometry_equivalence.py` to regenerate
`verification.json`.  Cells are identical and the largest Cartesian
coordinate difference is at binary floating-point roundoff.  Structure
conversion is therefore excluded as the source of the small remaining
energy-path residual.
