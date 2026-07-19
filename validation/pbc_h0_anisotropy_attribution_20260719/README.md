# Attribution of the final `pbc` provider-energy difference

The diagnostic calculation starts from `cp2k-integration` source
`15915c9435644eb257178ca8f8bf7220c38b1a84` and changes only two historical
`pbc` behaviours:

1. periodic g-xTB H0 anisotropy is evaluated from atom pairs inside the explicit
   cell rather than from the periodic H0 neighbour list; and
2. the final-`pbc` equal Wigner--Seitz weights are used.  A preceding A/B test
   showed that this weight change moves phase VII by less than `1e-12` hartree.

For the identical DMC-ICE13 phase-VII `2 x 2 x 2` explicit Born--von Karman
supercell the diagnostic energy is `-7352.468520192105` hartree.  The independently
built final-`pbc` executable gives `-7352.468520192124` hartree, a difference of
only `1.91e-11` hartree.  The unchanged current provider gives
`-7352.465349096680` hartree.  Consequently the historical H0-anisotropy treatment
accounts for more than 99.999999% of this provider-energy gap.

This is a source attribution, not evidence that the historical treatment is
preferable.  It omits periodic-image pairs crossing the explicit cell boundary;
the current integration uses the periodic neighbour list and is intended to be
translation- and cell-partition invariant.

The `invariance` calculation represents exactly the same periodic crystal twice:
in the second POSCAR one hydrogen atom is displaced by a complete lattice vector.
With periodic-neighbour H0 anisotropy the two energies differ by only
`2.73e-12` hartree per supercell.  With historical central-cell anisotropy they
differ by `5.60e-8` hartree.  The latter is small in chemical units but is a real
representation dependence, more than four orders of magnitude above the current
path's numerical residual.  This supports retaining the image-complete periodic
H0 treatment and reporting final benchmark data from that definition.

The diagnostic executable SHA-256 is
`5380387260fc6ec6968d8e430c866ff39004832875a8bcfc46bb63b9680f5874`.
Run `python3 verify_h0_attribution.py` to verify the raw result and regenerate
`verification.json`.
