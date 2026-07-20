# Wigner--Seitz self-image attribution

This directory isolates the dominant difference between the historical
`mstore-inorganic` and later `pbc` source states of `save_tblite`.  All energy
comparisons use the same archived DMC-ICE13 Ih and VII Cartesian structures,
the same explicit 2x2x2 Born--von Karman construction, 96 water molecules per
supercell, `--acc 0.1`, and independently reconverged SCC solutions.

The source change under test is the `get_pairs` fix in
`src/tblite/wignerseitz.f90`.  After the zero-distance origin is removed, the
distance array is compacted.  The historical routine returned the position in
that compacted array as though it were the original translation index.  For a
periodic self-image, compacted entry 1 therefore pointed back to translation 1,
the zero vector, instead of the actual nearest periodic image.  Downstream
exchange code either skipped or misassigned that image.  The correction keeps
an `orig` map and returns `orig(pos)`.

Two reciprocal controlled tests were performed:

1. The exact `pbc` source/dependency build was evaluated with the old routine
   and then rebuilt with only the corrected routine.
2. The exact historical `mstore-inorganic` source and pinned dependencies were
   rebuilt with only the corrected routine.

The corrected same-build `pbc` executable reproduces the archived author-`pbc`
Ih and VII energies to 2e-12 hartree or better.  The old-to-correct WSC shifts
from the two reciprocal builds agree within 0.5 kJ mol-1 per H2O.  More than
95% of the full historical branch gap is therefore causally assigned to this
single indexing correction.  The remaining few kJ mol-1 per H2O contain the
other source and dependency differences and are deliberately left
unattributed.

This result also changes the interpretation of the apparently better sparse
k-mesh behavior of the historical branch: most of that improvement resulted
from missing or misassigned periodic self-image exchange, not from a more
accurate periodic model and not from an error in the CP2K-native interface.

Run `python3 evaluate_wsc_attribution.py` to regenerate
`wsc_relative_energies.csv` and `verification.json`.  The verifier checks every
input and executable record, SCC thresholds, exit state, the reciprocal shifts,
and agreement with the archived author-`pbc` energies.  Raw executable files
are retained in the external audit workspace; the repository copy records
their SHA-256 values without versioning the binaries.
