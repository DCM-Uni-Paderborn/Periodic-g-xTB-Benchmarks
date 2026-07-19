# Author-archive structure equivalence

This gate compares only the DMC-ICE13 portion of the externally supplied
author archive with the primitive structures used for the direct save_tblite
and CP2K-native calculations.  It does not read or analyse the unrelated
material contained in that archive.

For every one of the 13 ice structures, the verifier checks the cell matrix,
composition, and a species-resolved periodic atom bijection.  Atomic positions
are compared after wrapping fractional displacements into the minimum image,
so harmless choices of atom order and periodic image do not affect the result.
The archived `report.json` records all source hashes and phase-resolved
residuals.

Reproduce the comparison with:

```bash
python3 compare_structures.py \
  /path/to/gamma_only_dmc_ice13_x23b_exchange_20260606.zip \
  ../../../seidler_dmc13_recalculation/structures/primitive \
  --output-json report.json
```

Passing this gate excludes a cell, coordinate-unit, periodic-image, atom-order,
or structure-selection difference between the supplied DMC-ICE13 structures
and the production inputs.  It does not establish which save_tblite source
history produced any previously published energy.
