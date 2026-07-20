# Independent relative-energy post-processing audit

This gate independently rebuilds the DMC-ICE13 relative energies and error
statistics from raw absolute outputs.  It deliberately does not import the
table-assembly implementation.

The verifier uses decimal arithmetic to check:

- the H$_2$O count derived from each primitive POSCAR;
- the explicit CLI BvK-supercell atom count and division by `N^3`;
- Hartree-to-kJ/mol conversion;
- same-mesh ice-Ih referencing;
- signed and absolute DMC errors; and
- ME, MAE, RMSE, and MaxAE values in the generated tables.

It covers every complete CP2K-native mesh from Gamma through `5 x 5 x 5`
and every complete current-provider CLI mesh.  The direct-CLI matrix is now
complete for all thirteen structures on each mesh from `1 x 1 x 1` through
`4 x 4 x 4`; no partial mesh is allowed to contribute a statistic.

Run from the repository root with

```bash
python3 validation/relative_energy_postprocessing_20260720/verify_relative_energy_postprocessing.py
```
