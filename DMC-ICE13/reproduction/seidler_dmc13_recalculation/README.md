# DMC-ICE13 recalculation package for the save_tblite authors

This compact package is intended for an independent rerun with any selected
`save_tblite` executable.  It contains all 13 DMC-ICE13 cells in absolute
Cartesian Angstrom coordinates, the published Ih-referenced DMC energies, and
the absolute energies obtained with the final author `pbc` provider, the
current CP2K-integration provider, and CP2K-native Bloch k points.

The primitive structures live under `structures/primitive`.  Explicit
Gamma-centred Born--von Karman cells are generated without changing atom
order by `scripts/build_bvk_from_poscar.py`; its output is checked against the
archived `1 x 1 x 1` through `4 x 4 x 4` cells with exact species/order and a
maximum cell/coordinate tolerance of `5e-12` Angstrom.  This avoids
shipping additional multi-gigabyte dense-mesh copies while retaining a fully
deterministic structure definition.

Run one phase with a chosen executable, for example:

```text
python3 scripts/run_save_tblite.py /path/to/tblite VII 3 results/my-build \
  --require-binary-sha256 <sha256>
```

The command is exactly equivalent to:

```text
tblite run --method gxtb --acc 0.1 --iterations 300 --no-restart \
  --json result.json POSCAR
```

`tables/current_absolute_energies_by_mesh.csv` contains every available
absolute current-CLI and CP2K-native value from `1^3` through `4^3`.
`tables/author_pbc_absolute_energies.csv` contains the complete final author
`pbc` direct-CLI series at `2^3` and `3^3`.  The companion relative table and
the direct three-route `3^3` closure distinguish provider-model changes from
the much smaller CP2K-native integration residual.

All generated files are covered by `SHA256SUMS`.  Rebuild and verify the
package from its authoritative parent archive with:

```text
python3 prepare_package.py
shasum -a 256 -c SHA256SUMS
```
