# Request for an independent DMC-ICE13 branch comparison

Dear Leopold,

we would like to resolve the remaining discrepancy between the lower
DMC-ICE13 error obtained in your calculations and the periodic g-xTB results
from the CP2K-native implementation.  The attached package fixes the
structures, atom order, units, mesh convention, convergence settings, and DMC
references and provides absolute energies before any relative-energy
post-processing.

Three source states are kept separate in the package.  The first two revisions
were independently confirmed as the corresponding upstream branch tips on
2026-07-20:

- older-lineage `mstore-inorganic` at `be87ef681acd...`;
- the upstream author `pbc` snapshot at `c932120d258...`;
- the later pbc-derived CP2K-integration provider at `15915c943564...`, which
  is the provider linked into the qualified CP2K executable.

Our present evidence is:

1. With the exact current pbc-derived integration provider, direct `save_tblite` CLI and CP2K-native
   Bloch calculations are energetically identical within numerical noise.
   All 52 required points from `1 x 1 x 1` through `4 x 4 x 4` pass: the
   largest absolute total-energy difference is `1.12e-7 Eh` per primitive
   cell and the largest difference in an Ih-referenced relative energy is
   `2.14e-5 kJ mol-1 H2O-1`.
2. An explicit CP2K `2 x 2 x 2` Gamma supercell and the native reduced
   `2 x 2 x 2` Bloch mesh agree to `1.13e-11 Eh` per primitive cell for ice
   XVII.  This independently checks the k-to-BvK transformation and symmetry
   reduction.
3. The historical `mstore-inorganic` branch is not energetically identical to
   either `pbc` state.  With the same supplied cells, its DMC MAE changes from
   `48.7108` at `2 x 2 x 2` to `17.8306 kJ mol-1 H2O-1` at `3 x 3 x 3`, whereas
   the current `pbc` CLI gives `88.6814` and `34.0485 kJ mol-1 H2O-1`,
   respectively.  These sparse meshes are not converged accuracy estimates,
   but they show that the branch choice changes the model energies much more
   than the CP2K integration does.
4. Our independently rebuilt `mstore-inorganic` executable gives the same
   complete `3 x 3 x 3` absolute-energy matrix at `--acc 0.1` and `--acc 0.01`
   to within `6.5e-11 Eh` per explicit supercell, corresponding to only
   `2.26e-10 kJ mol-1 H2O-1` after same-mesh Ih referencing.  The large branch shift is
   therefore not caused by the SCC accuracy setting.

Could you please run the following two independent series with your own clean
builds?

- `lmseidler/save_tblite:mstore-inorganic`
- `lmseidler/save_tblite:pbc`

For the integration-parity check, please additionally build the exact
pbc-derived source revision `15915c9435644eb257178ca8f8bf7220c38b1a84`
recorded in `sources.json`.  This third build is important because the later
provider and the `c932120...` `pbc` snapshot are close but not bitwise or
energetically identical.

For each branch, please evaluate all thirteen supplied structures at
`2 x 2 x 2` and `3 x 3 x 3`, using `--method gxtb --acc 0.1`, and retain the
absolute `result.json` energy and the complete text output.  This is the exact
accuracy value used by the qualified CP2K-native inputs.  The package runner
is:

```text
python3 scripts/run_save_tblite.py /path/to/tblite <phase> <mesh> results/<branch> \
  --accuracy 0.1 --require-binary-sha256 <binary-sha256>
```

If the lower previously reported DMC-ICE13 error used `--acc 0.01`, please
repeat that setting as a separately labelled sensitivity matrix.  Do not mix
the two accuracies within one mesh: the complete text output is used to verify
the effective convergence thresholds independently.

Please also report:

- the exact source revision and SHA-256 of each executable;
- whether the calculation used the supplied explicit BvK POSCAR without
  additional cell replication or coordinate conversion;
- the precise settings used for the lower DMC-ICE13 error quoted previously;
- whether those values were produced by `mstore-inorganic`, `pbc`, or another
  source state.

The `c932120...` build should first be compared with
`tables/author_pbc_absolute_energies.csv`.  The `15915c...` build should then
be compared with the current-CLI column in
`tables/pbc_cli_vs_cp2k_native_absolute_parity.csv`; this is the strict
same-source CP2K-integration test.  Relative energies and DMC
statistics are secondary and are provided in the remaining tables and in
`comparison_workbook.xlsx`.

If your exact `15915c...` build reproduces the direct CLI energies in this package,
then the CP2K-native route is already consistent with its pbc-derived provider, and the lower
benchmark error must originate from a different model revision, input, or
post-processing convention.  If it does not reproduce them, the absolute
energy files will localize the first differing phase and mesh without relying
on cancellation against ice Ih.

Thank you for checking this independently.
