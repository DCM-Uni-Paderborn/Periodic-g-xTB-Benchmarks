# Post-#5582 DMC-ICE13 cell-dilation diagnostic

This directory records the controlled 0D-versus-3D g-xTB dilation test for ice Ih, VII, and XIV. The same twelve rigid water molecules are compared in a true nonperiodic 40 A cell and in a 3D-periodic cell with an unreduced, Gamma-centered 2x2x2 mesh. All reported cell-vector lengths are at most 40 A; VII at scale 4 is excluded because its longest vector would be 42.2194773541 A.

The full raw archive is stored persistently on Terok at:

`/home/kuehne88/work/gxtb-dmc-cell-dilation-post5582-20260716T0837Z`

The `analysis` directory contains the machine-readable full table, compact paper tables, LaTeX fragments, geometry/PBC audits, and the failure audit. `parameters` contains the exact exported g-xTB parameter file and each schema-supported diagnostic modification. `scripts` contains the deterministic parameter generator, campaign preparation/runner, and analysis code.

Build identity:

- CP2K revision: `28df9380abb327d56bbf216d2469a1fd8c953fc0`
- save_tblite revision: `257ba442684c39454175e5192c8a2342b4c6380f`
- CP2K executable SHA-256: `86949402a326e3118551fea079dd2b79df331087b47cb355b8803964f15deb89`
- save_tblite library SHA-256: `8ac8c98f462c6b29a2350ed341bf310addbc4692b0f6339a28ecb26c996c13a4`

Raw-campaign manifest SHA-256 values:

- full model: `aebbb9fcc9eac4f384989b55ee805e62b01d27784949702ca641f08613359034`
- no exchange: `44759b3e646796d8f2a64c22eb585316fe794f0a444bf0a70723ffd4a65475dd`
- frozen q-vSZP environment response: `c78130ea54f899a20e0933587af82adf60cc2ac43cf2a2a421e3300176f01163`
- no anisotropic multipole: `546bc609a85fcc9924b9b354b1d351540fcf4413fe0d6da03fbbf4e118d464e0`
- no ACP: `6950c1ce881cbdacfd8552e3e4e6d6fb5e58a5014f978ae5ff0d5d588274f1b4`
- no exchange plus no ACP interaction diagnostic: `6551bac51835933edc09e21bb04db0e5b70f1007073b6dbbb706a9674408dbef`
- successful component smokes: `b18288e3f5ac11420276d855dbd374b864612185ee1a8d2eedd74174eb25b684`
- archived long-path smoke failures: `ff9b719ff4867cc698a2e5174498a3b96bb213679c24b25b93b9a4da29227424`

The modified parameter files are component-deletion diagnostics, not reparameterized physical models. Removing the anisotropic multipole block does not disable all periodic electrostatics or periodic images. The combined exchange/ACP deletion exposes a CP2K integration-path abort before the first periodic SCC step and therefore yields no periodic model result; see `analysis/campaign_failure_audit.json`.
