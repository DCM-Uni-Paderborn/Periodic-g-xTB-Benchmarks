# DMC-ICE13 controlled 0D/3D cell-dilation diagnostic

This archive compares the same 12 rigid water molecules as a true nonperiodic 0D cluster and as a 3D-periodic primitive cell. The 0D inputs use `PERIODIC NONE`, the analytic Poisson solver and no k-point section; CP2K therefore passes three false periodic flags to save_tblite. The periodic inputs use `PERIODIC XYZ`, periodic Poisson, and an unreduced Gamma-centered 2x2x2 mesh.

All reported cell-vector lengths are <=40 A. VII at scale 4 is deliberately excluded because its largest periodic vector is 42.2194773541 A. The geometry audit aligns each periodic water by an integer lattice translation plus one global cluster translation and checks every O-H distance.

The three modified parameter files are component-deletion diagnostics, not physical replacements or reparameterized models. In particular, removing the anisotropic multipole block does not disable all periodic electrostatics or image interactions. The residual component in the CSV is the total energy minus all CP2K components printed separately; in the full model it contains exchange and other tblite interactions and must not be labeled as pure exchange.

Build: CP2K `28df9380abb327d56bbf216d2469a1fd8c953fc0`, save_tblite `257ba442684c39454175e5192c8a2342b4c6380f`.
