# g-xTB accuracy sensitivity

The archived DMC-ICE13 inputs associated with the earlier calculation use
`ACCURACY 0.01`, whereas the current Part-I production protocol uses
`ACCURACY 0.1`.  This directory tests that protocol difference explicitly
with the current `save_tblite` executable and the current CP2K-native
implementation.

Coverage:

- all 13 phases at Gamma with the direct CLI at `0.01`;
- all 13 phases at Gamma with CP2K-native g-xTB at `0.01`;
- Ih, VII, XI, and XIV on commensurate `2 x 2 x 2` BvK supercells with the
  direct CLI at `0.01`;
- the corresponding `0.1` raw results are retained in the package-level
  `results/current_save_tblite_cli` and `results/current_cp2k_native` trees.

At Gamma the largest change in a relative phase energy on changing `0.1` to
`0.01` is `5.6e-8 kJ mol-1` per water and the DMC-ICE13 MAE changes by only
`2.8e-8 kJ mol-1`.  On the selected `2 x 2 x 2` BvK cells the largest relative
energy change is `1.9e-8 kJ mol-1` per water.  The complete Gamma
CP2K-native/direct-CLI comparison at `0.01` has a maximum absolute-energy
difference of `3.08e-8` Hartree per primitive cell (RMS `1.60e-8` Hartree).
The CP2K-native Gamma energies at `0.01` and `0.1` are identical at the printed
precision.

Consequently, the input-level `ACCURACY` difference does not explain either
the reported reference deviation or the different apparent k-point
convergence.  Run `verify_accuracy.py` to recompute every statement above
from the archived raw outputs.
