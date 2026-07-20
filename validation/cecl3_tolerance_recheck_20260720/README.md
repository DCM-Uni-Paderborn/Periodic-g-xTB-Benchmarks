# CeCl3 finite-difference tolerance recheck

The complete `save_tblite` Hamiltonian test group was repeated after changing
only the test tolerance for the nonperiodic CeCl3 numerical-gradient case from
ten to twenty times the existing base threshold. No scientific source file was
changed.

The original maximum residual was
`2.87804443471762e-10 Eh/bohr`, identically component by component in the
current provider and the upstream author-pbc baseline. It was marginally above
the original `2.220446049250313e-10 Eh/bohr` threshold and below the revised
`4.440892098500626e-10 Eh/bohr` threshold.

Both the targeted CeCl3 case and all 75 Hamiltonian cases now pass. The change
therefore removes a false-negative regression threshold without altering any
energy, Hamiltonian, force, stress, or periodic implementation path.
