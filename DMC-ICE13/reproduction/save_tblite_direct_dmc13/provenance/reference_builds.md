# Reference-build provenance

The hashes below identify the executables that produced the supplied absolute
energies.  Source revisions are recorded here, in the benchmark repository,
and are deliberately not embedded in the manuscript or supporting information.

| role | source branch or revision | executable SHA-256 |
|---|---|---|
| current direct `save_tblite` reference | `acp-projector-cache`, `15915c9435644eb257178ca8f8bf7220c38b1a84` | `c4c6b31546e3da4bb906f08aeef7ae123a1eba1c71b93d13db31a8bac528190c` |
| unchanged Seidler `pbc` comparison build | `c932120d2580811901de6a1fe3f89b943c251766` | `795ba8516910892dddb97fbccc319c7b14bf0ba46cc6e2a06c72d08deec41f5c` |
| CP2K-native reference executable | `symmetry-fused-exchange`, `8520b2e592cd04d35081ab4ad46d92c606071e23` | `e034824111011b1177ed78f77f6049eeae5aca56dd7d96dfa923af0e29495b8d` |

The CP2K executable was built from the committed Part-I native BvK Coulomb and
ACP-response implementation.  The direct recalculation does not depend on
CP2K; it only requires one of the explicit Cartesian POSCAR files and a
`save_tblite` CLI executable.
