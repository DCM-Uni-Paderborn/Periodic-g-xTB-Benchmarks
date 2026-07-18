# Reference-build provenance

The hashes below identify the executables that produced the supplied absolute
energies.  Source revisions are recorded here, in the benchmark repository,
and are deliberately not embedded in the manuscript or supporting information.

| role | source branch or revision | executable SHA-256 |
|---|---|---|
| current direct `save_tblite` reference | `acp-projector-cache`, `15915c9435644eb257178ca8f8bf7220c38b1a84` | `c4c6b31546e3da4bb906f08aeef7ae123a1eba1c71b93d13db31a8bac528190c` |
| unchanged Seidler `pbc` comparison build | `c932120d2580811901de6a1fe3f89b943c251766` | `795ba8516910892dddb97fbccc319c7b14bf0ba46cc6e2a06c72d08deec41f5c` |
| CP2K-native reference executable | `symmetry-fused-exchange`, base `fc1a7e79a1c256ec4e5b555009ad6751b08243b7` plus the Part-I native BvK Coulomb correction | `946cb0755078d2b609872da58faf56862dea3f08dbdb1e4ee8b901bfef648f73` |

The CP2K executable hash is the authoritative identifier until the native BvK
correction is committed on the public integration branch.  The direct
recalculation does not depend on CP2K; it only requires one of the explicit
Cartesian POSCAR files and a `save_tblite` CLI executable.
