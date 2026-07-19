# Reference-build provenance

The hashes below identify the executables that produced the supplied absolute
energies.  Source revisions are recorded here, in the benchmark repository,
and are deliberately not embedded in the manuscript or supporting information.

| role | source branch or revision | executable SHA-256 |
|---|---|---|
| current direct `save_tblite` reference | `acp-projector-cache`, `15915c9435644eb257178ca8f8bf7220c38b1a84` | `c4c6b31546e3da4bb906f08aeef7ae123a1eba1c71b93d13db31a8bac528190c` |
| unchanged Seidler `pbc` comparison build | `c932120d2580811901de6a1fe3f89b943c251766` | `795ba8516910892dddb97fbccc319c7b14bf0ba46cc6e2a06c72d08deec41f5c` |
| authors' periodic-exchange diagnostic, macOS | `mstore-inorganic`, `be87ef681acd880705d83b8b1f7c19b58ca5ea85` | `324c2c1e4968eab579fae1bd8571a467d62a8eaf372f2b88906bb0d9f7ba7549` |
| authors' periodic-exchange diagnostic, Linux | `mstore-inorganic`, `be87ef681acd880705d83b8b1f7c19b58ca5ea85` | `4fa6fd99e1b0de2d0aa76b80cc9089a0ceeefdaf1bc787042221c7fb63479ffd` |
| post-March molecular g-xTB diagnostic | `1d06f6d0a973e78ae3522ddfd30b8b9056d5cdc1` | `c87471101170b506dae7f54700d5724aad9ce3dc5923e48d5317a4fd8f6cac60` |
| DCM `main` diagnostic | `4d614699849616d3fdd855ef9d886d34873b6759` | `2af03fdc70875df823038e49319f69751ae4a94dada58ce2960d09d358884bf0` |
| CP2K-native pre-response reference executable | `symmetry-fused-exchange`, `8520b2e592cd04d35081ab4ad46d92c606071e23` | `e034824111011b1177ed78f77f6049eeae5aca56dd7d96dfa923af0e29495b8d` |
| response-corrected CP2K-native executable | base `8520b2e592cd04d35081ab4ad46d92c606071e23`, effective source `10f4e246967de34cfcb2d26b9d587ecc00bbb2f4` | `b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f` |

The response-corrected CP2K executable reports the listed base revision.  Its
archived build tree additionally contained the two hardening edits recorded by
the source patch with SHA-256
`6fd50608c4837866032df666e33b39f43a6334e6cee536f3cabf5211119db279`;
those edits were subsequently committed unchanged in the listed effective
source revision.  The direct recalculation does not depend on CP2K; it only
requires one of the explicit Cartesian POSCAR files and a `save_tblite` CLI
executable.

The additional g-xTB builds are deliberately labelled as diagnostics.  Their
coarse-grid DMC-ICE13 values identify source-history effects but do not replace
the converged native reference or establish a preferred model revision.
