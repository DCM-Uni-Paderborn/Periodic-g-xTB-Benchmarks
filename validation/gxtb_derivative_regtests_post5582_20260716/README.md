# Post-#5582 periodic g-xTB derivative qualification

This archive freezes the analytical-force and analytical-stress qualification used by the periodic g-xTB paper. The production data were generated on `terok.casus.science` on 2026-07-16 with the post-CP2K-#5582 build described in
`../../campaigns/gxtb-pbc-v1-post5582-20260714/build_manifest.json`.

## Qualified build

- CP2K revision: `28df9380abb327d56bbf216d2469a1fd8c953fc0`
- CP2K executable SHA-256: `86949402a326e3118551fea079dd2b79df331087b47cb355b8803964f15deb89`
- save_tblite revision: `257ba442684c39454175e5192c8a2342b4c6380f`
- save_tblite static-library SHA-256: `8ac8c98f462c6b29a2350ed341bf310addbc4692b0f6339a28ecb26c996c13a4`
- CP2K finite-difference displacement and strain step: `1.0e-4`

`derivative_summary.csv` is the machine-readable paper table. Every production case ended with return code zero and `PROGRAM ENDED AT`. CP2K prints the summed absolute force-component residual with a resolution of `1e-8 Eh/a0`; a printed zero therefore means below that resolution, not mathematical identity. The stress residual is the summed absolute virial-component residual divided by the cell volume.

The production matrix covers implicit Gamma, explicit Gamma, a full `1x1x1` mesh, a full `2x2x2` mesh, K290 and SPGLIB reduction, a native `3x1x1` Bloch mesh, time-reversal reduction, and the commensurate Gamma supercell. The complete inputs and outputs are under `inputs/` and `runs/`.

## CPU-affinity incident

The first launch at 08:04 used eleven concurrent `mpirun -np 4` commands. Open MPI bound each independent job to the same four logical CPUs (0--3), so 44 ranks contended for four CPUs and each rank obtained only about 9--10% CPU. That campaign was stopped and is retained under `discarded_affinity_collision_0804/` only as diagnostic provenance; it is not a numerical or timing source for the paper.

The production relaunch at 08:30 added `--bind-to none`. Its ranks were distributed across the host and sustained roughly one full CPU per rank. All eleven cases then completed successfully. `SHA256SUMS` covers the entire archive, while every production run also contains its own initial and final SHA-256 manifests.
