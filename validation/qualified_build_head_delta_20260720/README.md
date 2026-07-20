# Qualified-build versus branch-head audit

The DMC-ICE13 production binary was built from the qualified CP2K and
`save_tblite` revisions recorded in `metadata.json`.  Each integration branch
subsequently gained one commit.  This gate archives the exact two diffs and
checks whether either successor commit can change the production energy path.

The CP2K successor changes two source paths:

- explicitly listed `GENERAL` k-point meshes now infer their regular
  Cartesian product before constructing the g-xTB BvK Coulomb cell; and
- inconsistent redundant aliases in a written k-point restart are projected
  onto the canonical BvK subspace instead of aborting the checkpoint.

The archived DMC-ICE13 production inputs use `MACDONALD` and cold-start
`SCF_GUESS MOPAC`; none contains `EXT_RESTART`, `SCF_GUESS RESTART`, or a
restart-file input.  A `PRINT/RESTART ON` section is allowed because it only
writes a checkpoint and does not activate the restart-transfer energy path.
Consequently, neither new CP2K branch enters those energy evaluations.  The
`save_tblite` successor changes only the numerical tolerance of the inherited
nonperiodic CeCl3 unit test; no provider source file changes.

This proves that updating the two branches does not by itself require a full
DMC-ICE13 energy rerun.  The new restart behavior still receives separate
same-mesh, cross-mesh, and cold-start tests before it is used operationally.

Run from the repository root with

```bash
python3 validation/qualified_build_head_delta_20260720/verify_qualified_build_head_delta.py
```
