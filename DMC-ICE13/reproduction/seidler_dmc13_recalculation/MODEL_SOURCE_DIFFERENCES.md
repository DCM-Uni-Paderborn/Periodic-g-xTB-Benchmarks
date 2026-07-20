# Why `mstore-inorganic` and `pbc` are separate model states

The older-lineage `mstore-inorganic` revision and the later `pbc` revision are
not two front ends for an otherwise identical Hamiltonian.  The source history
between the exact revisions recorded in `sources.json` contains model-relevant
changes in all of the following areas:

- several successive g-xTB implementations and a later synchronization with
  the tblite main-line data structures;
- q-vSZP basis construction, including coordination-number regularization and
  removal of a basis filter;
- effective-charge Coulomb terms, self-interaction removal, and the
  minimum-image form of the second-order term;
- nonlocal exchange, including its periodic implementation and the periodic
  strain derivative;
- H0 periodic strain response;
- three-dimensional ACPs and their gradients;
- Wigner--Seitz self-image indexing and other periodic image handling.

Representative model-path commits in chronological order are listed below.
They are identifiers for source provenance, not manuscript content.

| Revision | Source-history subject |
|---|---|
| `fe7479d094c6` | g-xTB implementation 040426 |
| `4bb996902209` | g-xTB implementation 180426 |
| `f9d230950fc5` | g-xTB implementation 210426 |
| `b4b2b79506e3` | regularize square-root coordination numbers in basis construction |
| `c141b1f6c5fb` | g-xTB implementation consistent with the tblite main branch |
| `694c80e4d3d4` | periodic effective-charge derivative correction |
| `975a19218be8` | self-interaction-removal update |
| `a88c76fc2d0f` | periodic exchange implementation |
| `1e5fb8bde39c` | remove basis filter and alternative smoothing functions |
| `4a3bc7a8488b` | H0 strain-response correction |
| `593d55321cb3` | three-dimensional ACP implementation |
| `69efcea732c0` | periodic model update |
| `30b04691e0af` | Wigner--Seitz self-image indexing correction |
| `083f22030f0b` | minimum-image variant of the second-order term |
| `143046c5185a` | three-dimensional exchange strain-response correction |

The complete list can be regenerated in a repository containing both source
revisions with

```sh
git log --reverse --format='%H %s' \
  be87ef681acd880705d83b8b1f7c19b58ca5ea85..c932120d2580811901de6a1fe3f89b943c251766 \
  -- src/tblite/xtb/gxtb.f90 src/tblite/xtb/h0.f90 \
     src/tblite/exchange src/tblite/acp src/tblite/coulomb/charge \
     src/tblite/basis/q-vszp.f90 src/tblite/xtb/singlepoint.f90 \
     src/tblite/wignerseitz.f90
```

The literal output for the branch tips observed on 2026-07-20 is retained in
`model_path_commits_mstore_to_pbc.txt`; the independently queried remote tips
are in `author_branch_heads_20260720.txt`.

The numerical evidence in this package is decisive: the complete `2^3` and
`3^3` relative-energy matrices differ by many kJ mol-1 per water molecule,
whereas the current pbc-derived CLI and CP2K-native implementation agree near
the numerical convergence threshold.  Consequently, a lower
`mstore-inorganic` DMC error cannot be used as evidence for an error in the
CP2K interface to the later pbc-derived provider.

| Source state | `2^3` MAE | `3^3` MAE |
|---|---:|---:|
| older-lineage `mstore-inorganic` | 48.7108 | 17.8306 |
| current pbc-derived CLI | 88.6814 | 34.0485 |
| CP2K-native with the same current provider | 88.6814 | 34.0485 |

All entries are in kJ mol-1 per H2O.  Both branches have now been evaluated at
CLI accuracy `0.1` for the complete `2^3` and `3^3` matrices.  In addition, an
independently rebuilt `mstore-inorganic` executable was evaluated
at `3^3` with both `0.1` and `0.01`; the largest absolute total-energy change
was only `6.5e-11 Eh` per explicit supercell, or `2.26e-10 kJ mol-1 H2O-1`
after same-mesh Ih referencing.  The corresponding controlled
pbc-derived sensitivity test changes total energies by less than `2e-10 Eh`.
The branch separation is therefore neither a mixed-setting artifact nor an
SCC stopping-threshold effect.  The settings of every row are exposed in the
absolute-energy and difference tables, while the complete sensitivity matrix
is archived under `validation/mstore_accuracy_equivalence_20260720/` in the
repository root.

The largest individual relative-energy
shift from `mstore-inorganic` to the current pbc-derived source is 148.0272 at
`2^3` and 76.4515 kJ mol-1 per H2O at `3^3` (ice VII).  By contrast, the
author-`pbc` snapshot and the later pbc-derived integration provider differ by
at most 0.1355 and 0.1005 kJ mol-1 per H2O on those meshes.  The source history
changes several coupled self-consistent terms, so the large `mstore`/`pbc`
shift must not be assigned to exchange, ACP, H0, or the basis alone without a
dedicated fixed-source ablation.
