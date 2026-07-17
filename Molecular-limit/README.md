# H2O molecular/periodic-limit check for periodic g-xTB

This directory archives a reproducible comparison of standalone `save_tblite`,
CP2K with no periodic boundary conditions, and CP2K with an implicit-Gamma
three-dimensional cubic cell.  It does not modify either source tree or the
manuscript.

## Frozen provenance

- Host: `terok.casus.science`; one MPI singleton process and
  `OMP_NUM_THREADS=1` per case.
- CP2K source: branch `g-xTB-pbc`, commit
  `18d37c946413dba1b848f57563c46d16b866ce20`, clean worktree.
- `save_tblite` source: branch `cp2k-integration`, commit
  `1449febde312874cd0fac4227919f5ba4e4b69b8`, clean worktree.
- CP2K executable SHA-256:
  `c6b51be7e356170dcb39a597d0e389bd701586e6131365ba317da3968c36eea7`.
- `tblite` executable SHA-256:
  `e590bce468964fb3d6ab8a4ffe9d5313f0b6a4ce87be9cf1456c6b6a8be1c690`
  (`tblite version 0.6.0`).
- Static `libtblite.a` SHA-256:
  `19dc000c604529048a99d7590193e68e026b58d2326e995369be7ebb0e65b577`.
- Main input SHA-256:
  `bfd0bca70f5ed011f4f9182cb893539967838541e604e212e46c3949cef1537c`.
- Tight input SHA-256:
  `f0383b66bf7e6c49cf513364c676564c9b602eff3ad2fc7f0480f4afaff9e4ca`.
- Shifted input SHA-256:
  `5e46e48b2c61c4d1f5caab4b055497e2f4a3428d2a3b7d2ab40e191938eba3cb`.
- The complete raw-file manifest is `SHA256SUMS.raw`, whose own SHA-256 is
  `fcbf1f73694a4d29526a5a9cb400a0862ce67e695ccef98ef1bd5b8038040ae7`.

The geometry (angstrom) is

```text
O   0.000000  0.000000  0.000000
H   0.758602  0.000000  0.504284
H  -0.758602  0.000000  0.504284
```

The production sequence uses g-xTB, `ACCURACY 0.1`, native
`SCC_MIXER TBLITE`, 300 mixer iterations, `EPS_DEFAULT=1e-12`, and
`EPS_SCF=1e-10`.  Every CP2K job runs `ENERGY_FORCE` and independently launches
the frozen standalone CLI through `REFERENCE_CLI`, retaining the generated GEN,
gradient, JSON, and log files.  The CLI command uses the same method, accuracy,
charge, spin, solver, iteration limit, and an explicit electronic temperature
of 300 K.  `run_on_terok.sh` records the full command sequence.

## Main results

`results.csv` contains all cells (8, 10, 12, 15, 20, 30, 40, 50, 60, 80,
100, 150, and 200 angstrom).  Selected values are:

| Boundary/cell | CP2K energy / Eh | standalone CLI / Eh | CP2K-CLI / Eh | 3D-0D / kJ mol-1 | max component force change from 0D / eV A-1 |
|---|---:|---:|---:|---:|---:|
| 0D | -76.432502147208922 | -76.432502146455420 | 7.535e-10 | 0 | 0 |
| 3D, 10 A | -76.492086050520996 | -76.492086050556040 | 3.504e-11 | -156.437517 | 2.929e-1 |
| 3D, 20 A | -76.435005436312565 | -76.435005435640278 | 6.723e-10 | -6.572385 | 1.211e-2 |
| 3D, 30 A | -76.432560702365919 | -76.432560701615145 | 7.508e-10 | -0.153737 | 5.552e-4 |
| 3D, 40 A | -76.432502176456893 | -76.432502175703362 | 7.535e-10 | -0.0000768 | 5.323e-4 |
| 3D, 50 A | -76.432499889813869 | -76.432499889060239 | 7.536e-10 | +0.005927 | 4.871e-4 |
| 3D, 100 A | -76.432498875647767 | -76.432498874894122 | 7.536e-10 | +0.008589 | 4.865e-4 |
| 3D, 200 A | -76.432498752282740 | -76.432498751529081 | 7.537e-10 | +0.008913 | 4.864e-4 |

Across the complete main sequence, the largest CP2K-versus-CLI absolute energy
difference is `7.537e-10 Eh`; the largest component-wise analytical-gradient
difference is `2.137e-8 Eh/a0` (`1.10e-6 eV/A`).  Thus the CP2K provider and
standalone CLI agree much more closely than the molecular-versus-periodic
finite-size difference.

The periodic energy is not monotonic: it crosses the molecular value near
40 A.  For 60--200 A, a linear fit in `L^-3` has a maximum residual of
`2.52e-10 Eh` and extrapolates to an offset of about `+0.00896 kJ mol-1`
relative to the molecular calculation.  The force difference likewise tends
toward about `4.86e-4 eV/A` in the sampled range.  These small residual offsets
are therefore properties of the molecular versus periodic paths in
`save_tblite`, not evidence of a CP2K/CLI binding mismatch.  They should be
reported as observed numerical/model-path offsets, not silently called exact
convergence.

## Robustness checks

- Tightening g-xTB from `ACCURACY 0.1` to `0.001` and CP2K from
  `EPS_SCF=1e-10` to `1e-12` changes the energies by at most `2.84e-14 Eh` and
  the gradients by at most `1.72e-8 Eh/a0`.  It reduces the tight CP2K/CLI
  maximum gradient difference to about `5.04e-9 Eh/a0`, but does not remove the
  molecular/periodic offset.
- Translating the molecule by `(10,10,10) A` changes energies by at most
  `1.42e-14 Eh` and gradients by at most `2.36e-15 Eh/a0`; the placement across
  a cell boundary is not responsible for the offset.
- Replaying the 300 K molecular CLI command by hand produces a gradient file
  and JSON file bitwise identical to the CP2K-triggered CLI run.  Their hashes
  are respectively
  `06ee4494a8bbc4676de74c018e536734ddec5a2dbc8641a8b42c46c43ecc0309`
  and
  `de9780e582ecc37b9e68bdcd7d185c7b6a28f6eb16462bf34eb81c15412a1ad3`.
- Omitting `--etemp 300` and thereby using the native g-xTB default changes this
  closed-shell H2O energy by only `2.15e-11 Eh` and the largest gradient
  component by `2.41e-12 Eh/a0`.
- A diagnostic 300 A CP2K cell aborts before the first SCF update because CP2K
  supplies a non-finite Fock matrix to the native potential mixer.  It is not
  included in the numerical sequence or extrapolation.

For the tight molecular standalone calculation, the analytical forces
(`-gradient`) in eV/A are approximately

```text
O   ( 0.000000000, 0.000000000, -3.421390784)
H1  ( 1.952615109, 0.000000000,  1.710695392)
H2  (-1.952615109, 0.000000000,  1.710695392)
```

At 200 A their largest component-wise change is `4.864e-4 eV/A`, while the
tight CP2K-versus-periodic-CLI maximum component difference is only
`2.59e-7 eV/A`.

## Meaning of `H2O_gxtb_reference_cli.inp`

The CP2K regression input named `H2O_gxtb_reference_cli.inp` is a periodic
Gamma-parity test, not a molecular/0D reference.  It explicitly sets a 10 A
cell with `PERIODIC XYZ` and has no `KPOINTS` section, hence implicit Gamma.
`tb_write_reference_gen` writes a periodic DFTB+ GEN structure (`3 S`) plus all
three lattice vectors whenever any periodic direction is active.  The
reference implementation also rejects a multi-k-point mesh.  Consequently,
that test compares CP2K implicit-Gamma periodic g-xTB to standalone periodic
Gamma g-xTB.  The separate `H2O_0D` case in this archive is the actual
molecular comparison.
