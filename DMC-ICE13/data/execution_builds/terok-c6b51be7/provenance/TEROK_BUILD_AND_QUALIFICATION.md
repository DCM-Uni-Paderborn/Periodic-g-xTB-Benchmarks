# terok g-xTB/CP2K build and qualification provenance

Date: 2026-07-14 (Europe/Berlin)

Remote host: `kuehne88@terok.casus.science`

Workspace:

`/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb`

This workspace is isolated from existing user installations. No production DMC13
jobs were launched during build qualification.

## Host

- Debian GNU/Linux 13.5 (trixie), x86_64
- 240 online CPUs: four sockets, 60 cores/socket, one thread/core
- CPU model reported by the VM: AMD EPYC Processor, AVX2
- one NUMA node
- 503 GiB RAM
- `/home`: 1 TB, 711 GB free at qualification time
- `/tmp`: 252 GB tmpfs
- no Slurm, PBS, or LSF scheduler detected

Build parallelism was deliberately limited to 16 jobs. Qualification calculations
used one CPU core per process.

## Exact sources

| component | branch / role | revision | clean `git archive HEAD` SHA-256 |
|---|---|---|---|
| CP2K | `g-xTB-pbc` | `18d37c946413dba1b848f57563c46d16b866ce20` | `c17a5fa907d91ed7f1f2cc5c3d3f506cb4402aecd6d1ae7ffe613a3990403542` |
| save_tblite | `cp2k-integration` | `1449febde312874cd0fac4227919f5ba4e4b69b8` | `cd6daa00087da6af709e55f566a4276b98a988fe892c0aed076aa3adaad71677` |
| DFTD4 dependency | vendored clean source | `99d64ee83832fda03df2c7dd7b7fd9f4c9ebb098` | `0f4edd511a2b377a49da281a046024575c61da2abeec1b179e5bc04a276969da` |

The DFTD4 source was transferred from the already validated local dependency
checkout because the upstream FetchContent URL requested GitHub credentials on
the remote host. No source patch was applied; both the vendor and build checkouts
are clean.

## Toolchain and environment

The build uses a private micromamba prefix at `env/`; no system packages or
pre-existing environments were modified.

- GCC/GFortran 15.2.0 (conda-forge)
- Open MPI 5.0.10
- CMake 4.3.3
- Ninja 1.13.2
- Python 3.12.13
- OpenBLAS 0.3.33
- ScaLAPACK 2.2.0
- DBCSR 2.9.1 (`mpi_openmpi`)
- FFTW 3.3.11
- LibXC 7
- Spglib 2.7

Exact environment lock:

`manifests/conda-explicit-linux-64.txt`

SHA-256: `df77c7f116b85ff97e0f2760eb39d355ddc29ef70322fe91381ad80a3ecebc97`

## save_tblite configuration

- `Release`, static library, PIC
- tests enabled
- OpenMP enabled
- ddX, HDF5, and TREXIO disabled
- DFTD4 supplied through
  `FETCHCONTENT_SOURCE_DIR_DFTD4=build/dftd4-fetch-source`
- install prefix: `install/save_tblite`

Installed artifacts:

| artifact | SHA-256 |
|---|---|
| `install/save_tblite/bin/tblite` | `e590bce468964fb3d6ab8a4ffe9d5313f0b6a4ce87be9cf1456c6b6a8be1c690` |
| `install/save_tblite/lib/libtblite.a` | `19dc000c604529048a99d7590193e68e026b58d2326e995369be7ebb0e65b577` |

CTest: 98/100 passed. All g-XTB, periodic, ACP, and exchange tests passed.
The two failures are solely tests that require the intentionally disabled ddX
feature (`C-API` ddX cases and `xtbml-energy-sum-up-gfn1`). A g-XTB water CLI
smoke converged to `-76.437387426460 Eh`.

## CP2K configuration

- `Release`, shared CP2K library
- MPI enabled, MPI Fortran-2008 bindings disabled
- FFTW3, LibXC, Spglib, ScaLAPACK and DBCSR enabled
- save_tblite selected through `CP2K_TBLITE_PROVIDER=SAVE`
- pinned `CP2K_TBLITE_REVISION=1449febde312874cd0fac4227919f5ba4e4b69b8`
- g-XTB compile probe passed (`CP2K_TBLITE_HAS_GXTB=1`)
- CP2K's separate DFTD4, Libint2, LIBXS, ELPA, COSMA and GreenX paths disabled
- `DBCSR_USE_MPI=ON` explicitly supplied because the conda DBCSR CMake config
  links MPI but does not itself export that CP2K cache flag
- install prefix: `install/cp2k`
- install library directory: `install/cp2k/lib`

An initially relative CMake library-directory argument was canonicalized against
the login directory and created a new `~/lib`. The complete newly generated tree
was moved without overwriting anything to
`build/cp2k-misplaced-home-lib`, CMake was reconfigured with the absolute isolated
library directory, and CP2K was reinstalled. `~/lib` did not exist before and does
not exist after this correction.

Installed artifacts:

| artifact | SHA-256 |
|---|---|
| `install/cp2k/bin/cp2k.psmp` | `c6b51be7e356170dcb39a597d0e389bd701586e6131365ba317da3968c36eea7` |
| `install/cp2k/lib/libcp2k.so.2026.1` | `7813f8b2afddf5c355cbf41a06109d59f2872d5dd0be6645093c2f04489c1975` |

`cp2k.psmp --version` reports source `18d37c9` and flags
`omp fftw3 libxc parallel scalapack spglib s_dftd3 mctc-lib tblite tblite_gxtb`.
`ldd` contains no unresolved libraries and resolves `libcp2k.so.2026.1` from the
isolated install tree.

Qualification tests passed:

- `libcp2k_unittest.psmp`
- `grid_unittest.psmp <cp2k-source-root>`
- `memory_utilities_unittest.psmp`
- `orbital_transformation_matrices_unittest.psmp`

The generated `gx_ac_unittest.psmp` is not a relevant enabled test: it deliberately
aborts when CP2K is configured without the unrelated GreenX library.

## DMC13 cross-host qualification

All calculations use byte-identical validated inputs, one process, one OpenMP
thread, and one OpenBLAS thread in isolated `benchmark/DMC-ICE13` directories.

Initial Ih gate:

| mesh | input SHA-256 | macOS energy / Eh | terok energy / Eh | absolute delta / Eh | SCF steps |
|---|---|---:|---:|---:|---:|
| `1x1x1` | `d779ddf0934fe2b3198101cae51a81c099f98ef9ca483432585630560210ce25` | -918.879863282956762 | -918.879863282956762 | 0.0 at printed precision | 13 |
| `2x2x2` | `841882c10b4a811d639085e6c8a03da6e0a9ecf287c5f991f683c9a9c0cfc5c3` | -917.672960087289766 | -917.672960087289766 | 0.0 at printed precision | 13 |

The complete printed SCF energy trajectories also agree; only timing columns and
two sub-resolution last-step delta displays differ. CP2K wall times were 72.429 s
and 109.504 s on terok versus 42.531 s and 61.933 s on the Mac.

Final `6x6x6` qualification:

| phase | input SHA-256 | macOS energy / Eh | terok energy / Eh | terok - macOS / Eh | SCF steps (both) |
|---|---|---:|---:|---:|---:|
| Ih | `d508d1051a3e0d688b3b97d4f8196215693133f532023c77e1183265bd9e95fd` | -917.464649947946100 | -917.464649947946100 | 0.000000000000000 | 13 |
| VII | `25f89520f3384f2e28f0cc34b4963b32ae9650cef44d46850c95f22baa066893` | -917.472187146001261 | -917.472187146001147 | +0.000000000000114 | 12 |

Using 12 water molecules in each cell and
`1 Eh = 2625.4996394798254 kJ/mol`, the same-mesh VII-minus-Ih relative
energies are:

- macOS: `-1.6490758980427705549 kJ mol^-1 H2O^-1`
- terok: `-1.6490758980178283083 kJ mol^-1 H2O^-1`
- absolute cross-host delta:
  `0.0000000000249422466 kJ mol^-1 H2O^-1`

The requested tolerance of `0.001 kJ mol^-1 H2O^-1` is therefore passed by
about eight orders of magnitude. Full machine-readable values are in
`manifests/dmc13-k666-qualification.json`.

## Logs

Configuration, build, install, test, CLI smoke, version and `ldd` logs are retained
under `logs/`. Reference Mac outputs and qualification inputs/outputs are retained
under `benchmark/DMC-ICE13/`.

## Production-use note

There is no batch scheduler on terok. The current local DMC13 runner also contains
macOS-specific binary checks (`otool`/`.dylib`) and must not be pointed at this
Linux build unchanged. A Linux-aware provenance check and a coordinated process
pool are required before production. The 240 cores are shared with other users;
start only after checking load and agreeing on a reserve. No production jobs were
started by this qualification.
