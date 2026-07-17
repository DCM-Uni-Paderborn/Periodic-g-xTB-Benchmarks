# KGROUP_PARTIAL_DISTRIBUTED_IMAGES independent oracle harness

This temporary, uncommitted harness compares the selectable CP2K
`KGROUP_PARTIAL_DISTRIBUTED_IMAGES` bridge with the unchanged `DENSE` bridge at identical MPI
rank counts.  It freezes the already-qualified KGROUP_OWNER inputs and covers:

- CH4 full `2x2x2`, P=1/2/4;
- CH4 K290 and SPGLIB reductions, P=1/2;
- H2 and UKS O2 time-reversal `3x1x1`, P=1/2/4, including P>nfull
  and therefore a deliberately empty importer;
- Ar2 1D and Ar4 2D periodicity, P=1/2;
- shifted-grid Si, P=1/2/4.

Every distributed-image output must contain the in-process complete-mesh forward
oracle (`dE`, `dVsh`, folded Fock response) and reverse oracle (overlap
adjoint, direct force, direct stress) with each residual at most `1e-10`.
Final energy, atomic forces, and analytical stress are parsed independently
from DENSE and distributed-image outputs.  Return code zero and exactly one
`PROGRAM ENDED` marker are mandatory.

The internal forward/reverse oracle ceiling is `1e-10`.  The independently
printed DENSE-vs-distributed-image observables use the SI acceptance limits:
`1e-9` Ha for the total energy, `1e-7` Ha/bohr for atomic forces, and
`1e-5` GPa = `0.1` bar for analytical stress.  Actual residuals are retained
at full parsed precision even when they are far below these limits.  All
cases use nontrivial full meshes, so the forward and reverse replicated
provider paths are exercised directly rather than by a one-point fallback.

Threading is forced to one for OpenMP, OpenBLAS, MKL, BLIS, GotoBLAS, and
Accelerate.  On terok the intended invocation pins every sequential MPI case
to one singleton mask per local MPI rank inside the reserved CPU set 96--127;
no cases are run concurrently.  The campaign uses `mpiexec --bind-to none`
and the fail-closed `gxtb_rank_taskset_96_127.sh` rank wrapper.

Example:

```sh
./run_matrix.py --cp2k /path/to/cp2k.psmp \
  --mpiexec-arg=--bind-to --mpiexec-arg=none \
  --rank-prefix '/home/kuehne88/work/gxtb_rank_taskset_96_127.sh'
./verify_matrix.py
./run_negative.py --cp2k /path/to/cp2k.psmp \
  --mpiexec-arg=--bind-to --mpiexec-arg=none \
  --rank-prefix '/home/kuehne88/work/gxtb_rank_taskset_96_127.sh'
./verify_negative.py
./freeze_manifest.py
sha256sum -c SHA256SUMS
```

The runners refuse to overwrite a nonempty run directory and the verifiers
fail on missing pairs, missing markers, metadata/hash drift, non-finite
values, a numerical gate violation, a fault-test timeout, a zero return code,
or a missing fail-closed diagnostic.
