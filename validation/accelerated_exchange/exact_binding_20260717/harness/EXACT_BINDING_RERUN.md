# Clean exact-binding qualification rerun

The original schema-1 `taskset`/shared-mask outputs, if restored under `runs/`,
remain immutable numerical evidence only. Their wall times are classified
`legacy_timing_non_scaling`. The new runner writes to
`runs_v2_exact_binding/` by default and refuses every nonempty per-run
directory, so it cannot overwrite a historical run.

The fixed matrix contains 48 runs: 16 case/rank combinations times the
`DENSE`, `STREAMED`, and `QUALIFY` selectors. Before launch, provide exactly
32 distinct, currently available, single-PU logical CPUs in the desired order
(eight worker slots of four CPUs). No range notation is accepted.
The supplied inputs are read directly from `test_inputs/`; no staging or copy
step is required.  Every CP2K invocation runs in its own new per-run directory,
which isolates all project-relative restart, auxiliary, and temporary files.

```bash
export CP2K_EXE=/path/to/qualified/cp2k.psmp
export CP2K_LIB=/path/to/qualified/libcp2k.so
export MPIEXEC_EXE=/path/to/the/same/qualified/mpiexec
export RUN_ROOT=runs_v2_exact_binding
export ORDERED_PE_RESERVATION=96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127

python3 run_test_matrix.py
RUN_ROOT=runs_v2_exact_binding python3 verify_test_matrix.py
```

The verifier writes `runs_v2_exact_binding_summary.tsv` beside, not inside,
the raw run directory. Likewise, verification of restored schema-1 data writes
`runs_summary.tsv` outside `runs/`, leaving every historical raw file untouched.

The runner holds a host-local `flock` for every reserved CPU and a writer lock
for the new evidence root.  After those locks are acquired and again before
each individual launch, Linux `/proc` is scanned for live non-zombie CP2K or
MPI-rank masks that overlap the requested CPUs.  The exact launcher contract is
`mpiexec --map-by pe-list=<literal-list>:ordered --bind-to core
--report-bindings -np <ranks> <cp2k> -i <absolute-input>`; no user-supplied MPI
arguments are accepted.  Each rank is monitored by `OMPI_COMM_WORLD_RANK`
throughout its lifetime and cross-checked against the complete Open MPI binding
report.  Sequential same-rank/same-mask PID generations are aggregated.  Rank
migration, changed successor masks, unranked CP2K processes, and concurrently
live duplicate-rank PIDs remain sticky failures.  Raw child histories and any
temporal duplicate-rank samples are persisted; the verifier reconstructs every
aggregate before accepting timing evidence. All inherited `OMPI_MCA_*` and `PRTE_MCA_*`
variables, including indirect MCA parameter-file selectors, are removed before
launch. On a proof failure the complete process group receives `SIGTERM`, then
`SIGKILL` after the timeout; CPU locks remain held until no live group member
remains and the launcher has been reaped.
