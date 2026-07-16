# Terok Linux MPI qualification

This directory is the Linux/x86-64 cross-platform qualification of the frozen
CP2K `KGROUP_OWNER` precursor.  The same CH4 `2x2x2` energy/force/stress input
was run serially with one, two, and four MPI ranks, one OpenMP thread per rank,
and Open MPI `--bind-to core --map-by core`.

All three executions used an in-process dense complete-mesh provider oracle in
SCF iteration 1.  They ended normally and reproduced the oracle with maximum
residuals of `2.553513e-15 Ha` in exchange energy, exactly zero in the shell
potential at printed precision, and `6.941217e-18 Ha` in the folded Fock
matrix.  Final total energies agree across rank counts to `7.11e-15 Ha`,
forces to `1.92e-16 Ha/bohr`, and printed stresses to `8.1e-7 bar`.

The timings are deliberately **not** evidence of scalable provider exchange.
The coupled save_tblite kernel remains a single-state calculation on one
global source, while the unchanged DBCSR overlap path is still entered by all
ranks.  Wall time is therefore flat within startup/noise (4.43--4.60 s), and
aggregate child CPU time increases from 3.47 to 11.49 s.  This is the expected
signature of a correctness/ownership precursor, not a speedup.

`summary.tsv` contains the numerical results.  `runs/*/cp2k.out` and
`runs/*/time.txt` are the raw outputs.  CP2K's own per-process peak-memory
estimate is 218--220 MiB; `child_maxrss_kb` is the Linux `RUSAGE_CHILDREN`
maximum recorded around the complete `mpirun` command, not a sum over ranks.

The Release binary linked successfully against the isolated save_tblite/ddX
install; `ldd` reported no missing libraries.  Exact source, provider, binary,
and input hashes plus the dependency resolution are in `provenance.txt`.
