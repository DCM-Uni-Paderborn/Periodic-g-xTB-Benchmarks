# Seidler-package tracking gate

The recalculation package has its own `SHA256SUMS`, but a valid local hash
manifest alone does not prove that ignored `*.out` or `*.log` files were
published. This repository-level gate therefore checks five independent
conditions:

- every manifest entry exists;
- every package file other than the manifest itself is listed;
- every listed file is in the Git index;
- every digest matches; and
- all raw CP2K/direct-CLI text outputs are present and tracked; and
- every nonfailed run directory identified by an `exit_status` or `result.json`
  marker contains a published `.out` or `.log` file.

The expected run set is derived from the package tree. This avoids a brittle
fixed file count while still failing if a newly added qualified endpoint lacks
its raw text output.

Run from the repository root with

```bash
python3 validation/seidler_package_tracking_20260720/verify_package_tracking.py
```

This check is intentionally Git-aware so that a fresh clone cannot silently
lose raw outputs through repository-wide ignore patterns while retaining a
manifest that only passes in the original working tree.
