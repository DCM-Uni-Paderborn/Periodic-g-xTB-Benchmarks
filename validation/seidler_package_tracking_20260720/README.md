# Seidler-package tracking gate

The recalculation package has its own `SHA256SUMS`, but a valid local hash
manifest alone does not prove that ignored `*.out` or `*.log` files were
published. This repository-level gate therefore checks five independent
conditions:

- every manifest entry exists;
- every package file other than the manifest itself is listed;
- every listed file is in the Git index;
- every digest matches; and
- all 176 raw CP2K/direct-CLI text outputs are present and tracked.

Run from the repository root with

```bash
python3 validation/seidler_package_tracking_20260720/verify_package_tracking.py
```

This check is intentionally Git-aware so that a fresh clone cannot silently
lose raw outputs through repository-wide ignore patterns while retaining a
manifest that only passes in the original working tree.
