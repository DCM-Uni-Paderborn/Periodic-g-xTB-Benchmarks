# Binary/provider identity gate

This gate records the build-level identity behind the direct `save_tblite`
CLI/CP2K-native energy comparison.  It supplements numerical parity by
checking that both executables were built from the same clean provider source
state and that the exact same static `libtblite.a` archive entered both link
paths.

The snapshot was captured on the Linux production host while the qualified
DMC-ICE13 calculations were running.  The relevant facts are:

- the direct CLI source tree is clean at the provider revision recorded by
  the Part-I provenance contract;
- the CP2K CMake cache selects the `SAVE` provider and records the same
  provider revision;
- the CLI link rule depends on the build-tree `libtblite.a`;
- the build-tree, installed, and CP2K-linked `libtblite.a` SHA-256 values are
  identical;
- both executables resolve the same compiler runtime and OpenBLAS libraries.

The qualified CP2K revision predates one later hardening commit.  Its source
diff was audited explicitly: the interface change applies only to explicitly
listed `GENERAL` meshes and the restart change only to restart-density alias
projection.  The production inputs use `MACDONALD` meshes, `SCF_GUESS MOPAC`,
and `RESTART OFF`, so neither later branch enters the DMC-ICE13 energy path.

Run the internal consistency check with

```bash
python3 validation/binary_provider_identity_20260720/verify_binary_provider_identity.py
```

The gate does not replace the energy comparison.  It excludes a different
provider binary as a hidden explanation for the observed numerical parity.
