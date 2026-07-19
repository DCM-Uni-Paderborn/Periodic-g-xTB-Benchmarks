# Provenance of the supplied Gamma-only exchange archive

The externally supplied file
`gamma_only_dmc_ice13_x23b_exchange_20260606.zip` has SHA-256
`716aeda1d664d6d71d56b3ce1ff9a412d9fac7aab9cea428caa54bed6a9bd600`.
It is not an independent direct-`save_tblite` reference calculation by the
model authors.

The classification follows directly from the archive contents:

- the embedded README states that the package was created on 6 June 2026 from
  the local `/Users/tkuehne/dropbox/cp2k-g-xTB` source tree;
- it describes all DMC-ICE13 entries as Gamma-only CP2K calculations;
- `dmc_ice13/gamma_total_and_relative_energies.csv` names local
  `dmc-ice13-work/runs/...after_gxtb_periodic_fix_fullgrid_nosym_20260603`
  output paths as its energy sources;
- the representative Ih g-XTB input has SHA-256
  `e11e132fa70978e99039081bd5954501e1b168cdcb20cb860bc6c8a8e0c0521f`
  and requests a full, unreduced `1 x 1 x 1` CP2K mesh with
  `ACCURACY 0.01`.

The archive remains useful for two narrowly defined purposes: it preserves a
historical pre-qualification CP2K Gamma snapshot, and it supplies a separately
packaged copy of the DMC-ICE13 geometries.  The periodic pair-distance audit
shows that these geometries and the current production inputs agree to better
than `1.2e-9` Angstrom.  No energy from this archive is used as evidence for
direct-CLI/native parity or for the final adaptive DMC-ICE13 MAE.
