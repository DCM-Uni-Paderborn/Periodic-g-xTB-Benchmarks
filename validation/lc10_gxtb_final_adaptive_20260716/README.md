# LC10 first-passing adaptive g-xTB endpoint (2026-07-16)

This snapshot freezes the ten-solid LC10 selection under the confirmed
one-step convergence rule. Starting from `3^3`, the first adjacent interval is
accepted when both absolute changes satisfy

- `|delta a0| <= 0.025 A`, and
- `|delta Ecoh| <= 0.25 kJ mol-1 atom-1`.

Exactly one passing interval is sufficient and the denser endpoint is
retained. There is no extra RMS or two-step gate. Applying this rule to every
available `3^3` through `9^3` result gives three selected `7^3`, four selected
`8^3`, and three selected `9^3` meshes. The complete decisions and adjacent
deltas are recorded in `lc10_final_aggregate.json`; the selected values and
hash-bound output paths are in `lc10_final_selected_values.csv`.

For MgS, the `6^3 -> 7^3` energy change is
`-0.245234850626 kJ mol-1 atom-1`, just inside the accepted threshold. The
`7^3` value uses the user-approved cleaned physical SCC branch; the collapsed
point is explicitly excluded. LiF first passes only at `8^3 -> 9^3` because
its preceding lattice-constant step fails even though that step's cohesive
energy change passes.

`raw/` contains the copied Terok equilibrium inputs, outputs, job lineage, and
campaign stamps for the accepted pairs. `status/` contains the fit and branch
decisions used for the `7^3`, `8^3`, and `9^3` endpoints. Available `10^3`,
`11^3`, and `12^3` calculations are preserved in this archive as sensitivity
data only; they do not replace an earlier first-passing endpoint.

Build identity:

- campaign fingerprint: `8edfa0417321680d8f22b648863b3d3e008ed2368596fc31250328bcb807cf55`
- CP2K revision: `28df9380abb327d56bbf216d2469a1fd8c953fc0`
- save_tblite revision: `257ba442684c39454175e5192c8a2342b4c6380f`

Run `sha256sum -c SHA256SUMS` from this directory to validate the snapshot.
