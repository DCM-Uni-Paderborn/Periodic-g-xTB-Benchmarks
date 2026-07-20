# Current build-consistent DMC-ICE13 status

This is the qualified Part-I state after completion of ice XI at `8 x 8 x 8`.
Every retained phase value and its same-mesh ice-Ih reference terminated
normally and carries CP2K executable SHA-256
`b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f`.
The one-step phase criterion is
`|Delta R| <= 0.10 kJ mol-1 per H2O`; the denser endpoint is retained.

| Phase | Retained mesh | Absolute DMC error | Latest `|Delta R|` | Converged | State |
|---|---:|---:|---:|:---:|---|
| II | `6^3` | 1.9831 | 0.0781 | yes | complete |
| III | `6^3` | 0.4222 | 0.0347 | yes | complete |
| IV | `6^3` | 2.5517 | 0.0085 | yes | complete |
| VI | `7^3` | 1.1303 | 0.0193 | yes | complete |
| VII | `8^3` | 5.4315 | 1.3847 | no | `8^3` cap reached; no `9^3` planned |
| VIII | `8^3` | 2.4198 | 0.0278 | yes | complete |
| IX | `6^3` | 0.7492 | 0.0180 | yes | complete |
| XI | `8^3` | 0.3488 | 0.0884 | yes | complete |
| XIII | `5^3` | 2.0638 | 0.0669 | yes | complete |
| XIV | `8^3` | 1.6095 | 0.2015 | no | `8^3` cap reached; no `9^3` planned |
| XV | `7^3` | 1.1041 | 0.0199 | yes | complete |
| XVII | `7^3` | 1.0777 | 0.0148 | yes | complete |

Ten of twelve phases pass: II, III, IV, VI, VIII, IX, XI, XIII, XV, and
XVII. Ice VII and XIV remain explicitly unresolved at the production cap.
No native phase or reference calculation is running or waiting.

## Aggregate comparison

- previous accepted same-build mixed MAE (before XI `8^3`):
  `1.733618175732 kJ mol-1 per H2O`;
- current same-build retained/capped mixed MAE:
  `1.740983997847 kJ mol-1 per H2O`;
- change caused by replacing XI `7^3` with the denser passing XI `8^3`
  endpoint: `+0.007365822115 kJ mol-1 per H2O`;
- earlier not fully build-consistent comparator:
  `1.749770892708 kJ mol-1 per H2O`;
- earlier paper data evaluated on the now-current retained mesh mix:
  `3.949012851607 kJ mol-1 per H2O`;
- reduction relative to that same-mesh paper comparator:
  `2.208028853760 kJ mol-1 per H2O`, or `55.913438%`.

Because VII and XIV do not pass, the current aggregate is a qualified capped
statistic rather than a final fully phase-wise-converged accuracy result.
