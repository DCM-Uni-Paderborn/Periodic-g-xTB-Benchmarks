# Source and build identity

## CP2K source

- Local reviewed worktree:
  `/private/tmp/cp2k_gxtb_partial_root_bridge_clean_20260716`
- Terok copy:
  `/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_clean_20260716`
- Branch: `codex/gxtb-partial-root-bridge-clean`
- Base commit: `0a1f7e3329a3e6c2a6accff28617af53fb9943b4`
- Sole changed file: `src/tblite_interface.F`
- Changed-file SHA-256:
  `47c9b039b2e0d081f1ac3688f29f5c75ffed9a60acbb490f7bee1ae99593dd5d`
- Patience diff SHA-256:
  `8763981ef9c6ba7e9db26a11f4245e11fbcc52ac393a6b6421b09e8e7628156f`
- `git diff --check`: clean
- CP2K precommit status: OK (retained in `cp2k/source/precommit.txt`)

## `save_tblite` source and installs

- Qualified commit:
  `35e7942b60edd89bb407ab3da5768d3410af83f5`
- Provider source tar on Terok:
  `/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942.tar.gz`
- Provider source-tar SHA-256:
  `2535519767302bc851c30512852e4ca031fadc0f97d0ece860e983effedbfd28`
- Release install:
  `/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942_install_release`
- Release `libtblite.a` SHA-256:
  `20c74bf3272a229e125893956880757bca365c0dc1b54fa8e892bf99a67e7760`
- Debug install:
  `/home/kuehne88/work/save_tblite_partial_k_to_r_35e7942_install_debug`
- Debug `libtblite.a` SHA-256:
  `cf6aaf4b9174e9759ec7d6ec9ac91b0d660643318270c2229a9cd27fc6bc05df`

Both provider builds are static (`BUILD_SHARED_LIBS=OFF`).  Their exact CMake
caches and Ninja graphs are under `provider/build/`; both focused exchange
logs end with 31/31 passing tests.  The full-CTest base classification records
`WITH_DDX=OFF` explicitly.

## CP2K final builds

Release build:

- Path:
  `/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_35e7942_release_20260716`
- `bin/cp2k.psmp` SHA-256:
  `9b832420483ca29ac3adc05bd60316ac4431fb782a7259ac39ac4cdfcb39e75a`
- `src/libcp2k.so.2026.2` SHA-256:
  `b0b619a4b758b1baf35d417afbd9792703335c4a0177476dbab5adbc068b2818`
- Exported provider symbol:
  `__tblite_cp2k_compat_MOD_cp2k_exchange_partial_begin`

Debug build:

- Path:
  `/home/kuehne88/work/cp2k_gxtb_partial_root_bridge_35e7942_debug_20260716`
- `bin/cp2k.pdbg` SHA-256:
  `5c1074e41989b4b76ff348939a03b9b39ebe1796f01ed5ac988e00afd1dea42f`
- `src/libcp2k.so.2026.2` SHA-256:
  `f4fe80c2dcacf20d58546497bb04201765b03f615dcd7cda21662811c15bdf9a`
- Exported provider symbol:
  `__tblite_cp2k_compat_MOD_cp2k_exchange_partial_begin`

For each configuration, `CMakeCache.txt` records
`CP2K_TBLITE_PROVIDER=SAVE` and the configuration-specific provider install;
the `cp2k` link edge in `build.ninja` names that exact `libtblite.a`.  The
Release and Debug build logs both finish at `[4159/4159]` with the final CP2K
executable linked.

## Frozen runtime roots

- Final positive/fault/Gamma smokes:
  `/home/kuehne88/work/gxtb-partial-root-final-35e7942-20260716`
- Full provider CTest classification:
  `/home/kuehne88/work/save_tblite_ctest_classification_20260717_review`
- Final independent oracle matrix:
  `/home/kuehne88/work/gxtb_partial_root_oracle_35e7942_20260717`
- Oracle internal manifest SHA-256:
  `5eceecdaf3ec4bee2a2e25c47f062c8f62529f9d0a0cf4548edfb9b9143cf188`
- Provider-classification internal manifest SHA-256:
  `686fe504de6ce66bb4e5e8af8a189e8fde6b949fae5c4df394af69f7b3d9ec1a`
