# LC10 periodic g-xTB benchmark

This directory contains the fixed ten-solid Part-I subset:

```text
C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF, LiCl
```

MgO and LiH are not members of the reported benchmark and no longer appear in
the current data tree.

## Selection rule

Starting at 3x3x3, each solid is advanced independently. The denser endpoint
of the first adjacent pair satisfying both

```text
|Delta a0|   <= 0.025 A
|Delta Ecoh| <= 0.25 kJ mol-1 atom-1
```

is retained. One passing interval is sufficient; no aggregate or two-step gate
is applied. All ten solids pass by at most 9x9x9.

## Paper data

- `data/lc10_gxtb_final.csv`: final per-solid lattice constants, cohesive
  energies, reference values, signed errors, retained meshes, and the actual
  first-passing changes;
- `data/lc10_method_comparison.csv`: compact g-xTB/GFN1-xTB/GFN2-xTB
  aggregate comparison over exactly the common ten-solid set;
- `data/lc10_adaptive_progress_mae.csv`: uniform 3x3x3--6x6x6 and subsequent
  adaptive-endpoint plot values;
- `data/reference_goldzak2022.csv`: zero-point-corrected experimental and
  high-level reference values;
- `figures/lc10_gxtb_adaptive_mae.*`: paper convergence figure.

The final g-xTB MAEs are 0.1434 A for lattice constants and
0.2947 eV atom-1 for cohesive energies.

The plotting script regenerates the repository figure directly from the
versioned seven-stage table. Generated working directories and raw standard
output are intentionally not versioned. Complete GFN1-xTB/GFN2-xTB raw data
remain in the separate `Periodic-GFN2-Benchmarks` repository.
