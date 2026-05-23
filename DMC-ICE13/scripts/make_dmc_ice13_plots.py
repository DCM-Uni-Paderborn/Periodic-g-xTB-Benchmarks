#!/usr/bin/env python3
"""Create DMC-ICE13 tables and plots for the periodic GFN benchmark."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"

PHASES = ["Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"]

DMC_ABS = {
    "Ih": -59.45,
    "II": -59.14,
    "III": -58.20,
    "IV": -55.62,
    "VI": -57.67,
    "VII": -54.46,
    "VIII": -55.22,
    "IX": -58.85,
    "XI": -59.29,
    "XIII": -57.33,
    "XIV": -57.75,
    "XV": -57.71,
    "XVII": -57.70,
}

# Published single-point DFT values on the DMC-ICE13 structures, from
# Della Pia et al., J. Chem. Phys. 157, 134701 (2022), Table II.
PUBLISHED_ABS = {
    "B3LYP-D3atm": [-64.14, -64.02, -62.10, -60.89, -62.24, -59.89, -60.95, -63.39, -64.28, -62.99, -62.61, -62.16, -63.24],
    "optB86b-vdW": [-68.69, -67.89, -67.47, -65.74, -66.83, -62.84, -63.93, -68.22, -69.18, -67.68, -67.43, -66.63, -67.21],
    "SCAN+rVV10": [-68.26, -67.83, -66.33, -65.17, -66.27, -64.07, -65.53, -67.38, -68.52, -67.14, -66.90, -66.61, -67.29],
    "revPBE-D3": [-59.01, -57.75, -56.69, -55.04, -56.16, -54.83, -55.74, -57.41, -59.25, -56.71, -56.30, -56.07, -58.00],
    "RSCAN": [-61.40, -59.38, -58.67, -56.51, -57.35, -53.83, -55.04, -59.15, -61.68, -58.35, -57.90, -57.21, -60.75],
    "PBE-D4": [-69.62, -65.55, -66.60, -62.97, -62.71, -56.21, -57.35, -66.75, -70.10, -64.82, -63.96, -62.58, -68.61],
    "PBE": [-62.23, -55.54, -57.85, -52.85, -51.58, -43.79, -45.05, -57.56, -62.60, -54.54, -53.24, -51.56, -61.52],
}


def rel_from_abs(values: list[float]) -> dict[str, float]:
    ref = values[0]
    return {phase: value - ref for phase, value in zip(PHASES, values)}


def stats(errors: list[float]) -> dict[str, float]:
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(e) for e in errors) / len(errors),
        "RMSE": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "MaxAE": max(abs(e) for e in errors),
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_gnuplot(script: str) -> None:
    subprocess.run(["gnuplot"], input=script.encode(), check=True)


def main() -> None:
    DATA.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    results = json.loads((DATA / "results.json").read_text())

    dmc_rel = {phase: DMC_ABS[phase] - DMC_ABS["Ih"] for phase in PHASES}
    method_rel = {
        "GFN1-xTB": results["GFN1"]["relative_kjmol"],
        "GFN2-xTB": results["GFN2"]["relative_kjmol"],
    }
    method_rel.update({name: rel_from_abs(vals) for name, vals in PUBLISHED_ABS.items()})

    rows = []
    for index, phase in enumerate(PHASES, start=1):
        row = {
            "index": index,
            "phase": phase,
            "DMC_relative_kJmol": f"{dmc_rel[phase]:.4f}",
            "GFN1_relative_kJmol": f"{method_rel['GFN1-xTB'][phase]:.4f}",
            "GFN2_relative_kJmol": f"{method_rel['GFN2-xTB'][phase]:.4f}",
            "GFN1_error_kJmol": f"{method_rel['GFN1-xTB'][phase] - dmc_rel[phase]:.4f}",
            "GFN2_error_kJmol": f"{method_rel['GFN2-xTB'][phase] - dmc_rel[phase]:.4f}",
        }
        rows.append(row)
    write_csv(
        DATA / "dmc_ice13_relative_energies.csv",
        rows,
        ["index", "phase", "DMC_relative_kJmol", "GFN1_relative_kJmol", "GFN2_relative_kJmol", "GFN1_error_kJmol", "GFN2_error_kJmol"],
    )

    summary_rows = []
    for method, rel in method_rel.items():
        errors = [rel[phase] - dmc_rel[phase] for phase in PHASES if phase != "Ih"]
        row = {"method": method}
        row.update({key: f"{value:.4f}" for key, value in stats(errors).items()})
        summary_rows.append(row)
    summary_rows.sort(key=lambda row: float(row["MAE"]))
    write_csv(DATA / "dmc_ice13_relative_mae_comparison.csv", summary_rows, ["method", "ME", "MAE", "RMSE", "MaxAE"])

    rel_dat = DATA / "relative_energies_for_plot.dat"
    with rel_dat.open("w") as handle:
        handle.write("# index phase DMC GFN1 GFN2 revPBE-D3 optB86b-vdW SCAN+rVV10\n")
        for index, phase in enumerate(PHASES, start=1):
            handle.write(
                f"{index} {phase} {dmc_rel[phase]:.6f} "
                f"{method_rel['GFN1-xTB'][phase]:.6f} {method_rel['GFN2-xTB'][phase]:.6f} "
                f"{method_rel['revPBE-D3'][phase]:.6f} {method_rel['optB86b-vdW'][phase]:.6f} "
                f"{method_rel['SCAN+rVV10'][phase]:.6f}\n"
            )

    err_dat = DATA / "relative_errors_for_plot.dat"
    with err_dat.open("w") as handle:
        handle.write("# index phase GFN1_error GFN2_error\n")
        for index, phase in enumerate(PHASES, start=1):
            handle.write(
                f"{index} {phase} {method_rel['GFN1-xTB'][phase] - dmc_rel[phase]:.6f} "
                f"{method_rel['GFN2-xTB'][phase] - dmc_rel[phase]:.6f}\n"
            )

    mae_dat = DATA / "relative_mae_for_plot.dat"
    with mae_dat.open("w") as handle:
        handle.write("# index method MAE\n")
        for index, row in enumerate(summary_rows, start=1):
            handle.write(f"{index} \"{row['method']}\" {row['MAE']}\n")

    common = """
set terminal svg enhanced font 'Helvetica,12' size 980,500
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set key outside right center spacing 1.2 samplen 1.8
set style line 1 lc rgb '#111111' lw 2.2 pt 7 ps 0.75
set style line 2 lc rgb '#c44e52' lw 1.9 pt 5 ps 0.65
set style line 3 lc rgb '#4c72b0' lw 1.9 pt 9 ps 0.65
set style line 4 lc rgb '#55a868' lw 1.4 pt 13 ps 0.55
set style line 5 lc rgb '#8172b3' lw 1.4 pt 11 ps 0.55
set style line 6 lc rgb '#ccb974' lw 1.4 pt 15 ps 0.55
"""
    run_gnuplot(
        common
        + f"""
set output '{FIGURES / 'dmc_ice13_relative_energies.svg'}'
set ylabel 'Relative energy to ice Ih / kJ mol^{-1}'
set xlabel 'Ice polymorph'
set xrange [0.5:13.5]
set xtics rotate by -45
set yrange [-20:25]
plot '{rel_dat}' using 1:3:xtic(2) w lp ls 1 title 'DMC', \\
     '' using 1:4 w lp ls 2 title 'GFN1-xTB', \\
     '' using 1:5 w lp ls 3 title 'GFN2-xTB', \\
     '' using 1:6 w lp ls 4 title 'revPBE-D3', \\
     '' using 1:7 w lp ls 5 title 'optB86b-vdW', \\
     '' using 1:8 w lp ls 6 title 'SCAN+rVV10'
"""
    )

    run_gnuplot(
        common
        + f"""
set output '{FIGURES / 'dmc_ice13_relative_errors.svg'}'
set ylabel 'Relative-energy error vs DMC / kJ mol^{-1}'
set xlabel 'Ice polymorph'
set xrange [0.5:13.5]
set xtics rotate by -45
set yrange [-24:22]
set yzeroaxis lw 1.2 lc rgb '#333333'
plot '{err_dat}' using 1:3:xtic(2) w lp ls 2 title 'GFN1-xTB', \\
     '' using 1:4 w lp ls 3 title 'GFN2-xTB'
"""
    )

    run_gnuplot(
        """
set terminal svg enhanced font 'Helvetica,12' size 720,480
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set output '"""
        + str(FIGURES / "dmc_ice13_relative_mae_comparison.svg")
        + """'
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set style fill solid 0.85 border -1
set boxwidth 0.72
set ylabel 'MAE of relative energies / kJ mol^{-1}'
set xlabel 'Method'
set xrange [0.3:10.7]
set yrange [0:12]
set xtics rotate by -45
unset key
plot '"""
        + str(mae_dat)
        + """' using 1:3:xtic(2) with boxes lc rgb '#4c72b0'
"""
    )

    for svg in FIGURES.glob("*.svg"):
        png = svg.with_suffix(".png")
        subprocess.run(["rsvg-convert", str(svg), "-o", str(png)], check=True)


if __name__ == "__main__":
    main()
