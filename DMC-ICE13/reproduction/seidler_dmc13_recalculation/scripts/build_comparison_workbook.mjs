#!/usr/bin/env node
/** Build the review workbook from the machine-readable package CSV files. */

import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { FileBlob, SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const packageDir = process.env.GXTB_PACKAGE_DIR || path.dirname(scriptDir);
const tablesDir = path.join(packageDir, "tables");
const outputPath = process.argv[2] || path.join(packageDir, "comparison_workbook.xlsx");
const previewDir = process.argv[3] || "";

const imports = [
  ["Branch statistics", "branch_comparison_statistics.csv"],
  ["Relative energies", "all_branch_relative_energy_comparison.csv"],
  ["Branch differences", "mstore_vs_pbc_relative_differences.csv"],
  ["Absolute parity", "pbc_cli_vs_cp2k_native_absolute_parity.csv"],
  ["Native meshes", "cp2k_native_relative_energies_by_mesh.csv"],
  ["Mstore absolute", "mstore_inorganic_absolute_energies.csv"],
  ["References", "dmc_reference_relative_energies.csv"],
];

const workbook = Workbook.create();
const summary = workbook.worksheets.add("Summary");
for (const [sheetName, filename] of imports) {
  const csvText = await fs.readFile(path.join(tablesDir, filename), "utf8");
  const imported = await Workbook.fromCSV(csvText, { sheetName });
  const importedSheet = imported.worksheets.getItem(sheetName);
  const importedRange = importedSheet.getUsedRange();
  const sheet = workbook.worksheets.add(sheetName);
  if (importedRange.rowCount > 0 && importedRange.columnCount > 0) {
    sheet.getRangeByIndexes(0, 0, importedRange.rowCount, importedRange.columnCount).values = importedRange.values;
  }
}

const navy = "#17365D";
const blue = "#DCE6F1";
const teal = "#0F766E";
const orange = "#F4B183";
const paleGreen = "#E2F0D9";
const paleRed = "#FCE4D6";
const border = "#B4C6E7";

summary.showGridLines = false;
summary.getRange("A1:N1").merge();
summary.getRange("A1").values = [["DMC-ICE13 branch and CP2K-native comparison"]];
summary.getRange("A1:N1").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF", size: 16 },
  horizontalAlignment: "left",
  verticalAlignment: "center",
};
summary.getRange("A1:N1").format.rowHeight = 30;
summary.getRange("A3:N3").merge();
summary.getRange("A3").values = [[
  "Independent diagnostic package: historical mstore-inorganic, author pbc, current pbc CLI, and CP2K-native Bloch k points",
]];
summary.getRange("A3:N3").format = {
  fill: blue,
  font: { italic: true, color: navy },
  wrapText: true,
};

summary.getRange("A5:B5").values = [["Validated numerical check", "Result"]];
summary.getRange("A6:B11").values = [
  ["Maximum |CP2K-native - current pbc CLI| across 1^3-4^3 (Ha/primitive)", null],
  ["Maximum |CP2K-native - current pbc CLI| at 3^3 (Ha/primitive)", null],
  ["Maximum relative-energy difference at 3^3 (kJ mol-1 H2O-1)", null],
  ["CP2K native 2^3 vs explicit Gamma-BvK oracle (Ha/primitive)", 1.1255e-11],
  ["Complete absolute-parity points across 1^3-4^3", null],
  ["Interpretation", "Current pbc CLI and CP2K-native agree within numerical noise."],
];
summary.getRange("B6").formulas = [[
  "=MAX(MAX('Absolute parity'!F2:F53),-MIN('Absolute parity'!F2:F53))",
]];
summary.getRange("B7").formulas = [[
  "=MAX(MAX('Absolute parity'!F28:F40),-MIN('Absolute parity'!F28:F40))",
]];
summary.getRange("B8").formulas = [["=MAX(K21:K32)"]];
summary.getRange("B10").formulas = [["=COUNT('Absolute parity'!F2:F53)"]];
summary.getRange("A5:B5").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
};
summary.getRange("A5:B11").format.borders = { preset: "all", style: "thin", color: border };
summary.getRange("B6:B9").format.numberFormat = "0.000E+00";
summary.getRange("B10").format.numberFormat = "0";
summary.getRange("B11").format.wrapText = true;

const statsCsv = await fs.readFile(path.join(tablesDir, "branch_comparison_statistics.csv"), "utf8");
const statLines = statsCsv.trim().split(/\r?\n/).slice(1).map((line, index) => {
  const fields = line.split(",");
  return { method: fields[0], mesh: Number(fields[1]), mae: Number(fields[4]), excelRow: index + 2 };
});
const methods = [
  "CP2K-native pbc provider",
  "current pbc CLI",
  "author pbc CLI",
  "historical mstore-inorganic CLI",
];
const methodLabels = ["CP2K-native", "current pbc CLI", "author pbc", "mstore-inorganic"];
summary.getRange("A12:E12").values = [["Mesh", ...methodLabels]];
const helperValues = [];
const helperFormulas = [];
for (const mesh of [1, 2, 3]) {
  helperValues.push([mesh, null, null, null, null]);
  const formulas = [null];
  for (const method of methods) {
    const match = statLines.find((row) => row.method === method && row.mesh === mesh);
    formulas.push(match ? `='Branch statistics'!E${match.excelRow}` : "");
  }
  helperFormulas.push(formulas);
}
summary.getRange("A13:E15").values = helperValues;
summary.getRange("A13:E15").formulas = helperFormulas;
summary.getRange("A12:E12").format = {
  fill: navy,
  font: { bold: true, color: "#FFFFFF" },
};
summary.getRange("A12:E15").format.borders = { preset: "all", style: "thin", color: border };
summary.getRange("B13:E15").format.numberFormat = "0.0000";

const chart = summary.charts.add("line", summary.getRange("A12:E15"));
chart.title = "Low-mesh DMC-ICE13 branch comparison";
chart.hasLegend = true;
chart.xAxis = { axisType: "textAxis" };
chart.xAxis.title.text = "mesh n in n x n x n";
chart.yAxis = { numberFormatCode: "0.0" };
chart.yAxis.title.text = "MAE (kJ mol-1 H2O-1)";
chart.setPosition("G5", "N18");

summary.getRange("A20:E20").merge();
summary.getRange("A20").values = [["How to read the result"]];
summary.getRange("A20:E20").format = {
  fill: orange,
  font: { bold: true, color: "#000000" },
};
summary.getRange("A21:E24").merge(true);
summary.getRange("A21:A24").values = [
  ["1. Absolute parity isolates CP2K integration from model-revision effects."],
  ["2. Branch differences isolate mstore-inorganic from pbc at fixed inputs."],
  ["3. Relative energies always use ice Ih from the same mesh and method."],
  ["4. Blank entries are incomplete runs and are excluded from statistics."],
];
summary.getRange("A21:E24").format = { wrapText: true };
summary.getRange("A20:E24").format.borders = { preset: "outside", style: "thin", color: border };

const parityCsv = await fs.readFile(path.join(tablesDir, "pbc_cli_vs_cp2k_native_absolute_parity.csv"), "utf8");
const parityLines = parityCsv.trim().split(/\r?\n/).slice(1).map((line, index) => {
  const fields = line.split(",");
  return { mesh: Number(fields[0]), phase: fields[1], excelRow: index + 2 };
});
const parityPhases = ["II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"];
summary.getRange("J20:K20").values = [["3^3 phase", "|native - CLI| (kJ mol-1 H2O-1)"]];
summary.getRange("J20:K20").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF" },
  wrapText: true,
};
summary.getRange("J21:J32").values = parityPhases.map((phase) => [phase]);
summary.getRange("M20:M21").values = [["Hartree to kJ mol-1"], [2625.4996394798254]];
summary.getRange("M20").format = { fill: blue, font: { bold: true, color: navy }, wrapText: true };
summary.getRange("M21").format.numberFormat = "0.0000000000000";
const ihParity = parityLines.find((row) => row.mesh === 3 && row.phase === "Ih");
if (!ihParity) throw new Error("Missing 3^3 ice-Ih absolute parity source");
summary.getRange("K21:K32").formulas = parityPhases.map((phase) => {
  const point = parityLines.find((row) => row.mesh === 3 && row.phase === phase);
  if (!point) throw new Error(`Missing 3^3 absolute parity source for ${phase}`);
  return [
    `=ABS((('Absolute parity'!D${point.excelRow}/'Absolute parity'!C${point.excelRow})-` +
      `('Absolute parity'!D${ihParity.excelRow}/'Absolute parity'!C${ihParity.excelRow}))-` +
      `(('Absolute parity'!E${point.excelRow}/'Absolute parity'!C${point.excelRow})-` +
      `('Absolute parity'!E${ihParity.excelRow}/'Absolute parity'!C${ihParity.excelRow})))*$M$21`,
  ];
});
summary.getRange("J20:K32").format.borders = { preset: "all", style: "thin", color: border };
summary.getRange("K21:K32").format.numberFormat = "0.000000E+00";
summary.getRange("J:J").format.columnWidth = 13;
summary.getRange("K:K").format.columnWidth = 24;
summary.getRange("M:M").format.columnWidth = 22;
summary.getRange("A:A").format.columnWidth = 58;
summary.getRange("B:E").format.columnWidth = 19;
summary.freezePanes.freezeRows(3);

for (const [sheetName] of imports) {
  const sheet = workbook.worksheets.getItem(sheetName);
  const used = sheet.getUsedRange();
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  used.format.borders = { preset: "all", style: "thin", color: "#D9E2F3" };
  used.format.verticalAlignment = "center";
  const header = used.getRow(0);
  header.format = {
    fill: navy,
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
    verticalAlignment: "center",
  };
  header.format.rowHeight = 34;
  used.format.autofitColumns();
  const columnCount = used.columnCount;
  for (let col = 0; col < columnCount; col += 1) {
    const column = used.getColumn(col);
    const headerText = String(used.values[0][col] ?? "").toLowerCase();
    const longProvenanceText = headerText.includes("sha256") || headerText === "raw_result";
    if (column.format.columnWidth > 32) column.format.columnWidth = 32;
    if (column.format.columnWidth < 11) column.format.columnWidth = 11;
    if (longProvenanceText) {
      column.format.columnWidth = 24;
      column.format.wrapText = true;
    }
  }
  used.format.autofitRows();
  header.format.rowHeight = 34;
}

const statsSheet = workbook.worksheets.getItem("Branch statistics");
const statsUsed = statsSheet.getUsedRange();
statsSheet.getRange(`D2:G${statsUsed.rowCount}`).format.numberFormat = "0.0000";
statsSheet.getRange(`E2:E${statsUsed.rowCount}`).conditionalFormats.add("colorScale", {
  colors: [paleGreen, "#FFF2CC", paleRed],
});

const diffSheet = workbook.worksheets.getItem("Branch differences");
const diffUsed = diffSheet.getUsedRange();
if (diffUsed.rowCount > 1) {
  diffSheet.getRange(`D2:H${diffUsed.rowCount}`).format.numberFormat = "0.000000";
  diffSheet.getRange(`F2:F${diffUsed.rowCount}`).conditionalFormats.add("colorScale", {
    colors: ["#5B9BD5", "#FFFFFF", "#ED7D31"],
  });
}

for (const sheetName of ["Relative energies", "Native meshes"]) {
  const sheet = workbook.worksheets.getItem(sheetName);
  const used = sheet.getUsedRange();
  if (used.rowCount > 1) sheet.getRange(`D2:G${used.rowCount}`).format.numberFormat = "0.000000";
}
const paritySheet = workbook.worksheets.getItem("Absolute parity");
const parityUsed = paritySheet.getUsedRange();
if (parityUsed.rowCount > 1) paritySheet.getRange(`D2:F${parityUsed.rowCount}`).format.numberFormat = "0.000000000000E+00";

const inspection = await workbook.inspect({
  kind: "workbook,sheet",
  maxChars: 8000,
  tableMaxRows: 5,
  tableMaxCols: 8,
});
console.log(inspection.ndjson);
const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?",
  options: { useRegex: true, maxResults: 100 },
  maxChars: 4000,
});
console.log(formulaErrors.ndjson);

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
const exported = await FileBlob.load(outputPath);
const reopened = await SpreadsheetFile.importXlsx(exported);
const reopenedCheck = await reopened.inspect({
  kind: "sheet,region",
  sheetId: "Summary",
  range: "A1:E24",
  maxChars: 5000,
});
console.log(reopenedCheck.ndjson);

if (previewDir) {
  await fs.mkdir(previewDir, { recursive: true });
  for (const sheetName of ["Summary", ...imports.map(([name]) => name)]) {
    const preview = await workbook.render({ sheetName, autoCrop: "all", scale: 1, format: "png" });
    const safeName = sheetName.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
    await fs.writeFile(path.join(previewDir, `${safeName}.png`), new Uint8Array(await preview.arrayBuffer()));
  }
}

// artifact-tool may materialize a diagnostic inspect sidecar beside an
// imported workbook.  It is a transient QA artifact, not part of the
// scientific reproduction package.
await fs.rm(`${outputPath}.inspect.ndjson`, { force: true });

console.log(`workbook=${outputPath}`);
