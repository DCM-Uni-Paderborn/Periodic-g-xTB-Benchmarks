#!/usr/bin/env python3
"""Compare the archived author DMC-ICE13 cells with the production POSCARs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import Counter
from itertools import product
from pathlib import Path

import numpy as np


PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII")
ARCHIVE_ROOT = "gamma_only_dmc_ice13_x23b_exchange_20260606/dmc_ice13"
MATCH_TOLERANCE_ANGSTROM = 1.0e-8
LATTICE_METRIC_TOLERANCE_ANGSTROM2 = 1.0e-8


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_xyz(data: bytes) -> tuple[list[str], np.ndarray]:
    lines = data.decode("utf-8").splitlines()
    atom_count = int(lines[0].strip())
    records = [line.split() for line in lines[2:] if line.strip()]
    if len(records) != atom_count or any(len(record) < 4 for record in records):
        raise ValueError("malformed XYZ coordinate file")
    species = [record[0] for record in records]
    coordinates = np.array([[float(value) for value in record[1:4]] for record in records])
    return species, coordinates


def parse_poscar(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    scale = float(lines[1])
    if scale <= 0.0:
        raise ValueError(f"unsupported non-positive POSCAR scale in {path}")
    cell = scale * np.array([[float(value) for value in lines[index].split()[:3]] for index in range(2, 5)])
    elements = lines[5].split()
    counts = [int(value) for value in lines[6].split()]
    if len(elements) != len(counts):
        raise ValueError(f"element/count mismatch in {path}")
    coordinate_line = 7
    if lines[coordinate_line].lower().startswith("s"):
        coordinate_line += 1
    mode = lines[coordinate_line].lower()
    atom_count = sum(counts)
    values = np.array(
        [[float(value) for value in lines[index].split()[:3]] for index in range(coordinate_line + 1, coordinate_line + 1 + atom_count)]
    )
    if values.shape != (atom_count, 3):
        raise ValueError(f"coordinate-count mismatch in {path}")
    if mode.startswith("d"):
        coordinates = values @ cell
    elif mode.startswith(("c", "k")):
        coordinates = scale * values
    else:
        raise ValueError(f"unsupported coordinate mode in {path}: {lines[coordinate_line]}")
    species = [element for element, count in zip(elements, counts) for _ in range(count)]
    return cell, species, coordinates


def perfect_species_match(
    cell: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    tolerance: float,
) -> tuple[float, list[float]]:
    inverse = np.linalg.inv(cell)
    left_fractional = left @ inverse
    right_fractional = right @ inverse
    delta = left_fractional[:, None, :] - right_fractional[None, :, :]
    delta -= np.rint(delta)
    distances = np.linalg.norm(delta @ cell, axis=2)
    adjacency = [list(np.flatnonzero(row <= tolerance)) for row in distances]
    assignment = [-1] * len(right)

    def augment(left_index: int, seen: set[int]) -> bool:
        for right_index in sorted(adjacency[left_index], key=lambda index: distances[left_index, index]):
            if right_index in seen:
                continue
            seen.add(right_index)
            if assignment[right_index] < 0 or augment(assignment[right_index], seen):
                assignment[right_index] = left_index
                return True
        return False

    if any(not augment(index, set()) for index in range(len(left))):
        nearest = np.min(distances, axis=1)
        raise ValueError(
            "no periodic atom bijection within tolerance; "
            f"largest nearest distance={float(np.max(nearest)):.12e} Angstrom"
        )
    matched = [float(distances[left_index, right_index]) for right_index, left_index in enumerate(assignment)]
    return max(matched, default=0.0), matched


def lattice_mappings(
    archived_cell: np.ndarray, current_cell: np.ndarray
) -> list[tuple[np.ndarray, np.ndarray, float, float]]:
    archived_metric = archived_cell @ archived_cell.T
    current_metric = current_cell @ current_cell.T
    integer_vectors = [
        np.array(vector, dtype=int)
        for vector in product(range(-3, 4), repeat=3)
        if vector != (0, 0, 0)
    ]
    candidates = []
    for current_row in current_cell:
        target = float(current_row @ current_row)
        rows = [
            vector
            for vector in integer_vectors
            if abs(float(vector @ archived_metric @ vector) - target)
            <= LATTICE_METRIC_TOLERANCE_ANGSTROM2
        ]
        if not rows:
            raise ValueError("no equivalent archived lattice vector found")
        candidates.append(rows)

    mappings = []
    for row_a, row_b, row_c in product(*candidates):
        transform = np.vstack((row_a, row_b, row_c))
        determinant = float(np.linalg.det(transform))
        if abs(abs(determinant) - 1.0) > 1.0e-8:
            continue
        transformed_metric = transform @ archived_metric @ transform.T
        metric_residual = float(np.max(np.abs(transformed_metric - current_metric)))
        if metric_residual > LATTICE_METRIC_TOLERANCE_ANGSTROM2:
            continue
        transformed_cell = transform @ archived_cell
        rotation = np.linalg.solve(transformed_cell, current_cell)
        orthogonality_residual = float(
            np.max(np.abs(rotation @ rotation.T - np.eye(3)))
        )
        mappings.append(
            (
                transform,
                rotation,
                metric_residual,
                orthogonality_residual,
            )
        )
    if not mappings:
        raise ValueError("cells are not related by a unimodular basis change and rotation")
    return sorted(mappings, key=lambda row: row[2] + row[3])


def structure_mapping(
    archived_cell: np.ndarray,
    archived_species: list[str],
    archived_coordinates: np.ndarray,
    current_cell: np.ndarray,
    current_species: list[str],
    current_coordinates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float, float, list[float]]:
    inverse = np.linalg.inv(current_cell)
    current_fractional = current_coordinates @ inverse
    anchor_species = min(Counter(archived_species), key=Counter(archived_species).get)
    archived_anchors = [
        index for index, value in enumerate(archived_species) if value == anchor_species
    ]
    current_anchors = [
        index for index, value in enumerate(current_species) if value == anchor_species
    ]
    best = None
    for transform, rotation, metric_residual, orthogonality_residual in lattice_mappings(
        archived_cell, current_cell
    ):
        mapped_fractional = (archived_coordinates @ rotation) @ inverse
        for archived_index, current_index in product(archived_anchors, current_anchors):
            translation = current_fractional[current_index] - mapped_fractional[archived_index]
            shifted_coordinates = (mapped_fractional + translation) @ current_cell
            matched_distances = []
            try:
                for species in sorted(set(archived_species)):
                    archived_indices = [
                        index for index, value in enumerate(archived_species) if value == species
                    ]
                    current_indices = [
                        index for index, value in enumerate(current_species) if value == species
                    ]
                    _, distances = perfect_species_match(
                        current_cell,
                        shifted_coordinates[archived_indices],
                        current_coordinates[current_indices],
                        MATCH_TOLERANCE_ANGSTROM,
                    )
                    matched_distances.extend(distances)
            except ValueError:
                continue
            maximum_distance = max(matched_distances, default=0.0)
            if best is None or maximum_distance < best[0]:
                best = (
                    maximum_distance,
                    transform,
                    rotation,
                    metric_residual,
                    orthogonality_residual,
                    translation - np.floor(translation),
                    matched_distances,
                )
    if best is None:
        raise ValueError("no global periodic translation produces an atom bijection")
    return (
        best[1],
        best[2],
        best[3],
        best[4],
        best[5],
        best[6],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("structure_root", type=Path)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    rows = []
    with zipfile.ZipFile(args.archive) as archive:
        cell_member = f"{ARCHIVE_ROOT}/cell_vectors_angstrom.csv"
        cell_data = archive.read(cell_member)
        cells = {}
        for row in csv.DictReader(io.StringIO(cell_data.decode("utf-8"))):
            phase = row["system"]
            cells[phase] = np.array(
                [
                    [float(row[f"{vector}_{axis}_angstrom"]) for axis in "xyz"]
                    for vector in "abc"
                ]
            )

        for phase in PHASES:
            xyz_member = f"{ARCHIVE_ROOT}/coordinates_angstrom/{phase}.xyz"
            xyz_data = archive.read(xyz_member)
            archived_species, archived_coordinates = parse_xyz(xyz_data)
            poscar = args.structure_root / phase / "POSCAR"
            current_cell, current_species, current_coordinates = parse_poscar(poscar)
            archived_cell = cells[phase]
            if Counter(archived_species) != Counter(current_species):
                raise ValueError(f"species-count mismatch for {phase}")
            volume_difference = abs(float(np.linalg.det(archived_cell)) - float(np.linalg.det(current_cell)))
            try:
                (
                    lattice_transform,
                    rotation,
                    metric_residual,
                    orthogonality_residual,
                    translation,
                    matched_distances,
                ) = structure_mapping(
                    archived_cell,
                    archived_species,
                    archived_coordinates,
                    current_cell,
                    current_species,
                    current_coordinates,
                )
            except ValueError as error:
                raise ValueError(f"{phase}: {error}") from error
            rows.append(
                {
                    "phase": phase,
                    "atom_count": len(archived_species),
                    "species_counts": dict(sorted(Counter(archived_species).items())),
                    "lattice_unimodular_transform": lattice_transform.tolist(),
                    "cartesian_rotation": rotation.tolist(),
                    "fractional_origin_translation": translation.tolist(),
                    "lattice_metric_max_abs_residual_angstrom2": metric_residual,
                    "rotation_orthogonality_max_abs_residual": orthogonality_residual,
                    "volume_abs_difference_angstrom3": volume_difference,
                    "maximum_matched_atom_distance_angstrom": max(matched_distances, default=0.0),
                    "archive_xyz_sha256": digest_bytes(xyz_data),
                    "production_poscar_sha256": digest_file(poscar),
                }
            )

    maximum_metric = max(row["lattice_metric_max_abs_residual_angstrom2"] for row in rows)
    maximum_orthogonality = max(row["rotation_orthogonality_max_abs_residual"] for row in rows)
    maximum_volume = max(row["volume_abs_difference_angstrom3"] for row in rows)
    maximum_atom = max(row["maximum_matched_atom_distance_angstrom"] for row in rows)
    payload = {
        "status": "PASS",
        "phase_count": len(rows),
        "archive_sha256": digest_file(args.archive),
        "archive_cell_table_sha256": digest_bytes(cell_data),
        "tolerances": {
            "lattice_metric_angstrom2": LATTICE_METRIC_TOLERANCE_ANGSTROM2,
            "matched_atom_angstrom": MATCH_TOLERANCE_ANGSTROM,
        },
        "summary": {
            "maximum_lattice_metric_abs_residual_angstrom2": maximum_metric,
            "maximum_rotation_orthogonality_abs_residual": maximum_orthogonality,
            "maximum_volume_abs_difference_angstrom3": maximum_volume,
            "maximum_matched_atom_distance_angstrom": maximum_atom,
        },
        "rows": rows,
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_json is not None:
        args.output_json.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
