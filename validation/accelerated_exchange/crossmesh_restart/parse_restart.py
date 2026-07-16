#!/usr/bin/env python3
"""Read and compare self-describing CP2K k-point density restart files.

This intentionally understands only the version-2 records needed by the
independent cross-mesh restart audit.  It rejects malformed record markers,
duplicate payload keys, and trailing partial records.
"""

from __future__ import annotations

import argparse
import pathlib
import struct

import numpy as np


MAGIC = b"CP2K_KPOINT_RESTART_V2"


def records(path: pathlib.Path) -> list[bytes]:
    raw = path.read_bytes()
    out: list[bytes] = []
    pos = 0
    while pos < len(raw):
        if pos + 4 > len(raw):
            raise ValueError(f"partial record marker at byte {pos}")
        (size,) = struct.unpack_from("=i", raw, pos)
        pos += 4
        if size < 0 or pos + size + 4 > len(raw):
            raise ValueError(f"invalid record length {size} at byte {pos - 4}")
        body = raw[pos : pos + size]
        pos += size
        (tail,) = struct.unpack_from("=i", raw, pos)
        pos += 4
        if tail != size:
            raise ValueError(f"record marker mismatch {size} != {tail}")
        out.append(body)
    return out


def parse(path: pathlib.Path) -> dict:
    rec = records(path)
    if len(rec) < 3:
        raise ValueError("too few records")
    version = struct.unpack("=i", rec[0])[0]
    natom, nao, nset_max, nshell_max = struct.unpack("=4i", rec[1])
    magic_indices = [i for i, value in enumerate(rec) if value.rstrip(b" \x00") == MAGIC]
    if len(magic_indices) != 1:
        raise ValueError(f"expected one magic record, got {magic_indices}")
    m = magic_indices[0]
    if m + 11 >= len(rec):
        raise ValueError("truncated metadata")
    scheme = rec[m + 1].decode("ascii", errors="replace").rstrip(" \x00")
    grid = np.frombuffer(rec[m + 2][:12], dtype=np.dtype("=i4")).astype(int)
    shift = np.frombuffer(rec[m + 2][12:36], dtype=np.dtype("=f8")).copy()
    nrecords = struct.unpack_from("=i", rec[m + 10], 0)[0]
    payload: dict[tuple[int, tuple[int, int, int]], np.ndarray] = {}
    pos = m + 11
    for _ in range(nrecords * struct.unpack_from("=i", rec[m + 8], 8)[0]):
        if pos >= len(rec):
            raise ValueError("truncated payload header")
        ispin, irecord, cx, cy, cz = struct.unpack("=5i", rec[pos])
        pos += 1
        real = np.empty((nao, nao), dtype=float)
        imag = np.empty((nao, nao), dtype=float)
        for j in range(nao):
            if pos >= len(rec):
                raise ValueError("truncated real payload")
            real[:, j] = np.frombuffer(rec[pos], dtype=np.dtype("=f8"))
            pos += 1
        for j in range(nao):
            if pos >= len(rec):
                raise ValueError("truncated imaginary payload")
            imag[:, j] = np.frombuffer(rec[pos], dtype=np.dtype("=f8"))
            pos += 1
        key = (ispin, (cx, cy, cz))
        if key in payload:
            raise ValueError(f"duplicate payload key {key}")
        payload[key] = real + 1j * imag
    if pos != len(rec):
        raise ValueError(f"unexpected records after payload: {len(rec) - pos}")
    charge, multiplicity, nspin, ne_a, ne_b = struct.unpack("=5i", rec[m + 8])
    return {
        "path": path,
        "version": version,
        "natom": natom,
        "nao": nao,
        "nset_max": nset_max,
        "nshell_max": nshell_max,
        "scheme": scheme,
        "grid": grid,
        "shift": shift,
        "charge": charge,
        "multiplicity": multiplicity,
        "nspin": nspin,
        "nelectron": (ne_a, ne_b),
        "payload": payload,
    }


def compare(left: dict, right: dict) -> None:
    fields = ("version", "natom", "nao", "scheme", "charge", "multiplicity", "nspin", "nelectron")
    for field in fields:
        if left[field] != right[field]:
            raise ValueError(f"metadata mismatch {field}: {left[field]} != {right[field]}")
    if not np.array_equal(left["grid"], right["grid"]):
        raise ValueError(f"grid mismatch {left['grid']} != {right['grid']}")
    keys_left = set(left["payload"])
    keys_right = set(right["payload"])
    if keys_left != keys_right:
        raise ValueError(f"payload-key mismatch: left-only={keys_left - keys_right}, right-only={keys_right - keys_left}")
    max_abs = 0.0
    sum_sq = 0.0
    count = 0
    ref_sq = 0.0
    for key in sorted(keys_left):
        delta = left["payload"][key] - right["payload"][key]
        max_abs = max(max_abs, float(np.max(np.abs(delta))))
        sum_sq += float(np.vdot(delta, delta).real)
        ref_sq += float(np.vdot(right["payload"][key], right["payload"][key]).real)
        count += delta.size
    rms = (sum_sq / count) ** 0.5
    rel_fro = sum_sq**0.5 / max(ref_sq**0.5, np.finfo(float).tiny)
    print(f"left={left['path']}")
    print(f"right={right['path']}")
    print(f"grid={'x'.join(map(str, left['grid']))} nspin={left['nspin']} nao={left['nao']} records={len(keys_left)}")
    print(f"max_abs={max_abs:.17e}")
    print(f"rms={rms:.17e}")
    print(f"relative_frobenius={rel_fro:.17e}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=pathlib.Path)
    parser.add_argument("right", type=pathlib.Path)
    args = parser.parse_args()
    compare(parse(args.left), parse(args.right))


if __name__ == "__main__":
    main()
