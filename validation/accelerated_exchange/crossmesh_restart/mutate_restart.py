#!/usr/bin/env python3
"""Create targeted corruptions of a CP2K sequential-unformatted restart."""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path


def records(blob: bytearray) -> list[tuple[int, int, int]]:
    result: list[tuple[int, int, int]] = []
    pos = 0
    while pos < len(blob):
        if pos + 4 > len(blob):
            raise ValueError("truncated leading marker")
        size = struct.unpack_from("<i", blob, pos)[0]
        data = pos + 4
        end = data + size
        if size < 0 or end + 4 > len(blob):
            raise ValueError("invalid record size")
        if struct.unpack_from("<i", blob, end)[0] != size:
            raise ValueError("record marker mismatch")
        result.append((pos, data, size))
        pos = end + 4
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument(
        "mode",
        choices=("nan_shift", "nan_hmat", "nan_checksum", "huge_payload", "append_record"),
    )
    args = parser.parse_args()

    blob = bytearray(args.source.read_bytes())
    recs = records(blob)
    magic = b"CP2K_KPOINT_RESTART_V2".ljust(32)
    magic_index = next(i for i, (_, data, size) in enumerate(recs) if size == 32 and blob[data : data + size] == magic)

    if args.mode == "nan_shift":
        _, data, size = recs[magic_index + 2]
        assert size == 40
        struct.pack_into("<d", blob, data + 12, math.nan)
    elif args.mode == "nan_hmat":
        _, data, size = recs[magic_index + 5]
        assert size == 72
        struct.pack_into("<d", blob, data, math.nan)
    elif args.mode == "nan_checksum":
        _, data, size = recs[magic_index + 10]
        assert size == 20
        struct.pack_into("<d", blob, data + 4, math.nan)
    elif args.mode == "huge_payload":
        _, header_data, header_size = recs[magic_index + 11]
        assert header_size == 20
        _, data, size = recs[magic_index + 12]
        assert size % 8 == 0
        struct.pack_into("<d", blob, data, 1.0e200)
    elif args.mode == "append_record":
        payload = b"TRAILING_GARBAGE"
        marker = struct.pack("<i", len(payload))
        blob.extend(marker + payload + marker)

    args.target.write_bytes(blob)


if __name__ == "__main__":
    main()
