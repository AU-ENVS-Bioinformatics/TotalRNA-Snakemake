#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Set


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# -------------------------
# CIGAR parsing
# -------------------------
def cigar_aligned_bases(cigar: str) -> int:
    return sum(
        int(length)
        for length, op in CIGAR_RE.findall(cigar)
        if op in {"M", "=", "X"}
    )


# -------------------------
# NM extraction
# -------------------------
def extract_nm(fields: list[str]) -> int:
    for f in fields:
        if f.startswith("NM:i:"):
            return int(f.split(":")[-1])
    return 0


# -------------------------
# Mate detection (robust)
# -------------------------
def get_mate_from_flag(flag: int) -> str:
    if flag & 64:
        return "R1"
    if flag & 128:
        return "R2"
    return "U"


# -------------------------
# PASS logic
# -------------------------
def passes(fraction: float, nm: int, min_fraction: float, max_nm: int | None) -> bool:
    if fraction < min_fraction:
        return False
    if max_nm is not None and nm > max_nm:
        return False
    return True


# -------------------------
# MAIN PROCESS
# -------------------------
def process(
    input_path: Path,
    min_fraction: float,
    max_nm: int | None,
    prefix: str,
    read_out_path: Path | None,
    sam_out_path: Path | None,
) -> None:

    read_map: Dict[str, Dict[str, bool]] = {}

    def init_read(read_id: str) -> None:
        if read_id not in read_map:
            read_map[read_id] = {
                "prefix_R1": False,
                "prefix_R2": False,
                "pass": False,
            }

    headers: list[str] = []
    all_lines: list[str] = []

    # -------------------------
    # READ SAM (robust)
    # -------------------------
    with open(input_path) as fh:
        for i, line in enumerate(fh, 1):

            if line.startswith("@"):
                headers.append(line)
                all_lines.append(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            read_id = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            cigar = fields[5]
            seq = fields[9]

            # ✅ Track all alignments for later writing
            all_lines.append(line)

            # prefix filtering
            if not rname.startswith(prefix):
                continue

            init_read(read_id)

            # detect mate
            mate = get_mate_from_flag(flag)

            # fallback: SortMeRNA → treat as both mates
            if mate == "U":
                mates = ["R1", "R2"]
            else:
                mates = [mate]

            # mark prefix
            for m in mates:
                read_map[read_id][f"prefix_{m}"] = True

            # compute alignment stats
            if not seq:
                continue

            aligned = cigar_aligned_bases(cigar)
            fraction = aligned / len(seq)

            nm = extract_nm(fields)

            if passes(fraction, nm, min_fraction, max_nm):
                read_map[read_id]["pass"] = True

            if i % 100000 == 0:
                print(f"[INFO] Processed {i:,} lines")

    print(f"[INFO] Total reads seen: {len(read_map):,}")

    # -------------------------
    # FILTER READS
    # -------------------------
    passing_reads: Set[str] = set()

    for read_id, flags in read_map.items():

        # must map to prefix
        if not (flags["prefix_R1"] or flags["prefix_R2"]):
            continue

        # disallow mixed SSU/LSU assignments
        if flags["prefix_R1"] != flags["prefix_R2"]:
            continue

        # must pass at least once
        if not flags["pass"]:
            continue

        passing_reads.add(read_id)

    print(f"[INFO] Passing reads: {len(passing_reads):,}")

    # -------------------------
    # OUTPUT READ IDS
    # -------------------------
    if read_out_path:
        with open(read_out_path, "w") as out:
            for rid in sorted(passing_reads):
                out.write(rid + "\n")

    # -------------------------
    # OUTPUT SAM FILE
    # -------------------------
    if sam_out_path:
        print(f"[INFO] Writing filtered SAM: {sam_out_path}")

        with open(sam_out_path, "w") as out:

            # write header
            for h in headers:
                out.write(h)

            # write only passing reads
            for line in all_lines:
                if line.startswith("@"):
                    continue
                read_id = line.split("\t")[0]
                if read_id in passing_reads:
                    out.write(line)


# -------------------------
# CLI
# -------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()

    p.add_argument("--sam_in", required=True)
    p.add_argument("--fraction", type=float, required=True)
    p.add_argument("--matches", type=int, default=None)
    p.add_argument("--prefix", required=True)

    p.add_argument("--read_out", default=None)
    p.add_argument("--sam_out", default=None)

    return p


def main() -> None:
    args = build_parser().parse_args()

    process(
        input_path=Path(args.sam_in),
        min_fraction=args.fraction,
        max_nm=args.matches,
        prefix=args.prefix,
        read_out_path=Path(args.read_out) if args.read_out else None,
        sam_out_path=Path(args.sam_out) if args.sam_out else None,
    )


if __name__ == "__main__":
    main()