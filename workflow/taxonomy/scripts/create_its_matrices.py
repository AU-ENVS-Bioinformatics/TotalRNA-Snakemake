#!/usr/bin/env python3
"""Create cross-sample count matrices from combined ITS assignment tables."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


VALID_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", required=True, type=Path)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--rank", choices=VALID_RANKS, default="genus")
    parser.add_argument("--assigned-output", required=True, type=Path)
    parser.add_argument("--unassigned-output", required=True, type=Path)
    parser.add_argument("--status-output", required=True, type=Path)
    args = parser.parse_args()
    if len(args.tables) != len(args.samples):
        parser.error("--tables and --samples must contain the same number of values")
    if len(args.samples) != len(set(args.samples)):
        parser.error("--samples contains duplicate sample names")
    return args


def read_rows(path: Path, expected_sample: str) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "status", "classification_reason"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing required columns in {path}")
        rows = list(reader)
    for row in rows:
        if row["sample"] != expected_sample:
            raise ValueError(
                f"Table {path} contains sample '{row['sample']}', expected '{expected_sample}'"
            )
    return rows


def write_matrix(path: Path, row_header: str, samples: list[str], counts: dict[str, Counter[str]]) -> None:
    labels = sorted({label for sample_counts in counts.values() for label in sample_counts})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_header, *samples])
        for label in labels:
            writer.writerow([label, *(counts[sample][label] for sample in samples)])


def main() -> None:
    args = arguments()
    assigned: dict[str, Counter[str]] = defaultdict(Counter)
    unassigned: dict[str, Counter[str]] = defaultdict(Counter)
    statuses: dict[str, Counter[str]] = defaultdict(Counter)
    for sample in args.samples:
        assigned[sample]
        unassigned[sample]
        statuses[sample]

    for path, sample in zip(args.tables, args.samples):
        rows = read_rows(path, sample)
        if rows and args.rank not in rows[0]:
            raise ValueError(f"Missing rank column '{args.rank}' in {path}")
        for row in rows:
            status = row["status"] or "unassigned"
            statuses[sample][status] += 1
            if status == "assigned":
                taxon = row.get(args.rank, "") or f"Unclassified_{args.rank}"
                assigned[sample][taxon] += 1
            else:
                reason = row["classification_reason"] or "unspecified_reason"
                unassigned[sample][reason] += 1

    write_matrix(args.assigned_output, "taxon", args.samples, assigned)
    write_matrix(args.unassigned_output, "classification_reason", args.samples, unassigned)
    write_matrix(args.status_output, "status", args.samples, statuses)


if __name__ == "__main__":
    main()
