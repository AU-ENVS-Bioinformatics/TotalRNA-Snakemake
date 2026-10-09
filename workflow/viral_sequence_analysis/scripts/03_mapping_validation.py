#!/usr/bin/env python3
"""Validate overlap between TransDecoder mapping IDs and eggNOG query IDs.

Functionality
-------------
The script reproduces the logic of sorting both ID sets, applying ``comm -12``,
and previewing the first ten shared IDs. It fails when no identifiers overlap,
preventing silent production of empty viral annotation tables.

Inputs
------
* ``protein_to_transcript.tsv``.
* ``eggnog.emapper.annotations``.

Outputs
-------
* All overlapping protein IDs.
* A ten-ID preview file.
* A TSV summary with mapping counts and overlap percentage.
* A validation log.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Set


def read_mapping_ids(path: str | Path) -> Set[str]:
    """Read protein identifiers from the mapping TSV."""

    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "protein_id" not in reader.fieldnames:
            raise ValueError("Mapping TSV lacks the required 'protein_id' column")
        return {row["protein_id"] for row in reader if row.get("protein_id")}


def read_eggnog_ids(path: str | Path) -> Set[str]:
    """Read query identifiers from non-comment eggNOG annotation rows."""

    identifiers: Set[str] = set()
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            identifiers.add(line.split("\t", 1)[0])
    return identifiers


def write_ids(path: str | Path, identifiers: list[str]) -> None:
    """Write one protein identifier per line."""

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w") as handle:
        for protein_id in identifiers:
            handle.write(f"{protein_id}\n")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line inputs and validation output paths."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--eggnog", required=True)
    parser.add_argument("--overlap", required=True)
    parser.add_argument("--preview", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> None:
    """Calculate ID overlap, write validation outputs, and enforce a match."""

    args = parse_arguments()
    mapping_ids = read_mapping_ids(args.mapping)
    eggnog_ids = read_eggnog_ids(args.eggnog)
    overlapping_ids = sorted(mapping_ids & eggnog_ids)

    write_ids(args.overlap, overlapping_ids)
    write_ids(args.preview, overlapping_ids[:10])

    overlap_percentage = (
        100.0 * len(overlapping_ids) / len(mapping_ids) if mapping_ids else 0.0
    )
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.summary).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["transdecoder_mapped_proteins", len(mapping_ids)])
        writer.writerow(["eggnog_annotated_proteins", len(eggnog_ids)])
        writer.writerow(["overlapping_protein_ids", len(overlapping_ids)])
        writer.writerow(["mapping_overlap_percent", f"{overlap_percentage:.4f}"])

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.log).open("w") as handle:
        handle.write(f"TransDecoder mapped proteins: {len(mapping_ids)}\n")
        handle.write(f"eggNOG annotated proteins: {len(eggnog_ids)}\n")
        handle.write(f"Overlapping protein IDs: {len(overlapping_ids)}\n")

    if not overlapping_ids:
        raise ValueError(
            "No TransDecoder protein IDs overlap eggNOG query IDs; inspect "
            "the GFF3 and eggNOG identifier formats"
        )


if __name__ == "__main__":
    main()
