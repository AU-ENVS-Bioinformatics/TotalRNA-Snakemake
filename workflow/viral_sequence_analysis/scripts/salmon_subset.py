#!/usr/bin/env python3
"""Subset an existing transcript-level Salmon quantification to viral contigs.

Functionality
-------------
The script retains rows from a Salmon ``quant.sf`` whose ``Name`` occurs in the
viral parent-contig ID list. It does not recalculate abundance, so TPM and
NumReads remain directly comparable with the full transcriptome quantification.

Inputs
------
* ``viral_candidate_ids.txt`` containing one transcript ID per line.
* A transcript-reference Salmon ``quant.sf`` file.

Outputs
-------
* ``viral_quant.sf`` with the original Salmon columns and viral rows only.
* A log reporting matched and unmatched candidate counts.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set


def read_candidate_ids(path: str | Path) -> Set[str]:
    """Read one candidate identifier per line into a set."""

    with Path(path).open() as handle:
        return {line.strip() for line in handle if line.strip()}


def subset_quantification(
    quant_path: str | Path,
    candidate_ids: Set[str],
    output_path: str | Path,
) -> Set[str]:
    """Write viral Salmon rows and return the identifiers that were matched."""

    matched: Set[str] = set()

    with Path(quant_path).open(newline="") as source, Path(output_path).open(
        "w", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None or "Name" not in reader.fieldnames:
            raise ValueError("Salmon quant.sf lacks the required 'Name' column")

        writer = csv.DictWriter(
            destination,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            sequence_id = row["Name"]
            if sequence_id in candidate_ids:
                writer.writerow(row)
                matched.add(sequence_id)

    return matched


def parse_arguments() -> argparse.Namespace:
    """Parse command-line paths for Salmon subsetting."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", required=True, help="Viral candidate ID list")
    parser.add_argument("--quant", required=True, help="Input Salmon quant.sf")
    parser.add_argument("--output", required=True, help="Output viral quant.sf")
    parser.add_argument("--log", required=True, help="Output log file")
    return parser.parse_args()


def main() -> None:
    """Subset Salmon and record identifier matching statistics."""

    args = parse_arguments()
    candidate_ids = read_candidate_ids(args.ids)
    matched_ids = subset_quantification(args.quant, candidate_ids, args.output)

    with Path(args.log).open("w") as handle:
        handle.write(f"Viral candidate IDs: {len(candidate_ids)}\n")
        handle.write(f"Candidates found in Salmon: {len(matched_ids)}\n")
        handle.write(
            f"Candidates absent from Salmon: {len(candidate_ids - matched_ids)}\n"
        )

if __name__ == "__main__":
    main()
