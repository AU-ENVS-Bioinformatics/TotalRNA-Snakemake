#!/usr/bin/env python3
"""Create a TransDecoder protein-to-parent-transcript mapping.

Functionality
-------------
The script reads TransDecoder's GFF3 and extracts every ``mRNA`` feature. GFF3
column one supplies the original RNA-SPAdes transcript ID, while the feature's
``ID`` attribute supplies the TransDecoder protein/ORF identifier used by the
protein FASTA and eggNOG-mapper.

Input
-----
* TransDecoder ``ORF_genes.gff3``.

Outputs
-------
* ``protein_to_transcript.tsv`` with ``protein_id`` and ``contig_id`` columns.
* A log reporting the number of unique mappings.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Set, Tuple


def parse_gff3_attributes(value: str) -> Dict[str, str]:
    """Parse a semicolon-separated GFF3 attribute field."""

    attributes: Dict[str, str] = {}
    for entry in value.split(";"):
        entry = entry.strip()
        if "=" in entry:
            key, item = entry.split("=", 1)
            attributes[key] = item
    return attributes


def extract_protein_mapping(gff_path: str | Path) -> Set[Tuple[str, str]]:
    """Return unique ``(protein_id, contig_id)`` pairs from mRNA features."""

    mappings: Set[Tuple[str, str]] = set()
    with Path(gff_path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Invalid GFF3 record at line {line_number}: expected 9 columns"
                )
            if fields[2] != "mRNA":
                continue

            protein_id = parse_gff3_attributes(fields[8]).get("ID")
            if protein_id:
                mappings.add((protein_id, fields[0]))
    return mappings


def parse_arguments() -> argparse.Namespace:
    """Parse command-line paths for protein mapping."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", required=True, help="TransDecoder GFF3 input")
    parser.add_argument("--output", required=True, help="Output mapping TSV")
    parser.add_argument("--log", required=True, help="Output log file")
    return parser.parse_args()


def main() -> None:
    """Extract, sort, and write TransDecoder protein mappings."""

    args = parse_arguments()
    mappings = sorted(extract_protein_mapping(args.gff))

    with Path(args.output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "contig_id"])
        writer.writerows(mappings)

    with Path(args.log).open("w") as handle:
        handle.write(f"Mapped proteins: {len(mappings)}\n")


if __name__ == "__main__":
    main()
