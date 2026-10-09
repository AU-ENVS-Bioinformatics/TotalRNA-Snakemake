#!/usr/bin/env python3
"""Extract eggNOG annotations belonging to geNomad-positive transcripts.

Functionality
-------------
The script uses the TransDecoder protein-to-transcript map to identify proteins
encoded by viral parent transcripts, then selects their existing eggNOG rows.
It adds explicit ``contig_id`` and ``protein_id`` columns to the result.

Inputs
------
* ``viral_candidate_ids.txt``.
* ``protein_to_transcript.tsv``.
* ``eggnog.emapper.annotations``.

Outputs
-------
* ``viral_eggnog.annotations.tsv``.
* A log with selected transcript, protein, and annotation counts.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple


MISSING = "NA"


def read_candidate_ids(path: str | Path) -> Set[str]:
    """Read viral parent transcript identifiers."""

    with Path(path).open() as handle:
        return {line.strip() for line in handle if line.strip()}


def read_protein_map(path: str | Path) -> Dict[str, str]:
    """Read and validate the protein-to-transcript mapping."""

    mapping: Dict[str, str] = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"protein_id", "contig_id"}.issubset(
            reader.fieldnames
        ):
            raise ValueError("Mapping TSV requires protein_id and contig_id columns")

        for row in reader:
            protein_id = row["protein_id"]
            contig_id = row["contig_id"]
            previous = mapping.get(protein_id)
            if previous is not None and previous != contig_id:
                raise ValueError(
                    f"Protein '{protein_id}' maps to both '{previous}' and "
                    f"'{contig_id}'"
                )
            mapping[protein_id] = contig_id
    return mapping


def read_eggnog(path: str | Path) -> Tuple[List[str], List[List[str]]]:
    """Read the commented eggNOG header and all annotation rows."""

    header: List[str] | None = None
    rows: List[List[str]] = []

    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#"):
                candidate = line.lstrip("#").rstrip("\n").split("\t")
                if candidate and candidate[0] in {"query", "query_name"}:
                    header = candidate
                continue
            if line.strip():
                rows.append(line.rstrip("\n").split("\t"))

    if header is None:
        raise ValueError("Could not find the '#query' eggNOG header")
    return header, rows


def parse_arguments() -> argparse.Namespace:
    """Parse viral ID, mapping, eggNOG, and output paths."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--eggnog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main() -> None:
    """Select viral eggNOG rows and write the standardized annotation table."""

    args = parse_arguments()
    candidate_ids = read_candidate_ids(args.ids)
    protein_map = read_protein_map(args.mapping)
    eggnog_header, eggnog_rows = read_eggnog(args.eggnog)

    output_header = ["contig_id", "protein_id", *eggnog_header[1:]]
    selected_rows: List[List[str]] = []
    selected_proteins: Set[str] = set()

    for values in eggnog_rows:
        values += [MISSING] * (len(eggnog_header) - len(values))
        protein_id = values[0]
        contig_id = protein_map.get(protein_id)
        if contig_id not in candidate_ids:
            continue
        selected_rows.append([contig_id, protein_id, *values[1:]])
        selected_proteins.add(protein_id)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_header)
        writer.writerows(selected_rows)

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.log).open("w") as handle:
        handle.write(f"Viral parent transcripts: {len(candidate_ids)}\n")
        handle.write(f"Viral proteins with eggNOG rows: {len(selected_proteins)}\n")
        handle.write(f"Viral eggNOG annotation rows: {len(selected_rows)}\n")


if __name__ == "__main__":
    main()
