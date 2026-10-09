#!/usr/bin/env python3
"""Create the final viral candidate table from geNomad and CheckV.

Functionality
-------------
The script joins geNomad's virus summary with CheckV's quality summary. It also
creates ``parent_contig_id`` by removing geNomad's provirus-coordinate suffix,
allowing the candidates to be linked to the original RNA-SPAdes transcripts.

Inputs
------
* A geNomad ``*_virus_summary.tsv`` file.
* A CheckV ``quality_summary.tsv`` file.

Outputs
-------
* ``viral_candidates.tsv``: one row per geNomad viral candidate.
* ``viral_candidate_ids.txt``: unique parent transcript IDs, one per line.
* A log containing candidate and CheckV match counts.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


MISSING = "NA"

def read_tsv(path: str | Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read a headered TSV file and return its columns and rows."""

    input_path = Path(path)
    if not input_path.exists() or input_path.stat().st_size == 0:
        return [], []

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [], []
        return list(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
    path: str | Path,
    columns: Sequence[str],
    rows: Iterable[Mapping[str, str]],
) -> None:
    """Write dictionaries to a headered TSV file."""

    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, MISSING) for column in columns})


def get_parent_contig_id(sequence_id: str) -> str:
    """Remove a terminal ``|provirus_START_END`` suffix from an ID."""

    return re.sub(r"\|provirus_\d+_\d+$", "", sequence_id)


def build_candidate_table(
    genomad_path: str | Path,
    checkv_path: str | Path,
) -> Tuple[List[str], List[Dict[str, str]]]:
    """Join geNomad candidates with corresponding CheckV quality rows."""

    genomad_columns, genomad_rows = read_tsv(genomad_path)
    if "seq_name" not in genomad_columns:
        raise ValueError("geNomad summary lacks the required 'seq_name' column")

    checkv_columns, checkv_rows = read_tsv(checkv_path)
    checkv_by_id = {
        row["contig_id"]: row
        for row in checkv_rows
        if row.get("contig_id")
    }

    renamed_genomad_columns = [
        "genomad_taxonomy" if column == "taxonomy" else column
        for column in genomad_columns
    ]
    output_columns = [*renamed_genomad_columns, "parent_contig_id"]

    checkv_output_names: Dict[str, str] = {}
    for column in checkv_columns:
        if column == "contig_id":
            continue
        output_name = (
            f"{column}_checkv" if column in output_columns else column
        )
        checkv_output_names[column] = output_name
        output_columns.append(output_name)

    output_rows: List[Dict[str, str]] = []
    for genomad_row in genomad_rows:
        sequence_id = genomad_row["seq_name"]
        combined = {
            ("genomad_taxonomy" if key == "taxonomy" else key): value
            for key, value in genomad_row.items()
        }
        combined["parent_contig_id"] = get_parent_contig_id(sequence_id)

        checkv_row = checkv_by_id.get(sequence_id, {})
        for source_name, output_name in checkv_output_names.items():
            combined[output_name] = checkv_row.get(source_name, MISSING)
        output_rows.append(combined)

    return output_columns, output_rows


def parse_arguments() -> argparse.Namespace:
    """Parse command-line paths for the candidate-table workflow."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomad", required=True, help="geNomad virus summary TSV")
    parser.add_argument("--checkv", required=True, help="CheckV quality summary TSV")
    parser.add_argument("--output", required=True, help="Combined candidate TSV")
    parser.add_argument("--ids", required=True, help="Output parent-contig ID list")
    parser.add_argument("--log", required=True, help="Output log file")
    return parser.parse_args()


def main() -> None:
    """Run candidate-table construction and write all declared outputs."""

    args = parse_arguments()
    columns, rows = build_candidate_table(args.genomad, args.checkv)
    write_tsv(args.output, columns, rows)

    parent_ids = sorted(
        {row["parent_contig_id"] for row in rows if row.get("parent_contig_id")}
    )
    with Path(args.ids).open("w") as handle:
        for contig_id in parent_ids:
            handle.write(f"{contig_id}\n")

    checkv_matches = sum(
        row.get("checkv_quality", MISSING) != MISSING for row in rows
    )
    with Path(args.log).open("w") as handle:
        handle.write(f"geNomad candidates: {len(rows)}\n")
        handle.write(f"Unique parent RNA-SPAdes transcripts: {len(parent_ids)}\n")
        handle.write(f"Candidates with CheckV quality: {checkv_matches}\n")


if __name__ == "__main__":
    main()
