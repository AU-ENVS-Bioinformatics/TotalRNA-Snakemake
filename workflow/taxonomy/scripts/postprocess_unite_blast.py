#!/usr/bin/env python3
"""Filter one ITS-region BLAST table and assign UNITE taxonomy by LCA."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "qlen", "sstart", "send", "slen", "qcovs",
    "qcovhsp", "evalue", "bitscore", "stitle",
]
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
PREFIX_TO_RANK = {
    "k__": "kingdom", "p__": "phylum", "c__": "class", "o__": "order",
    "f__": "family", "g__": "genus", "s__": "species",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter BLASTN hits by identity and query coverage, retain hits near "
            "the best bit score, and calculate their lowest common taxonomy."
        )
    )
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-query-coverage", type=float, default=80.0)
    parser.add_argument("--top-bitscore-fraction", type=float, default=0.99)
    parser.add_argument(
        "--query-fasta",
        type=Path,
        help="Optional ITSx FASTA used to retain queries without BLAST hits.",
    )
    args = parser.parse_args()
    if not 0 < args.top_bitscore_fraction <= 1:
        parser.error("--top-bitscore-fraction must be greater than 0 and at most 1")
    return args


def fasta_ids(path: Path | None) -> list[str]:
    if path is None:
        return []
    identifiers = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].strip().split()[0])
    return identifiers


def parse_taxonomy(subject: str) -> dict[str, str]:
    """Parse pipe metadata and semicolon-delimited UNITE lineage fields."""
    taxonomy = {rank: "" for rank in RANKS}
    fields = []
    for pipe_field in subject.split("|"):
        fields.extend(pipe_field.split(";"))
    for field in fields:
        field = field.strip()
        for prefix, rank in PREFIX_TO_RANK.items():
            if field.startswith(prefix):
                taxonomy[rank] = field[len(prefix):].strip()
                break
    return taxonomy


def common_taxonomy(hits: list[dict]) -> tuple[dict[str, str], str, str]:
    consensus = {rank: "" for rank in RANKS}
    assigned_rank = "unassigned"
    assigned_name = "unassigned"
    for rank in RANKS:
        values = [hit["taxonomy"].get(rank, "") for hit in hits]
        if any(not value for value in values) or len(set(values)) != 1:
            break
        value = values[0]
        consensus[rank] = value
        assigned_rank = rank
        assigned_name = value
    return consensus, assigned_rank, assigned_name


def read_hits(path: Path) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    if not path.exists() or path.stat().st_size == 0:
        return grouped
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=COLUMNS)
        for line_number, row in enumerate(reader, start=1):
            if None in row or any(row[name] is None for name in COLUMNS):
                raise ValueError(f"Malformed BLAST row {line_number} in {path}")
            try:
                row["pident"] = float(row["pident"])
                row["qcovs"] = float(row["qcovs"])
                row["qcovhsp"] = float(row["qcovhsp"])
                row["bitscore"] = float(row["bitscore"])
                row["evalue"] = float(row["evalue"])
            except ValueError as error:
                raise ValueError(f"Invalid numeric value on BLAST row {line_number}") from error
            row["taxonomy"] = parse_taxonomy(row["sseqid"])
            grouped[row["qseqid"]].append(row)
    return grouped


def fill_taxonomy(
    taxonomy: dict[str, str],
    assigned_rank: str,
    assigned_name: str,
) -> dict[str, str]:
    """Ensure that no taxonomy column is empty in the output TSV."""
    if assigned_rank not in RANKS:
        return {rank: "Unassigned" for rank in RANKS}

    assigned_index = RANKS.index(assigned_rank)
    unclassified_label = f"Unclassified_{assigned_name}"
    return {
        rank: taxonomy.get(rank, "")
        or (unclassified_label if index > assigned_index else "Unassigned")
        for index, rank in enumerate(RANKS)
    }


def main() -> None:
    args = arguments()
    grouped = read_hits(args.blast)
    query_ids = list(dict.fromkeys(fasta_ids(args.query_fasta) + list(grouped)))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "sample", "region", "query_id", "status", "classification_reason", "assigned_rank",
        "assigned_name", *RANKS, "n_raw_hits", "n_passing_hits",
        "n_near_top_hits", "best_subject", "best_identity",
        "best_query_coverage", "best_evalue", "best_bitscore",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for query_id in query_ids:
            raw = grouped.get(query_id, [])
            passing = [
                hit for hit in raw
                if hit["pident"] >= args.min_identity
                and hit["qcovs"] >= args.min_query_coverage
            ]
            result = {
                "sample": args.sample,
                "region": args.region,
                "query_id": query_id,
                "status": "no_blast_hits" if not raw else "no_passing_hits",
                "classification_reason": "no_blast_hits" if not raw else "no_hit_passed_both_thresholds",
                "assigned_rank": "unassigned",
                "assigned_name": "unassigned",
                **{rank: "" for rank in RANKS},
                "n_raw_hits": len(raw),
                "n_passing_hits": len(passing),
                "n_near_top_hits": 0,
                "best_subject": "",
                "best_identity": "",
                "best_query_coverage": "",
                "best_evalue": "",
                "best_bitscore": "",
            }
            if raw and not passing:
                identity_passes = [
                    hit["pident"] >= args.min_identity
                    for hit in raw
                ]
                coverage_passes = [
                    hit["qcovs"] >= args.min_query_coverage
                    for hit in raw
                ]
                if not any(identity_passes) and not any(coverage_passes):
                    result["classification_reason"] = "hits_below_identity_and_query_coverage"
                elif not any(identity_passes):
                    result["classification_reason"] = "hits_below_identity"
                elif not any(coverage_passes):
                    result["classification_reason"] = "hits_below_query_coverage"
            if passing:
                best_score = max(hit["bitscore"] for hit in passing)
                near_top = [
                    hit for hit in passing
                    if hit["bitscore"] >= best_score * args.top_bitscore_fraction
                ]
                near_top.sort(key=lambda hit: (-hit["bitscore"], hit["evalue"], -hit["pident"]))
                best = near_top[0]
                consensus, assigned_rank, assigned_name = common_taxonomy(near_top)
                result.update(consensus)
                result.update({
                    "status": "assigned" if assigned_rank != "unassigned" else "ambiguous",
                    "classification_reason": (
                        "passed_filters_and_assigned"
                        if assigned_rank != "unassigned"
                        else "near_top_hits_taxonomically_ambiguous"
                    ),
                    "assigned_rank": assigned_rank,
                    "assigned_name": assigned_name,
                    "n_near_top_hits": len(near_top),
                    "best_subject": best["sseqid"],
                    "best_identity": best["pident"],
                    "best_query_coverage": best["qcovs"],
                    "best_evalue": best["evalue"],
                    "best_bitscore": best["bitscore"],
                })

            result.update(
                fill_taxonomy(
                    taxonomy={rank: result[rank] for rank in RANKS},
                    assigned_rank=result["assigned_rank"],
                    assigned_name=result["assigned_name"],
                )
            )

            for field in [
                "best_subject",
                "best_identity",
                "best_query_coverage",
                "best_evalue",
                "best_bitscore",
            ]:
                if result[field] == "":
                    result[field] = "None"

            writer.writerow(result)


if __name__ == "__main__":
    main()
