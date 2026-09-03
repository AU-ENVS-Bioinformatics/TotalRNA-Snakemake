#!/usr/bin/env python3
"""Combine full, ITS1, ITS2 and 5.8S UNITE assignments per parent contig."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
INFORMATIVE = ["full", "ITS2", "ITS1"]


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignments", nargs="+", required=True, type=Path)
    parser.add_argument("--assigned-output", required=True, type=Path)
    parser.add_argument("--unassigned-output", required=True, type=Path)
    parser.add_argument("--all-output", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--assembly-fasta",
        required=True,
        type=Path,
        help="Original rnaSPAdes ITS-candidate assembly before ITSx.",
    )
    return parser.parse_args()


def fasta_ids(path: Path) -> list[str]:
    identifiers = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].strip().split()[0])
    return identifiers


def parent_id(query_id: str) -> str:
    return query_id.split("|F", 1)[0]


def depth(row: dict) -> int:
    rank = row.get("assigned_rank", "unassigned")
    return RANKS.index(rank) if rank in RANKS else -1


def metric(row: dict, name: str) -> float:
    try:
        return float(row.get(name, ""))
    except (TypeError, ValueError):
        return float("-inf")


def is_placeholder(rank: str, value: str) -> bool:
    """Return True for UNITE labels that do not resolve the stated rank."""
    lowered = value.lower()
    if not value or "incertae_sedis" in lowered or lowered.startswith("unclassified"):
        return True
    if rank == "species" and (lowered.endswith("_sp") or lowered.endswith("_sp.")):
        return True
    return False


def normalize_assignment(row: dict) -> dict:
    """Stop an assignment before the first missing or placeholder rank."""
    normalized = dict(row)
    deepest_rank = "unassigned"
    deepest_name = "unassigned"

    for index, rank in enumerate(RANKS):
        value = normalized.get(rank, "").strip()
        if is_placeholder(rank, value):
            for lower_rank in RANKS[index:]:
                normalized[lower_rank] = ""
            break
        deepest_rank = rank
        deepest_name = value

    normalized["assigned_rank"] = deepest_rank
    normalized["assigned_name"] = deepest_name
    normalized["status"] = "assigned" if deepest_rank in RANKS else "unassigned"
    return normalized


def best_row(rows: list[dict]) -> dict | None:
    assigned = [
        normalize_assignment(row)
        for row in rows
        if row.get("status") == "assigned"
    ]
    assigned = [row for row in assigned if depth(row) >= 0]
    if not assigned:
        return None
    return max(
        assigned,
        key=lambda row: (
            depth(row),
            metric(row, "best_identity"),
            metric(row, "best_query_coverage"),
            metric(row, "best_bitscore"),
        ),
    )


def lineage(row: dict) -> dict[str, str]:
    maximum = depth(row)
    return {
        rank: row.get(rank, "") if index <= maximum else ""
        for index, rank in enumerate(RANKS)
    }


def lineages_conflict(first: dict, second: dict) -> bool:
    for rank in RANKS:
        left = first.get(rank, "")
        right = second.get(rank, "")
        if left and right and left != right:
            return True
    return False


def lca(rows: list[dict]) -> tuple[dict[str, str], str, str]:
    consensus = {rank: "" for rank in RANKS}
    assigned_rank = "unassigned"
    assigned_name = "unassigned"
    for rank in RANKS:
        values = [row.get(rank, "") for row in rows]
        if any(not value for value in values) or len(set(values)) != 1:
            break
        consensus[rank] = values[0]
        assigned_rank = rank
        assigned_name = values[0]
    return consensus, assigned_rank, assigned_name


def fill_unclassified_ranks(
    taxonomy: dict[str, str],
    assigned_rank: str,
    assigned_name: str,
) -> dict[str, str]:
    """Fill ranks below the deepest supported rank with an explicit label."""
    if assigned_rank not in RANKS:
        return {rank: "Unassigned" for rank in RANKS}

    filled = dict(taxonomy)
    assigned_index = RANKS.index(assigned_rank)
    label = f"Unclassified_{assigned_name}"

    for index, rank in enumerate(RANKS):
        if filled.get(rank):
            continue
        if index > assigned_index:
            filled[rank] = label
        else:
            # This should not normally occur because assignments are required
            # to be consecutive, but guarantees a fully populated TSV.
            filled[rank] = "Unassigned"
    return filled


def display_list(values: list[str]) -> str:
    """Write an explicit value instead of an empty TSV field."""
    return ",".join(values) if values else "None"


def combined_reason(
    method: str,
    conflicts: list[str],
    region_groups: dict[str, list[dict]],
) -> str:
    """Create an explicit explanation for the final classification status."""
    if conflicts:
        return "taxonomic_conflict_between_" + "_and_".join(conflicts)

    if method == "no_itsx_hits":
        return "no_itsx_hits"

    if method != "ITSx_fungal_unclassified":
        return f"classified_using_{method}"

    details = []
    for region in ["full", "ITS1", "ITS2", "5_8S"]:
        rows = region_groups.get(region, [])
        if not rows:
            continue
        reasons = sorted({
            row.get("classification_reason") or row.get("status") or "unknown_reason"
            for row in rows
        })
        details.append(f"{region}:{'|'.join(reasons)}")
    return ";".join(details) if details else "no_region_classification_available"


def choose_assignment(
    region_rows: dict[str, dict | None],
    detected_regions: list[str],
):
    if not detected_regions:
        return (
            {rank: "" for rank in RANKS},
            "unassigned",
            "unassigned",
            "no_itsx_hits",
            [],
            [],
            [],
            [],
        )

    evidence = [region for region, row in region_rows.items() if row is not None]
    informative = [region for region in INFORMATIVE if region_rows.get(region) is not None]
    conflict_regions = []

    if informative:
        selected_rows = [lineage(region_rows[region]) for region in informative]
        conflict = any(
            lineages_conflict(selected_rows[i], selected_rows[j])
            for i in range(len(selected_rows))
            for j in range(i + 1, len(selected_rows))
        )
        if conflict:
            final_lineage, final_rank, final_name = lca(selected_rows)
            conflict_regions.extend(informative)
            method = "informative_regions_LCA"
        else:
            chosen_region = max(
                informative,
                key=lambda region: (depth(region_rows[region]), -INFORMATIVE.index(region)),
            )
            chosen = region_rows[chosen_region]
            final_lineage = lineage(chosen)
            final_rank = chosen["assigned_rank"]
            final_name = chosen["assigned_name"]
            method = f"{chosen_region}_priority"

        supporting = region_rows.get("5_8S")
        if supporting is not None:
            supporting_lineage = lineage(supporting)
            if lineages_conflict(final_lineage, supporting_lineage):
                conflict_regions.extend(informative)
                conflict_regions.append("5_8S")
    else:
        supporting = region_rows.get("5_8S")
        if supporting is None:
            fungal = {rank: "" for rank in RANKS}
            fungal["kingdom"] = "Fungi"
            return (fungal, "kingdom", "Fungi", "ITSx_fungal_unclassified",
                    evidence, informative, [], detected_regions)
        chosen = supporting
        final_lineage = lineage(chosen)
        final_rank = chosen["assigned_rank"]
        final_name = chosen["assigned_name"]
        method = "5_8S_fallback"

    return (
        final_lineage,
        final_rank,
        final_name,
        method,
        evidence,
        informative,
        list(dict.fromkeys(conflict_regions)),
        detected_regions,
    )


def main():
    args = arguments()
    grouped = defaultdict(lambda: defaultdict(list))

    # Register every assembled ITS-candidate contig. Contigs absent from all
    # ITSx region files are retained with classification_reason=no_itsx_hits.
    for query_id in fasta_ids(args.assembly_fasta):
        grouped[parent_id(query_id)]

    for path in args.assignments:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                grouped[parent_id(row["query_id"])][row["region"]].append(row)

    fields = [
        "sample", "parent_contig", "status", "classification_reason", "assigned_rank", "assigned_name",
        *RANKS, "detected_regions", "evidence_regions", "assignment_method",
    ]
    records = []
    for parent in sorted(grouped):
        detected_regions = [
            region
            for region in ["full", "ITS1", "ITS2", "5_8S"]
            if grouped[parent].get(region)
        ]
        region_rows = {
            region: best_row(grouped[parent].get(region, []))
            for region in ["full", "ITS1", "ITS2", "5_8S"]
        }
        final, rank, name, method, evidence, informative, conflicts, detected = (
            choose_assignment(region_rows, detected_regions)
        )
        final = fill_unclassified_ranks(final, rank, name)
        reported_evidence = conflicts if conflicts else evidence
        reason = combined_reason(method, conflicts, grouped[parent])
        records.append({
            "sample": args.sample,
            "parent_contig": parent,
            "status": "conflicting" if conflicts else
                      "unclassified" if method == "ITSx_fungal_unclassified" else
                      "unassigned" if rank == "unassigned" else
                      "assigned",
            "classification_reason": reason,
            "assigned_rank": rank,
            "assigned_name": name,
            **final,
            "detected_regions": display_list(detected),
            "evidence_regions": display_list(reported_evidence),
            "assignment_method": method,
        })

    output_groups = [
        (args.all_output, records),
        (args.assigned_output, [row for row in records if row["status"] == "assigned"]),
        (args.unassigned_output, [row for row in records if row["status"] != "assigned"]),
    ]

    for output_path, output_records in output_groups:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            writer.writerows(output_records)


if __name__ == "__main__":
    main()
