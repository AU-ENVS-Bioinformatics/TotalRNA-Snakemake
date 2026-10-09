#!/usr/bin/env python3
"""Build a concise viral candidate summary from geNomad and CheckV.

No directories are created. Individual gene details remain in geNomad output.
The existing eggNOG subset and Salmon files are produced by separate rules.
"""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

NA = "NA"
CANDIDATE_COLUMNS = [
    "contig_id", "length", "virus_score", "n_genes", "n_hallmarks",
    "marker_enrichment", "genomad_taxonomy", "checkv_quality",
    "completeness", "contamination", "warnings", "hallmark_descriptions",
]


def read_table(path, required=()):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if not set(required).issubset(columns):
            raise ValueError(f"{path}: missing columns {sorted(set(required) - set(columns))}")
        if len(columns) != len(set(columns)):
            raise ValueError(f"{path}: duplicate column names")
        rows = []
        for number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}: wrong number of fields at line {number}")
            rows.append(row)
    return columns, rows


def write_table(path, columns, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", restval=NA)
        writer.writeheader()
        writer.writerows(rows)


def parent_id(candidate_id):
    return re.sub(r"\|provirus_\d+_\d+$", "", candidate_id)


def build_annotations(args):
    candidate_columns, candidates = read_table(args.candidates,
                                               ("seq_name", "parent_contig_id"))
    gene_columns, genes = read_table(args.genomad_genes,
                                    ("gene", "start", "end", "strand",
                                     "virus_hallmark", "annotation_description"))
    by_candidate = {}
    by_parent = defaultdict(list)
    for row in candidates:
        cid = row["seq_name"]
        if cid in by_candidate:
            raise ValueError(f"Duplicate candidate: {cid}")
        by_candidate[cid] = row
        by_parent[row["parent_contig_id"]].append(cid)

    # geNomad gene IDs end with an underscore and gene number. For excised
    # proviruses some outputs retain original parent-based gene identifiers.
    # Assign an exact candidate only when its ID is explicitly in the gene ID,
    # or the sole candidate is the whole parent contig. Otherwise keep a
    # parent-only association; do not guess provirus membership/coordinates.
    genes_by_candidate = defaultdict(list)
    seen_genes = set()
    unresolved = 0
    for gene in genes:
        gid = gene["gene"]
        if gid in seen_genes:
            raise ValueError(f"Duplicate geNomad gene: {gid}")
        seen_genes.add(gid)
        match = re.fullmatch(r"(.+)_\d+", gid)
        if not match:
            raise ValueError(f"Cannot extract contig from geNomad gene ID: {gid}")
        prefix = match.group(1)
        parent = parent_id(prefix)
        if parent not in by_parent:
            raise ValueError(f"geNomad gene has no candidate parent: {gid}")
        exact = prefix if prefix in by_candidate else NA
        if exact != NA:
            genes_by_candidate[exact].append(gene)
        else:
            unresolved += 1
    candidate_summary = []
    hallmark_mismatches = 0
    for row in candidates:
        cid = row["seq_name"]
        hallmarks = [gene for gene in genes_by_candidate[cid]
                     if gene["virus_hallmark"] == "1"]
        summary = {key: row.get(key, NA) for key in CANDIDATE_COLUMNS}
        summary["contig_id"] = cid
        summary["hallmark_descriptions"] = "|".join(
            gene["annotation_description"] or NA for gene in hallmarks
        ) or NA
        candidate_summary.append(summary)
        if row.get("n_hallmarks", NA) != str(len(hallmarks)):
            hallmark_mismatches += 1

    write_table(args.candidate_output, CANDIDATE_COLUMNS, candidate_summary)
    with Path(args.log).open("w") as handle:
        handle.write(f"Viral candidates: {len(candidates)}\n")
        handle.write(f"geNomad gene rows read: {len(genes)}\n")
        handle.write(f"geNomad genes with parent-only association: {unresolved}\n")
        handle.write(f"Candidates with hallmark count differing from assigned gene records: {hallmark_mismatches}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    annotations = commands.add_parser("annotations")
    for name in ("candidates", "genomad-genes", "candidate-output", "log"):
        annotations.add_argument(f"--{name}", required=True)
    annotations.set_defaults(run=build_annotations)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
