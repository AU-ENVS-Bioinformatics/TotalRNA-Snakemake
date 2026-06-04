#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import time

import pysam
from typing import Dict, List, Set, Tuple, Optional

from Bio import Entrez


# -------------------------
# ARGUMENTS
# -------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="NCBI organism lookup + taxonomy + CIGAR")
    parser.add_argument("--sam", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--prefix", default="SSU_")
    parser.add_argument("--taxonomy", default=None)
    parser.add_argument("--email", default="")
    return parser.parse_args()


# -------------------------
# ACCESSION EXTRACTION
# -------------------------
def extract_accession(rname: str, prefix: str) -> str:
    if rname == "*" or not rname:
        return "NA"

    if not rname.startswith(prefix):
        return "NA"

    core = rname.split(prefix, 1)[1]
    return core.split(".", 1)[0]


# -------------------------
# READ SAM/BAM
# -------------------------
def read_alignment(file_path: str, prefix: str) -> List[Tuple[str, str, str]]:
    records: List[Tuple[str, str, str]] = []

    if file_path.endswith(".bam"):
        with pysam.AlignmentFile(file_path, "rb") as handle:
            for i, read in enumerate(handle, 1):
                if read.is_unmapped:
                    continue

                rname = handle.get_reference_name(read.reference_id)
                acc = extract_accession(rname, prefix)

                if acc == "NA":
                    continue

                cigar = read.cigarstring or "*"

                records.append((read.query_name, acc, cigar))

                if i % 100000 == 0:
                    print(f"[INFO] Processed {i:,} reads")

    else:
        with open(file_path) as handle:
            for i, line in enumerate(handle, 1):
                if line.startswith("@"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue

                read_name = fields[0]
                rname = fields[2]
                cigar = fields[5]

                acc = extract_accession(rname, prefix)
                if acc == "NA":
                    continue

                records.append((read_name, acc, cigar))

                if i % 100000 == 0:
                    print(f"[INFO] Processed {i:,} lines")

    print(f"[INFO] Total records kept: {len(records):,}")
    return records


# -------------------------
# NCBI FETCH
# -------------------------
def fetch_organisms_ncbi(accessions: Set[str]) -> Dict[str, str]:
    results: Dict[str, str] = {}
    accession_list = list(accessions)

    batch_size = 50

    for i in range(0, len(accession_list), batch_size):
        batch = accession_list[i : i + batch_size]
        print(f"[INFO] Fetching batch {i} - {i+len(batch)}")

        try:
            handle = Entrez.efetch(
                db="nuccore",
                id=",".join(batch),
                rettype="gb",
                retmode="text",
            )

            from Bio import SeqIO
            records = list(SeqIO.parse(handle, "genbank"))
            handle.close()

            for rec in records:
                acc = rec.id.split(".")[0]
                organism = rec.annotations.get("organism", "Unknown")
                results[acc] = organism

        except Exception as e:
            print(f"[ERROR] Batch failed: {e}")
            for acc in batch:
                results[acc] = "Unknown"

        time.sleep(0.34)

    print(f"[INFO] Retrieved {len(results):,} accession annotations")
    return results


# -------------------------
# LOAD TAXONOMY
# -------------------------
def load_taxonomy(path: str) -> Tuple[Dict[str, List[str]], List[str]]:
    mapping: Dict[str, List[str]] = {}

    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)

        species_idx = header.index("species")

        for row in reader:
            if len(row) <= species_idx:
                continue

            mapping[row[species_idx]] = row

    return mapping, header


# -------------------------
# OUTPUT
# -------------------------
def write_tsv(
    records: List[Tuple[str, str, str]],
    accession_map: Dict[str, str],
    taxonomy_map: Optional[Dict[str, List[str]]],
    taxonomy_header: Optional[List[str]],
    output_path: str,
) -> None:

    with open(output_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")

        base_header = ["read_name", "accession", "organism_name", "cigar"]

        if taxonomy_header:
            writer.writerow(base_header + taxonomy_header)
        else:
            writer.writerow(base_header)

        for read, acc, cigar in records:
            organism = accession_map.get(acc, "Unknown")

            row = [read, acc, organism, cigar]

            if taxonomy_map and taxonomy_header:
                tax_row = taxonomy_map.get(organism)

                if tax_row:
                    row.extend(tax_row)
                else:
                    row.extend(["NA"] * len(taxonomy_header))

            writer.writerow(row)


# -------------------------
# MAIN
# -------------------------
def main() -> None:
    args = parse_args()

    if args.email:
        Entrez.email = args.email

    records = read_alignment(args.sam, args.prefix)

    accessions: Set[str] = {acc for _, acc, _ in records}
    print(f"[INFO] Unique accessions: {len(accessions):,}")

    accession_map = fetch_organisms_ncbi(accessions)

    taxonomy_map = None
    taxonomy_header = None

    if args.taxonomy:
        print(f"[INFO] Loading taxonomy: {args.taxonomy}")
        taxonomy_map, taxonomy_header = load_taxonomy(args.taxonomy)

    write_tsv(records, accession_map, taxonomy_map, taxonomy_header, args.out)

    print("[INFO] Done")


if __name__ == "__main__":
    main()