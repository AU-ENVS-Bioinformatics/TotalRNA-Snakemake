#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import time
from typing import Dict, List, Set, Tuple, Optional

from Bio import Entrez
from Bio import SeqIO


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# -------------------------
# ALIGNMENT UTILITIES
# -------------------------
def cigar_aligned_bases(cigar: str) -> int:
    return sum(
        int(length)
        for length, op in CIGAR_RE.findall(cigar)
        if op in {"M", "=", "X"}
    )


def extract_nm(nm_field: str) -> int:
    if nm_field.startswith("NM:i:"):
        return int(nm_field.split(":")[-1])
    return 0


def extract_accession(rname: str, prefix: str) -> str:
    if rname == "*" or not rname.startswith(prefix):
        return "NA"
    core = rname.split(prefix, 1)[1]
    return core.split(".", 1)[0]


def passes_filter(fraction: float, nm: int, min_fraction: float, max_nm: Optional[int]) -> bool:
    if fraction < min_fraction:
        return False
    if max_nm is not None and nm > max_nm:
        return False
    return True


# -------------------------
# TAXONOMY
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
# NCBI
# -------------------------
def fetch_organisms_ncbi(accessions: Set[str]) -> Dict[str, str]:
    results: Dict[str, str] = {}
    accession_list = list(accessions)

    for i in range(0, len(accession_list), 50):
        batch = accession_list[i:i+50]
        try:
            handle = Entrez.efetch(
                db="nuccore",
                id=",".join(batch),
                rettype="gb",
                retmode="text"
            )
            records = list(SeqIO.parse(handle, "genbank"))
            handle.close()

            for rec in records:
                acc = rec.id.split(".")[0]
                results[acc] = rec.annotations.get("organism", "Unknown")

        except Exception:
            for acc in batch:
                results[acc] = "Unknown"

        time.sleep(0.34)

    return results


# -------------------------
# PROCESS READS
# -------------------------
def process_reads(
    tsv_path: str,
    prefix: str,
    min_fraction: float,
    max_nm: Optional[int],
) -> Tuple[List[Tuple[str, str, str, str]], Set[str], Set[str], Set[str]]:

    read_map: Dict[str, Dict[str, object]] = {}
    records: List[Tuple[str, str, str, str]] = []

    with open(tsv_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue

            read_id, rname, cigar, _, nm_field, read_len = fields
            nm = extract_nm(nm_field)
            read_length = int(read_len)

            if read_id not in read_map:
                read_map[read_id] = {"all_ssu": True, "any_pass": False, "n_alignments": 0}

            read_map[read_id]["n_alignments"] += 1

            if not rname.startswith(prefix):
                read_map[read_id]["all_ssu"] = False

            aligned = cigar_aligned_bases(cigar)
            fraction = aligned / read_length

            if passes_filter(fraction, nm, min_fraction, max_nm):
                read_map[read_id]["any_pass"] = True

            acc = extract_accession(rname, prefix)
            if acc != "NA":
                records.append((read_id, rname, acc, cigar))

    passing, failed, discarded = set(), set(), set()

    for rid, flags in read_map.items():
        if not flags["all_ssu"] or flags["n_alignments"] < 2:
            discarded.add(rid)
        elif flags["any_pass"]:
            passing.add(rid)
        else:
            failed.add(rid)

    filtered_records = [r for r in records if r[0] in passing]
    return filtered_records, passing, failed, discarded


# -------------------------
# STATISTICS
# -------------------------
def write_statistics(
    path: str,
    passing: Set[str],
    failed: Set[str],
    discarded: Set[str],
    taxonomy_discarded: Set[str],
    failed_taxonomy: Set[str],
    args
):
    total = len(passing | failed | discarded)

    def pct(n): return (n / total * 100) if total else 0

    with open(path, "w") as f:
        f.write("n_reads\tpercentage\tfilename\tdescription\n")

        rows = [
            (passing, args.read_id, "kept read ids because of passing alignment filtering"),
            (taxonomy_discarded, args.taxonomy_discarded_id, "removed passing read ids based on taxonomical filtering"),
            (failed - failed_taxonomy, args.failed_id, "SSU reads that failed alignment filtering"),
            (failed_taxonomy, args.failed_taxonomy, "failed reads removed due to taxonomy filtering"),
            (discarded, args.discarded_id, "reads discarded because of prefix mismatch or incomplete pairing"),
        ]

        for s, fname, desc in rows:
            f.write(f"{len(s)}\t{pct(len(s)):.2f}%\t{fname}\t{desc}\n")

        f.write(f"{total}\t100.00%\tNA\tTOTAL\n")


# -------------------------
# OUTPUT
# -------------------------
def write_output(
    records,
    accession_map,
    taxonomy_map,
    taxonomy_header,
    discard_set,
    keep_set,
    output_path,
    read_id_path,
    taxonomy_discarded_path,
    failed_reads,
    failed_taxonomy_path,
):

    taxonomy_discarded, failed_taxonomy = set(), set()

    with open(output_path, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["read_name", "reference", "accession", "organism_name", "cigar"] + (taxonomy_header or []))

        seen = set()
        rid_out = open(read_id_path, "w") if read_id_path else None

        for read, rname, acc, cigar in records:
            if read in seen:
                continue
            seen.add(read)

            organism = accession_map.get(acc, "Unknown")
            tax_row = taxonomy_map.get(organism) if taxonomy_map else None

            reject = False
            if taxonomy_map and tax_row:
                tax_lower = [x.lower() for x in tax_row]

                if keep_set and not any(v in keep_set for v in tax_lower):
                    reject = True
                if discard_set and any(v in discard_set for v in tax_lower):
                    reject = True
            elif taxonomy_map:
                reject = True

            if reject:
                taxonomy_discarded.add(read)
                if read in failed_reads:
                    failed_taxonomy.add(read)
                continue

            row = [read, rname, acc, organism, cigar]
            if tax_row:
                row.extend(tax_row)

            writer.writerow(row)
            if rid_out:
                rid_out.write(read + "\n")

        if rid_out:
            rid_out.close()

    # write files
    if taxonomy_discarded_path:
        with open(taxonomy_discarded_path, "w") as f:
            f.write("\n".join(sorted(taxonomy_discarded)))

    if failed_taxonomy_path:
        with open(failed_taxonomy_path, "w") as f:
            f.write("\n".join(sorted(failed_taxonomy)))

    return taxonomy_discarded, failed_taxonomy


# -------------------------
# CLI
# -------------------------
def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--id_info", required=True)
    p.add_argument("--fraction", type=float, required=True)
    p.add_argument("--matches", type=int, default=None)
    p.add_argument("--prefix", required=True)

    p.add_argument("--out", required=True)

    p.add_argument("--read_id")
    p.add_argument("--failed_id")
    p.add_argument("--failed_taxonomy")
    p.add_argument("--discarded_id")
    p.add_argument("--taxonomy_discarded_id")

    p.add_argument("--statistics")

    p.add_argument("--email", required=True)
    p.add_argument("--taxonomy")

    p.add_argument("--discard")
    p.add_argument("--keep")

    return p.parse_args()


# -------------------------
# MAIN
# -------------------------
def main():
    args = parse_args()
    Entrez.email = args.email

    records, passing, failed, discarded = process_reads(
        args.id_info, args.prefix, args.fraction, args.matches
    )

    acc_map = fetch_organisms_ncbi({a for _, _, a, _ in records})

    tax_map, tax_header = (None, None)
    if args.taxonomy:
        tax_map, tax_header = load_taxonomy(args.taxonomy)

    discard_set = set(args.discard.lower().split(",")) if args.discard else set()
    keep_set = set(args.keep.lower().split(",")) if args.keep else set()

    taxonomy_discarded, failed_taxonomy = write_output(
        records,
        acc_map,
        tax_map,
        tax_header,
        discard_set,
        keep_set,
        args.out,
        args.read_id,
        args.taxonomy_discarded_id,
        failed,
        args.failed_taxonomy,
    )

    if args.failed_id:
        with open(args.failed_id, "w") as f:
            f.write("\n".join(sorted(failed)))

    if args.discarded_id:
        with open(args.discarded_id, "w") as f:
            f.write("\n".join(sorted(discarded)))

    if args.statistics:
        write_statistics(
            args.statistics,
            passing,
            failed,
            discarded,
            taxonomy_discarded,
            failed_taxonomy,
            args
        )


if __name__ == "__main__":
    main()