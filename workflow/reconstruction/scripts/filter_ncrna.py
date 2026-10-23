#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from pathlib import Path


def fasta_records(path):
    """Yield contig ID, sequence length and original FASTA lines."""
    with open(path) as handle:
        identifier = None
        lines = []
        length = 0

        for line in handle:
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, length, lines

                identifier = line[1:].split()[0]
                lines = [line]
                length = 0
            elif identifier is not None:
                lines.append(line)
                length += len("".join(line.split()))

        if identifier is not None:
            yield identifier, length, lines


def read_hits(path):
    """Read included hits from standard cmsearch --tblout output."""
    hits = defaultdict(list)

    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 18:
                raise ValueError(f"Malformed cmsearch row at line {number}.")

            if fields[16] != "!":
                continue

            start, end = sorted((int(fields[7]), int(fields[8])))

            if start < 1:
                raise ValueError(f"Invalid coordinates at line {number}.")

            hits[fields[0]].append({
                "model": fields[2],
                "accession": fields[3],
                "start": start,
                "end": end,
                "strand": fields[9],
                "score": fields[14],
                "evalue": fields[15],
            })

    return hits


def covered_bases(hits, length, identifier):
    """Count the union of hit intervals without counting overlaps twice."""
    regions = sorted((hit["start"], hit["end"]) for hit in hits)

    if not regions:
        return 0

    if any(end > length for _, end in regions):
        raise ValueError(f"Hit coordinates exceed contig length: {identifier}")

    start, end = regions[0]
    covered = 0

    for next_start, next_end in regions[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            covered += end - start + 1
            start, end = next_start, next_end

    return covered + end - start + 1


def output_path(path):
    """Resolve the parent without following an existing output symlink."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.parent.resolve() / path.name


def main():
    parser = argparse.ArgumentParser(
        description="Remove whole contigs exceeding the ncRNA coverage threshold."
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--removed_out", required=True)
    parser.add_argument("--removed_report", required=True)
    parser.add_argument("--max-ncrna-coverage", type=float, default=0.90)
    args = parser.parse_args()

    if not 0 <= args.max_ncrna_coverage <= 1:
        parser.error("--max-ncrna-coverage must be between 0 and 1.")

    source = Path(args.fasta).resolve()
    table = Path(args.tblout).resolve()

    for path in (source, table):
        if not path.is_file():
            parser.error(f"Missing input: {path}")

    destination = output_path(args.output)
    removed_out = output_path(args.removed_out)
    removed_report = output_path(args.removed_report)
    outputs = (destination, removed_out, removed_report)

    if len(set(outputs)) != len(outputs):
        parser.error("All output paths must be different.")

    for path in outputs:
        if path in {source, table}:
            parser.error("Output paths must differ from input paths.")

        if path.exists() and not path.is_symlink():
            if any(os.path.samefile(path, item) for item in (source, table)):
                parser.error(f"Output refers to an input file: {path}")

    hits = read_hits(table)
    excluded = {}
    seen = set()

    print("contig\tlength\tncrna_covered_bases\tncrna_fraction\taction")

    for identifier, length, _ in fasta_records(source):
        if identifier in seen:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")

        seen.add(identifier)

        covered = covered_bases(hits.get(identifier, []), length, identifier)
        fraction = covered / length if length else 0
        remove = fraction > args.max_ncrna_coverage

        if remove:
            excluded[identifier] = (length, covered, fraction)

        action = "removed" if remove else "retained"
        print(f"{identifier}\t{length}\t{covered}\t{fraction:.6f}\t{action}")

    unknown = set(hits) - seen

    if unknown:
        raise ValueError(
            "Hit table contains IDs absent from the FASTA: "
            + ", ".join(sorted(unknown)[:5])
        )

    # Remove old outputs without writing through existing symlinks.
    for path in outputs:
        if os.path.lexists(path):
            path.unlink()

    if not excluded:
        destination.symlink_to(source)
        removed_out.touch()
        print("\nNo contigs removed; retained FASTA is a symlink.")
    else:
        with open(destination, "w") as retained, open(removed_out, "w") as removed:
            for identifier, _, lines in fasta_records(source):
                if identifier in excluded:
                    removed.writelines(lines)
                else:
                    retained.writelines(lines)

        print("\nRetained and removed FASTAs written.")

    # One row per included Rfam hit on each removed contig.
    with open(removed_report, "w") as report:
        report.write(
            "contig\tlength\tncrna_covered_bases\tncrna_fraction\t"
            "rfam_name\trfam_accession\tstart\tend\tstrand\tbitscore\tevalue\n"
        )

        for identifier, (length, covered, fraction) in excluded.items():
            for hit in hits[identifier]:
                report.write(
                    f"{identifier}\t{length}\t{covered}\t{fraction:.6f}\t"
                    f"{hit['model']}\t{hit['accession']}\t"
                    f"{hit['start']}\t{hit['end']}\t{hit['strand']}\t"
                    f"{hit['score']}\t{hit['evalue']}\n"
                )

    print(f"Total contigs: {len(seen)}")
    print(f"Removed contigs: {len(excluded)}")
    print(f"Retained contigs: {len(seen) - len(excluded)}")


if __name__ == "__main__":
    main()