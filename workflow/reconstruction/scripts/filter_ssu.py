#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from pathlib import Path


def fasta_records(path):
    """Yield contig ID, sequence length and unchanged FASTA lines."""
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
                "truncation": fields[10],
                "score": float(fields[14]),
                "evalue": float(fields[15]),
            })

    return hits


def prepare_output(path):
    """Resolve the parent without following an existing output symlink."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.parent.resolve() / path.name


def main():
    parser = argparse.ArgumentParser(
        description="Keep whole contigs with an included SSU hit covering "
                    "at least the specified fraction of the contig."
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--removed_out", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--min-ssu-coverage", type=float, default=0.80)
    args = parser.parse_args()

    if not 0 < args.min_ssu_coverage <= 1:
        parser.error("--min-ssu-coverage must be greater than 0 and at most 1.")

    source = Path(args.fasta).resolve()
    table = Path(args.tblout).resolve()

    for path in (source, table):
        if not path.is_file():
            parser.error(f"Missing input: {path}")

    outputs = [
        prepare_output(args.output),
        prepare_output(args.removed_out),
        prepare_output(args.report),
    ]
    destination, removed_out, report_out = outputs

    if len(set(outputs)) != len(outputs):
        parser.error("All output paths must be different.")

    for path in outputs:
        if path in {source, table}:
            parser.error("Output paths must differ from input paths.")

        if path.exists() and not path.is_symlink():
            if any(os.path.samefile(path, item) for item in (source, table)):
                parser.error(f"Output refers to an input file: {path}")

    hits = read_hits(table)
    retained = set()
    seen = set()
    rows = []

    for identifier, length, _ in fasta_records(source):
        if identifier in seen:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")

        seen.add(identifier)
        candidates = hits.get(identifier, [])

        if any(hit["end"] > length for hit in candidates):
            raise ValueError(f"Hit exceeds contig length: {identifier}")

        if not candidates:
            rows.append(
                f"{identifier}\t{length}\t-\t-\t-\t-\t-\t-\t-\t-\t0.000000"
                "\tremoved\tno_included_SSU_hit\n"
            )
            continue

        # First require sufficient coverage by a single included hit.
        qualifying = [
            hit for hit in candidates
            if (hit["end"] - hit["start"] + 1) / length
            >= args.min_ssu_coverage
        ]

        if qualifying:
            # Select the highest bit score; break ties by covered length.
            best = max(
                qualifying,
                key=lambda hit: (
                    hit["score"],
                    hit["end"] - hit["start"] + 1,
                ),
            )
            keep = True
            retained.add(identifier)
        else:
            # For rejected contigs, report the most covering hit.
            best = max(
                candidates,
                key=lambda hit: (
                    hit["end"] - hit["start"] + 1,
                    hit["score"],
                ),
            )
            keep = False

        coverage = (best["end"] - best["start"] + 1) / length
        action = "retained" if keep else "removed"
        reason = "passes" if keep else "low_contig_coverage"

        rows.append(
            f"{identifier}\t{length}\t{best['model']}\t{best['accession']}\t"
            f"{best['start']}\t{best['end']}\t{best['strand']}\t"
            f"{best['truncation']}\t{best['score']}\t{best['evalue']}\t"
            f"{coverage:.6f}\t{action}\t{reason}\n"
        )

    unknown = set(hits) - seen

    if unknown:
        raise ValueError(
            "Hit table contains IDs absent from FASTA: "
            + ", ".join(sorted(unknown)[:5])
        )

    # Unlink old outputs so writes cannot follow an old output symlink.
    for path in outputs:
        if os.path.lexists(path):
            path.unlink()

    with open(report_out, "w") as report:
        report.write(
            "contig\tlength\tmodel\taccession\tstart\tend\tstrand\t"
            "truncation\tbitscore\tevalue\tcontig_coverage\taction\treason\n"
        )
        report.writelines(rows)

    if len(retained) == len(seen):
        destination.symlink_to(source)
        removed_out.touch()
        print("All contigs retained; output is a symlink.")
    else:
        with open(destination, "w") as kept, open(removed_out, "w") as removed:
            for identifier, _, lines in fasta_records(source):
                if identifier in retained:
                    kept.writelines(lines)
                else:
                    removed.writelines(lines)

    print(f"Total contigs: {len(seen)}")
    print(f"Retained contigs: {len(retained)}")
    print(f"Removed contigs: {len(seen) - len(retained)}")


if __name__ == "__main__":
    main()