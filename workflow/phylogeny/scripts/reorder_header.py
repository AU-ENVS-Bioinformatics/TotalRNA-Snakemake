#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


UNINFORMATIVE_PATTERNS = (
    re.compile(r"uncultured", re.IGNORECASE),
    re.compile(r"unclassified", re.IGNORECASE),
    re.compile(r"incertae[\s_-]*sedis", re.IGNORECASE),
    re.compile(r".*metagenome.*", re.IGNORECASE),
    re.compile(r".*Chloroplast.*", re.IGNORECASE),
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Rename FASTA headers using an optional accession-to-lineage "
            "TSV file."
        )
    )

    parser.add_argument(
        "--in_fasta",
        required=True,
        type=Path,
        help="Input reference FASTA file.",
    )

    parser.add_argument(
        "--out_fasta",
        required=True,
        type=Path,
        help="Output FASTA file with renamed headers.",
    )

    parser.add_argument(
        "--taxonomy_tsv",
        type=Path,
        help="Optional two-column TSV: accession and lineage. A header row is allowed."
    )

    return parser.parse_args()


def is_informative(taxon):
    """Return True when a taxonomic rank is informative."""

    taxon = taxon.strip()

    if not taxon:
        return False

    return not any(
        pattern.search(taxon)
        for pattern in UNINFORMATIVE_PATTERNS
    )


def load_taxonomy(taxonomy_tsv):
    """Read accession-to-lineage mappings from a two-column TSV file."""

    if not taxonomy_tsv.is_file():
        raise FileNotFoundError(
            f"Taxonomy TSV does not exist: {taxonomy_tsv}"
        )

    taxonomy = {}
    with taxonomy_tsv.open("r", encoding="utf-8", newline="") as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t")
        first_row = True

        for row in reader:
            if not row or not any(cell.strip() for cell in row):
                continue
            if row[0].lstrip().startswith("#"):
                continue

            if first_row and len(row) == 2 and (
                row[0].strip().lower() in {
                    "accession", "accession_id", "accession number"
                }
                and row[1].strip().lower() in {"lineage", "taxonomy", "taxon"}
            ):
                first_row = False
                continue

            first_row = False
            if len(row) != 2:
                raise ValueError(
                    f"Taxonomy TSV line {reader.line_num} must have exactly "
                    "two tab-separated columns."
                )

            accession, lineage = (cell.strip() for cell in row)
            if not accession or not lineage:
                raise ValueError(
                    f"Taxonomy TSV line {reader.line_num} has an empty "
                    "accession or lineage."
                )
            if accession in taxonomy:
                raise ValueError(
                    f"Duplicate accession in taxonomy TSV line "
                    f"{reader.line_num}: {accession}"
                )

            taxonomy[accession] = lineage

    if not taxonomy:
        raise ValueError(f"No accession mappings found in {taxonomy_tsv}")

    return taxonomy


def compile_accession_pattern(accessions):
    """Compile a literal, whole-accession matcher for FASTA headers."""

    alternatives = "|".join(
        re.escape(accession)
        for accession in sorted(accessions, key=len, reverse=True)
    )
    return re.compile(
        rf"(?<![A-Za-z0-9_.-])(?P<accession>{alternatives})"
        r"(?![A-Za-z0-9_.-])"
    )


def informative_taxon(lineage):
    """Select and sanitize the most specific informative lineage component."""

    for rank in reversed(lineage.split(";")):
        rank = re.sub(r"^[A-Za-z]+__", "", rank.strip())
        if is_informative(rank):
            return sanitize_taxon(rank)

    return ""


def sanitize_taxon(taxon):
    """
    Convert a taxonomic name into a FASTA-safe identifier.

    Examples
    --------

    Rule 1: Collapse any rank ending in 'sp. ...' to the name before 'sp.'

        Phormidium sp. D1 2
            -> Phormidium

        Nostoc sp. Peltigera malacea cyanobiont DB3992
            -> Nostoc

        Radiococcaceae sp. enrichment culture clone GOGrpn K43
            -> Radiococcaceae

        Parafrigoribacterium_sp.
            -> Parafrigoribacterium

    Rule 2: Keep species names but remove strain / isolate identifiers

        Phormidesmis priestleyi ULC007
            -> Phormidesmis_priestleyi

        Nostoc punctiforme ACSSI 037
            -> Nostoc_punctiforme

    Rule 3: Collapse genus + strain identifiers to genus

        Nostoc PCC-73102
            -> Nostoc

        Phormidesmis ANT.L52.6
            -> Phormidesmis

        Glaciihabitans sp.
            -> Glaciihabitans

    Rule 4: Convert remaining whitespace and special characters

        Glaciihabitans tibetensis
            -> Glaciihabitans_tibetensis

        Frigoribacterium endophyticum
            -> Frigoribacterium_endophyticum
    """

    taxon = taxon.strip()

    # -------------------------------------------------------------
    # Rule 1:
    # Collapse "Name sp. anything" -> "Name". This is applied to the
    # selected lineage rank, regardless of whether it is species or genus.
    #
    # Examples:
    #   Phormidium sp. D1 2 -> Phormidium
    #   Nostoc sp. cyanobiont UK53 -> Nostoc
    #   Parafrigoribacterium_sp. -> Parafrigoribacterium
    # -------------------------------------------------------------
    taxon = re.sub(
        r"^(.+?)(?:\s+|_+)sp\..*$",
        r"\1",
        taxon,
        flags=re.IGNORECASE,
    )

    # -------------------------------------------------------------
    # Convert whitespace before further processing.
    #
    # Example:
    #   Nostoc punctiforme ACSSI 037
    #       -> Nostoc_punctiforme_ACSSI_037
    # -------------------------------------------------------------
    taxon = re.sub(r"\s+", "_", taxon)

    # Replace problematic characters.
    taxon = re.sub(r"[^A-Za-z0-9_.-]", "_", taxon)

    # Collapse repeated underscores.
    taxon = re.sub(r"_+", "_", taxon)

    taxon = taxon.strip("_")

    # -------------------------------------------------------------
    # Rule 2:
    # If more than two underscore-separated parts exist,
    # keep only the first two.
    #
    # Examples:
    #   Phormidesmis_priestleyi_ULC007
    #       -> Phormidesmis_priestleyi
    #
    #   Nostoc_punctiforme_ACSSI_037
    #       -> Nostoc_punctiforme
    #
    #   Flexibacteraceae_bacterium_VUG-A32a
    #       -> Flexibacteraceae_bacterium
    # -------------------------------------------------------------
    parts = taxon.split("_")

    if len(parts) > 2:
        taxon = "_".join(parts[:2])

    # -------------------------------------------------------------
    # Rule 3:
    # Collapse "Genus_Strain" to "Genus" when the second
    # component looks like a culture / isolate / strain ID.
    #
    # Examples:
    #   Nostoc_PCC-73102
    #       -> Nostoc
    #
    #   Phormidesmis_ANT.L52.6
    #       -> Phormidesmis
    #
    #   Hymenobacter_JNLX01000069
    #       -> Hymenobacter
    # -------------------------------------------------------------
    parts = taxon.split("_")

    if len(parts) == 2:
        strain_pattern = re.compile(
            r"^[A-Z]{2,}[A-Z0-9.-]*\d",
            re.IGNORECASE,
        )

        if strain_pattern.match(parts[1]):
            taxon = parts[0]

    return taxon

def rename_header(header, line_number, accession_pattern=None, taxonomy=None):
    """Rename a header to its accession and, when available, taxon."""

    header_fields = header.removeprefix(">").strip().split(maxsplit=1)
    if not header_fields:
        raise ValueError(
            f"Line {line_number}: FASTA header is missing an accession: "
            f"{header.rstrip()}"
        )

    if accession_pattern is None or taxonomy is None:
        return f">{header_fields[0]}", False

    matches = {
        match.group("accession")
        for match in accession_pattern.finditer(header)
    }
    if len(matches) > 1:
        raise ValueError(
            f"Line {line_number}: FASTA header contains multiple mapped "
            f"accessions: {', '.join(sorted(matches))}"
        )

    if not matches:
        return f">{header_fields[0]}", False

    accession = matches.pop()
    taxon = informative_taxon(taxonomy[accession])
    if taxon:
        return f">{accession}_{taxon}", True

    return f">{accession}", True


def rename_fasta(input_fasta, output_fasta, taxonomy_tsv=None):
    """Rename all sequence headers while preserving the sequences."""

    if not input_fasta.is_file():
        raise FileNotFoundError(
            f"Input FASTA does not exist: {input_fasta}"
        )

    taxonomy = load_taxonomy(taxonomy_tsv) if taxonomy_tsv else None
    accession_pattern = (
        compile_accession_pattern(taxonomy)
        if taxonomy is not None
        else None
    )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    sequence_count = 0
    mapped_count = 0

    with input_fasta.open("r", encoding="utf-8") as infile, \
         output_fasta.open("w", encoding="utf-8") as outfile:

        for line_number, line in enumerate(infile, start=1):
            if line.startswith(">"):
                renamed_header, mapped = rename_header(
                    line,
                    line_number,
                    accession_pattern,
                    taxonomy,
                )
                outfile.write(renamed_header + "\n")
                sequence_count += 1
                mapped_count += mapped
            else:
                outfile.write(line)

    if sequence_count == 0:
        output_fasta.unlink(missing_ok=True)
        raise ValueError(
            f"No FASTA records were found in {input_fasta}"
        )

    print(f"Renamed {sequence_count} FASTA headers")
    if taxonomy_tsv:
        print(f"Matched taxonomy for {mapped_count} FASTA headers")
    print(f"Input:  {input_fasta}")
    print(f"Output: {output_fasta}")


def main():
    args = parse_arguments()
    rename_fasta(args.in_fasta, args.out_fasta, args.taxonomy_tsv)


if __name__ == "__main__":
    main()