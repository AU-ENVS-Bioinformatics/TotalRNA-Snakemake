#!/usr/bin/env python3

import argparse
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
            "Rename SILVA FASTA headers using the lowest informative "
            "taxonomic rank and accession number."
        )
    )

    parser.add_argument(
        "--in_fasta",
        required=True,
        type=Path,
        help="Input SILVA reference FASTA file.",
    )

    parser.add_argument(
        "--out_fasta",
        required=True,
        type=Path,
        help="Output FASTA file with renamed headers.",
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

def sanitize_taxon(taxon):
    """
    Convert a taxonomic name into a FASTA-safe identifier.

    Examples
    --------

    Rule 1: Collapse 'Genus sp. ...' to genus level

        Phormidium sp. D1 2
            -> Phormidium

        Nostoc sp. Peltigera malacea cyanobiont DB3992
            -> Nostoc

        Radiococcaceae sp. enrichment culture clone GOGrpn K43
            -> Radiococcaceae

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
    # Collapse "Genus sp. anything" -> "Genus"
    #
    # Examples:
    #   Phormidium sp. D1 2 -> Phormidium
    #   Nostoc sp. cyanobiont UK53 -> Nostoc
    # -------------------------------------------------------------
    taxon = re.sub(
        r"^(.+?)\s+sp\..*$",
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

def rename_header(header, line_number):
    """
    Convert a SILVA header such as:

    >AB113665.1.1450 Bacteria;...;Nostoc;Nostoc commune

    into:

    >Nostoc_commune(AB113665.1.1450)
    """

    header = header.removeprefix(">").strip()

    try:
        accession, taxonomy = header.split(maxsplit=1)
    except ValueError as error:
        raise ValueError(
            f"Line {line_number}: header does not contain both an "
            f"accession and taxonomy: >{header}"
        ) from error

    #accession = accession.split(".", 1)[0]
    accession = accession.replace(".", "_") #accession CP016282.767613.769121 -> CP016282_767613_769121

    ranks = [
        rank.strip()
        for rank in taxonomy.split(";")
        if rank.strip()
    ]

    if not ranks:
        selected_taxon = "Unresolved"
    else:
        selected_taxon = next(
            (
                rank
                for rank in reversed(ranks)
                if is_informative(rank)
            ),
            "Unresolved",
        )

    selected_taxon = sanitize_taxon(selected_taxon)

    return f">{selected_taxon}({accession})"


def rename_fasta(input_fasta, output_fasta):
    """Rename all sequence headers while preserving the sequences."""

    if not input_fasta.is_file():
        raise FileNotFoundError(
            f"Input FASTA does not exist: {input_fasta}"
        )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    sequence_count = 0

    with input_fasta.open("r", encoding="utf-8") as infile, \
         output_fasta.open("w", encoding="utf-8") as outfile:

        for line_number, line in enumerate(infile, start=1):
            if line.startswith(">"):
                outfile.write(rename_header(line, line_number) + "\n")
                sequence_count += 1
            else:
                outfile.write(line)

    if sequence_count == 0:
        output_fasta.unlink(missing_ok=True)
        raise ValueError(
            f"No FASTA records were found in {input_fasta}"
        )

    print(f"Renamed {sequence_count} FASTA headers")
    print(f"Input:  {input_fasta}")
    print(f"Output: {output_fasta}")


def main():
    args = parse_arguments()
    rename_fasta(args.in_fasta, args.out_fasta)


if __name__ == "__main__":
    main()