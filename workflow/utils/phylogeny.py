"""
Helpers for domain-specific SSU phylogeny.

Identifies domains with SSU candidates and filtered BLAST hits.
The Snakemake checkpoints object is passed explicitly.
"""

import csv
from typing import Any, List


def active_phylogeny_domains(*, checkpoints: Any) -> List[str]:
    """Return domains with nonempty final SSU FASTAs."""
    result = checkpoints.check_phylogeny_domains.get()

    with open(result.output.status) as handle:
        domains = [
            row["domain"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["status"] == "present"
        ]

    print(f"[Phylogeny] Domains with SSU candidates: {domains}")

    return domains


def domains_with_blast_hits(*, checkpoints: Any) -> List[str]:
    """Return active domains with filtered BLAST reference hits."""
    domains: List[str] = []

    for domain in active_phylogeny_domains(checkpoints=checkpoints):
        result = checkpoints.check_blast_top_ssu.get(domain=domain)

        with open(result.output.status) as handle:
            if handle.read().strip() == "present":
                domains.append(domain)

    print(f"[Phylogeny] Domains with filtered BLAST hits: {domains}")

    return domains