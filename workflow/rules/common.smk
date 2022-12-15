def choose_filepath(configuration: dict):
    """
    It returns a non-filtered mRNA fasta or a abundance_filtered fasta based
    on configuration file.
    """
    is_filter = configuration.get("filter_by_abundance_before_align", False)
    filtered = "results/mRNA/trinity/contigs_ncrna_filtered_AbundanceFiltered.fasta"
    non_filtered = "results/mRNA/trinity/contigs_ncrna_filtered.fasta"
    return filtered if is_filter else non_filtered
