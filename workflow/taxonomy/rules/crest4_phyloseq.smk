rule taxonomy_crest4_phyloseq:
    conda:
        "../envs/r_phyloseq.yaml"
    wildcard_constraints:
        marker=TAXONOMY_ASSEMBLY_MARKER_PATTERN
    message:
        "[phyloseq] Converting CREST4 assignments for {wildcards.marker}"
    input:
        assignments=rules.taxonomy_crest4.output.assignments,
        counts=rules.taxonomy_contig_count_matrix.output.matrix
    output:
        otu=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/crest4/mapped_reads_to_contigs.tsv",
        gg=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/crest4/mapped_reads_to_contigs_gg.tsv",
        rds=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/crest4/physeq.Rds"
    log:
        stdout=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/logs/crest4_phyloseq.log"
    benchmark:
        f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/benchmarks/crest4_phyloseq.txt"
    params:
        script=workflow.source_path("../scripts/crest4_to_phyloseq.R")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.rds}")" "$(dirname "{log.stdout}")"
        Rscript "{params.script}" \
            "{input.assignments}" \
            "{input.counts}" \
            "{output.otu}" \
            "{output.gg}" \
            "{output.rds}" \
            > "{log.stdout}" 2>&1
        """