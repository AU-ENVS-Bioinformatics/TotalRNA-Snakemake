rule taxonomy_crest4:
    conda:
        "../envs/crest4.yaml"
    wildcard_constraints:
        marker=TAXONOMY_ASSEMBLY_MARKER_PATTERN
    message:
        "[CREST4] Annotating the {wildcards.marker} co-assembly"
    input:
        fasta=lambda wc: TAXONOMY_ASSEMBLY_FASTAS[wc.marker],
        matrix=rules.taxonomy_contig_count_matrix.output.matrix
    output:
        otus_by_rank=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/crest4/otus_by_rank.tsv",
    log:
        stdout=f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/logs/crest4.log"
    benchmark:
        f"{TAXONOMY_ASSEMBLY_DIR}/summary/{{marker}}/benchmarks/crest4.txt"
    threads:
        lambda wc: config["taxonomy"][wc.marker].get("crest4", {}).get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        outdir=$(dirname {output.otus_by_rank})
        mkdir -p "$outdir"

        crest4 \
            -f {input.fasta} \
            -t {threads} \
            -u {input.matrix} \
            -o "$outdir" \
            -d silvamod138pr2 \
            --search_algo blast \
            > {log.stdout} 2>&1 || true
        """