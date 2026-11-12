rule taxonomy_bwamem2_index:
    conda:
        "../envs/bwamem2_samtools.yaml"
    wildcard_constraints:
        marker=TAXONOMY_ASSEMBLY_MARKER_PATTERN
    message:
        "[BWA-MEM2] Indexing the {wildcards.marker} co-assembly"
    input:
        fasta=lambda wc: TAXONOMY_ASSEMBLY_FASTAS[wc.marker]
    output:
        index=multiext(
            TAXONOMY_ASSEMBLY_INDEX_PREFIX,
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac",
        )
    log:
        stdout=f"{TAXONOMY_ASSEMBLY_DIR}/{{marker}}/logs/bwamem2_index.log"
    params:
        prefix=TAXONOMY_ASSEMBLY_INDEX_PREFIX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{params.prefix}")" "$(dirname "{log.stdout}")"
        bwa-mem2 index -p "{params.prefix}" "{input.fasta}" > "{log.stdout}" 2>&1
        """

rule taxonomy_bwamem2_mem:
    conda:
        "../envs/bwamem2_samtools.yaml"
    wildcard_constraints:
        marker=TAXONOMY_ASSEMBLY_MARKER_PATTERN
    message:
        "[BWA-MEM2] Aligning {wildcards.sample} {wildcards.marker} reads to the co-assembly"
    input:
        index=rules.taxonomy_bwamem2_index.output.index,
        read_1=lambda wc: assembly_taxonomy_read(wc, 1, results_dir=RESULTS_DIR),
        read_2=lambda wc: assembly_taxonomy_read(wc, 2, results_dir=RESULTS_DIR)
    output:
        bam=f"{TAXONOMY_ASSEMBLY_DIR}/{{sample}}/{{marker}}/alignment/{{sample}}.bam"
    log:
        stdout=f"{TAXONOMY_ASSEMBLY_DIR}/{{sample}}/{{marker}}/logs/{{sample}}.bwamem2.log"
    benchmark:
        f"{TAXONOMY_ASSEMBLY_DIR}/{{sample}}/{{marker}}/benchmarks/{{sample}}.bwamem2.txt"
    params:
        prefix=TAXONOMY_ASSEMBLY_INDEX_PREFIX
    threads:
        lambda wc: bwamem_threads(wc, config=config)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.bam}")" "$(dirname "{log.stdout}")"
        bwa_threads=$(( {threads} > 1 ? {threads} - 1 : 1 ))
        (
            bwa-mem2 mem -t "$bwa_threads" "{params.prefix}" \
                "{input.read_1}" \
                "{input.read_2}" |
                samtools sort -o "{output.bam}" -
        ) > "{log.stdout}" 2>&1
        """