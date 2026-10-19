rule sortmerna_staged:
    conda:
        "../envs/sortmerna.yaml"
    wildcard_constraints:
        stage="|".join(
            re.escape(stage)
            for stage in SORTMERNA_STAGES
        )
    message:
        "[SortMeRNA] Separating {wildcards.stage} reads for {wildcards.sample}"
    input:
        r1=lambda wc: sortmerna_r1(
            wc,
            results_dir=RESULTS_DIR,
            previous_stage=PREVIOUS_STAGE,
        ),
        r2=lambda wc: sortmerna_r2(
            wc,
            results_dir=RESULTS_DIR,
            previous_stage=PREVIOUS_STAGE,
        ),
        db=lambda wc: sortmerna_db(
            wc,
            config=config,
        ),
        db_idx=lambda wc: sortmerna_db_idx(
            wc,
            config=config,
        ),
    output:
        aligned_fwd=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.aligned_fwd.fq.gz",
        aligned_rev=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.aligned_rev.fq.gz",
        nonaligned_fwd=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.nonaligned_fwd.fq.gz",
        nonaligned_rev=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.nonaligned_rev.fq.gz",
    log:
        stdout=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/logs/sortmerna_{{stage}}.log"
    benchmark:
        f"{RNA_INTERMEDIATE_DIR}/{{sample}}/benchmarks/sortmerna_{{stage}}.txt"
    params:
        workdir=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}",
        aligned_prefix=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.aligned",
        nonaligned_prefix=f"{RNA_INTERMEDIATE_DIR}/{{sample}}/sortmerna/{{stage}}/{{sample}}_{{stage}}.nonaligned",
        options=config["RNA_classification"]["sortmerna_staged"]["options"]
    threads:
        config["RNA_classification"]["sortmerna_staged"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.workdir}"
        mkdir -p "$(dirname "{log.stdout}")"
 
        sortmerna \
            --ref "{input.db}" \
            --idx-dir "{input.db_idx}" \
            --workdir "{params.workdir}" \
            --threads {threads} \
            --reads "{input.r1}" \
            --reads "{input.r2}" \
            --aligned "{params.aligned_prefix}" \
            --other "{params.nonaligned_prefix}" \
            {params.options} \
            > "{log.stdout}" 2>&1

        rm -rf \
            "{params.workdir}/kvdb" \
            "{params.workdir}/readb"
        """