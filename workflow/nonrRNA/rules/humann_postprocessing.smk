HUMANN_DIR = f"{RESULTS_DIR}/nonrRNA/humann"
HUMAN_MERGED_DIR = f"{RESULTS_DIR}/nonrRNA/humann_merged"

# --------------------------
# Metrics
# --------------------------
METRICS = ["genefamilies", "pathabundance", "pathcoverage"]

# --------------------------
# 1. Join tables
# --------------------------
rule humann_join_tables:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        lambda wc: expand(
            f"{HUMANN_DIR}/{{sample}}_{wc.metric}.tsv",
            sample=SAMPLES,
        )
    output:
        f"{HUMAN_MERGED_DIR}/merged_{{metric}}.tsv"
    params:
        metric="{metric}"
    shell:  
        r"""
        set -euo pipefail
        mkdir -p {HUMAN_MERGED_DIR}

        humann_join_tables \
          --input {HUMANN_DIR} \
          --output {output} \
          --file_name {params.metric}
        """

# --------------------------
# 2. Normalize (all metrics)
# --------------------------
rule humann_renorm:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        humann_merged = expand(f"{HUMAN_MERGED_DIR}/merged_{{metric}}.tsv",metric=METRICS)
    output:
        humann_renorm = expand(f"{HUMAN_MERGED_DIR}/merged_renorm_{{metric}}.tsv",metric=METRICS)
    params:
        unit="cpm"
    shell:
        r"""
        set -euo pipefail

        humann_renorm_table \
          --input {input.humann_merged} \
          --output {output.humann_renorm} \
          --units {params.unit}
        """

# --------------------------
# 3. Regroup (ONLY genefamilies)
# --------------------------
rule humann_regroup_ko:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_cpm.tsv"
    output:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_ko.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_regroup_table \
          --input {input} \
          --output {output} \
          --groups uniref90_ko
        """

# Optional EggNOG regrouping
rule humann_regroup_eggnog:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_cpm.tsv"
    output:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_eggnog.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_regroup_table \
          --input {input} \
          --output {output} \
          --groups uniref90_eggnog
        """

# --------------------------
# 4. Split stratified (KO)
# --------------------------
rule humann_split_stratified_ko:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_ko.tsv"
    output:
        stratified=f"{HUMAN_MERGED_DIR}/merged_genefamilies_ko_stratified.tsv",
        unstratified=f"{HUMAN_MERGED_DIR}/merged_genefamilies_ko_unstratified.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_split_stratified_table \
          --input {input} \
          --output {HUMAN_MERGED_DIR}
        """