HUMAN_MERGED_DIR = f"{RESULTS_DIR}/nonrRNA/humann_merged"
HUMANN_JOIN_DIR = f"{RESULTS_DIR}/nonrRNA/humann_merged/humann_link_input"

METRICS = ["genefamilies", "pathabundance", "pathcoverage"]

# --------------------------
# 1. Link input files
# --------------------------
rule humann_link_data:
    input:
        f"{RESULTS_DIR}/nonrRNA/{{sample}}/humann/{{sample}}_{{metric}}.tsv"
    output:
        f"{HUMANN_JOIN_DIR}/{{sample}}_{{metric}}.tsv"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {HUMANN_JOIN_DIR}

        ln -sf $(realpath {input}) {output}
        """

# --------------------------
# 2. Join tables
# --------------------------
rule humann_join_tables:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        lambda wc: expand(
            f"{HUMANN_JOIN_DIR}/{{sample}}_{wc.metric}.tsv",
            sample=SAMPLES,
        )
    output:
        f"{HUMAN_MERGED_DIR}/merged_{{metric}}.tsv"
    params:
        file_name=lambda wc: wc.metric
    shell:
        r"""
        set -euo pipefail
        mkdir -p {HUMAN_MERGED_DIR}

        humann_join_tables \
            --input {HUMANN_JOIN_DIR} \
            --file_name {params.file_name} \
            --output {output}
        """

# --------------------------
# 3. Split stratified
# --------------------------
rule humann_split_stratified:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        f"{HUMAN_MERGED_DIR}/merged_{{metric}}.tsv"
    output:
        stratified=f"{HUMAN_MERGED_DIR}/merged_{{metric}}_stratified.tsv",
        unstratified=f"{HUMAN_MERGED_DIR}/merged_{{metric}}_unstratified.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_split_stratified_table \
          --input {input} \
          --output {HUMAN_MERGED_DIR}
        """

# --------------------------
# 4. Renormalize (relab)
# --------------------------
rule humann_renorm:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        f"{HUMAN_MERGED_DIR}/merged_{{metric}}_unstratified.tsv"
    output:
        f"{HUMAN_MERGED_DIR}/merged_{{metric}}_relab.tsv"
    params:
        unit="relab"
    shell:
        r"""
        set -euo pipefail

        humann_renorm_table \
          --input {input} \
          --output {output} \
          --units {params.unit}
        """

# --------------------------
# 5. Regroup (ONLY genefamilies)
# --------------------------
rule humann_regroup_xrn:
    conda:
        "/data_2/rasmus/conda-envs/humann_env"
    input:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_relab.tsv"
    output:
        f"{HUMAN_MERGED_DIR}/merged_genefamilies_xrn.tsv"
    shell:
        r"""
        set -euo pipefail

        humann_regroup_table \
          --input {input} \
          --groups uniref90_rxn \
          --output {output}
        """
