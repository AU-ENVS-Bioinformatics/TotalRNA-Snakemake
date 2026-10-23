################################################################################
# SEARCH SSU ASSEMBLIES WITH EACH CONFIGURED COVARIANCE MODEL
################################################################################

CMSEARCH_MODULES = list(config["reconstruction"]["cmsearch"]["models"])

rule cmsearch_ssu:
    conda:
        "../envs/infernal.yaml"
    wildcard_constraints:
        module_name="|".join(CMSEARCH_MODULES)
    message:
        "[cmsearch] Searching SSU assemblies using the {wildcards.module_name} model"
    input:
        fasta=f"{ASSEMBLY_DIR}/concatenated/SSU/cross_sample_SSU.fasta",
        model=lambda wc: config["reconstruction"]["cmsearch"]["models"][wc.module_name]
    output:
        tbl=f"{ASSEMBLY_DIR}/SSU_candidates/cmsearch/{{module_name}}/cross_sample_SSU.cmsearch.tblout",
        report=f"{ASSEMBLY_DIR}/SSU_candidates/cmsearch/{{module_name}}/cross_sample_SSU.cmsearch.txt"
    log:
        stdout=f"{ASSEMBLY_DIR}/SSU_candidates/cmsearch/{{module_name}}/logs/cmsearch.log",
        time=f"{ASSEMBLY_DIR}/SSU_candidates/cmsearch/{{module_name}}/benchmarks/cmsearch.time.txt"
    benchmark:
        f"{ASSEMBLY_DIR}/SSU_candidates/cmsearch/{{module_name}}/benchmarks/cmsearch.benchmark.tsv"
    params:
        options=config["reconstruction"]["cmsearch"].get("options", "")
    threads:
        config["reconstruction"]["cmsearch"].get("threads", 12)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.tbl:q})

        /usr/bin/time \
            -v \
            -o {log.time:q} \
            cmsearch \
                --cpu {threads} \
                {params.options} \
                --tblout {output.tbl:q} \
                -o {output.report:q} \
                {input.model:q} \
                {input.fasta} \
                2> {log.stdout}
        """