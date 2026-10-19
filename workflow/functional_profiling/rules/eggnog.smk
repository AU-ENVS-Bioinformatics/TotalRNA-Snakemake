rule eggnog_annotation:
    conda:
        "../envs/eggnog.yaml"
    message:
        "[eggNOG-mapper] Annotate predicted proteins with orthology and functional terms (GO, KEGG, COG)"
    input:
        proteins=f"{FUNCTION_DIR}/ORF_prediction/final_orf/ORF_proteins.faa"
    output:
        annotations=f"{FUNCTION_DIR}/eggnog/eggnog.emapper.annotations",
        hits=f"{FUNCTION_DIR}/eggnog/eggnog.emapper.hits",
        seed_orthologs=f"{FUNCTION_DIR}/eggnog/eggnog.emapper.seed_orthologs"
    log:
        stdout=f"{FUNCTION_DIR}/eggnog/logs/eggnog.log"
    benchmark:
        f"{FUNCTION_DIR}/eggnog/benchmarks/eggnog.txt"
    params:
        data_dir=config["functional_profiling"]["eggnog"]["database"],
        outdir=f"{FUNCTION_DIR}/eggnog",
        prefix="eggnog",
        options=config["functional_profiling"]["eggnog"].get("options", "")
    threads:
        config["functional_profiling"]["eggnog"].get("threads", 16)
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir} 

        emapper.py \
            -i {input.proteins} \
            --cpu {threads} \
            --data_dir {params.data_dir} \
            --output_dir {params.outdir} \
            --output {params.prefix} \
            {params.options} \
            > {log.stdout} 2>&1
        """