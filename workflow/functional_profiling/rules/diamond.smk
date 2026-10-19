rule diamond_blastp:
    conda:
        "../envs/diamond.yaml"
    message:
        "[DIAMOND] Search candidate TransDecoder ORFs against UniRef90"
    input:
        pep=f"{FUNCTION_DIR}/ORF_prediction/transdecoder/longorfs/longest_orfs.pep"
    output:
        diamond=f"{FUNCTION_DIR}/diamond/transdecoder/transdecoder.blastp.outfmt6"
    log:
        stdout=f"{FUNCTION_DIR}/diamond/transdecoder/logs/diamond.log"
    benchmark:
        f"{FUNCTION_DIR}/diamond/transdecoder/benchmarks/diamond.txt"
    params:
        db=config["functional_profiling"]["diamond_blastp"]["database"],
        options=config["functional_profiling"]["diamond_blastp"]["options"]
    threads:
        config["functional_profiling"]["diamond_blastp"].get("threads", 14)
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.diamond}) 

        diamond blastp \
            -q {input.pep} \
            -d {params.db} \
            -o {output.diamond} \
            --threads {threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """