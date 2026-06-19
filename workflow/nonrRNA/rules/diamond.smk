rule diamond:
    conda:
        "../envs/diamond.yaml"
    message:
        "[DIAMOND] sequence aligner for protein and translated DNA searches for {wildcards.sample}"
    input:
        nonrRNA_concatenate=rules.concatenate.output.nonrRNA_concatenate,
    output:
        diamond_uniref90=f"{RESULTS_DIR}/nonrRNA/{{sample}}/diamond/{{sample}}_uniref90.tsv",
    log:
        stdout=f"{RESULTS_DIR}/nonrRNA/{{sample}}/logs/diamond.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/{{sample}}/benchmarks/diamond.txt"
    params:
        protein_db=config["databases"]["diamond_proteindb"],
        options=config["nonrRNA"]["diamond"]["options"]
    threads:
        config["nonrRNA"]["diamond"].get("threads", 14)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.genefamilies})

        diamond blastx 
            -q {input.nonrRNA_concatenate} \
            -d {params.protein_db}
            -o {params.diamond_uniref90} \
            --threads {threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """

