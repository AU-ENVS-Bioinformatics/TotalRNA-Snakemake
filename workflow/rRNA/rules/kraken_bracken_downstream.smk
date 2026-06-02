rule kraken_biom:
    conda:
        "../envs/kraken_biom.yaml"
    message:
        "[Kraken-Biom] converts kraken report to BIOM format for for taxonomic classification"
    input:
        expand(
            f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.k2report",sample=SAMPLES
        )
    output:
        biom=f"{RESULTS_DIR}/rRNA/taxonomy/bracken_species.biom"
    log:
        f"{RESULTS_DIR}/rRNA/taxonomy/kraken_biom.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/taxonomy/bracken.txt"
    params:
        options=config["rRNA"]["kraken_biom"]["options"]
    shell:
        r"""
        mkdir -p $(dirname {output.biom})

        kraken-biom {input} \
            -o {output.biom} \
            {params.options} \
            > {log} 2>&1
        """

#kraken-biom *.k2report -o bracken_genus.biom --min S --max D