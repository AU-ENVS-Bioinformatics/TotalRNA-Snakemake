rule bracken:
    conda:
        "../envs/bracken.yaml"
    message:
        "[Bracken] resetimates abundance of kraken results for {wildcards.sample} for taxonomic classification"
    input:
        kraken_report=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.k2report",
    output:
        bracken_output=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.bracken.tsv",
        bracken_kreport_output=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{sample}}.bracken.k2report"
    log:
        stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/bracken.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/bracken.txt"
    params:
        db=config["databases"]["kraken_rRNA_db"],
        options=config["rRNA"]["bracken"]["options"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bracken_output})

        bracken \
          -d {params.db} \
          -i {input.kraken_report} \
          -o {output.bracken_output} \
          -w {output.bracken_kreport_output} \
          {params.options} \
          > {log.stdout} 2>&1
        """
#bracken -d /data_2/Databases/silva_kraken_db/SILVA_138_2_k2db -i ANN11.report.txt -o test.genus.tsv -w test.k2report -r 150 -l G -t 10
#kraken-biom *.k2report -o bracken_genus.biom --min S --max D