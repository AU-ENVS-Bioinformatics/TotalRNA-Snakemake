rule bracken:
    conda:
        "../envs/bracken.yaml"
    message:
        "[Bracken] re-estimate abundance for {wildcards.sample} for {wildcards.database}"
    input:
        kraken_report=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{database}}/{{sample}}.{{database}}.k2report",
    output:
        bracken_output=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{database}}/{{sample}}.{{database}}.bracken.tsv",
        bracken_kreport_output=f"{RESULTS_DIR}/rRNA/{{sample}}/classification/{{database}}/{{sample}}.{{database}}.bracken.k2report"
    log:
        stdout=f"{RESULTS_DIR}/rRNA/{{sample}}/logs/bracken_{{database}}.log"
    benchmark:
        f"{RESULTS_DIR}/rRNA/{{sample}}/benchmarks/bracken_{{database}}.txt"
    params:
        db=kraken_db,
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