rule bracken:
    conda:
        "../envs/bracken.yaml"
    message:
        "[Bracken] re-estimate abundance for {wildcards.sample} for {wildcards.database}"
    input:
        kraken_report=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/{{sample}}.{{database}}.k2report"
    output:
        bracken_output=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/{{sample}}.{{database}}.bracken.tsv",
        bracken_kreport_output=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/{{sample}}.{{database}}.bracken.k2report"
    log:
        stdout=f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/logs/{{sample}}.{{database}}.bracken.log"
    benchmark:
        f"{TAXONOMY_DIR}/{{sample}}/SSU/{{database}}/benchmarks/{{sample}}.{{database}}.bracken.txt"
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