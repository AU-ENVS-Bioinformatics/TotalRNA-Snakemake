rule transdecoder:
    input:
        fasta="{folder}/{file}.fasta",
    output:
        pep="{folder}/transdecoder/{file}.fasta.transdecoder.pep",
    log:
        "logs/{folder}/transdecoder/{file}.log",
    shadow:
        "minimal"
    conda:
        "../envs/transdecoder.yaml"
    shell:
        """
        cp {input.fasta} {wildcards.file}.fasta

        TransDecoder.LongOrfs -t {wildcards.file}.fasta &>> {log}

        TransDecoder.Predict -t {wildcards.file}.fasta &>> {log}

        mv {wildcards.file}.fasta.transdecoder.pep {output.pep}
        """


rule annotation_eggnog:
    input:
        pep="results/{folder}/transdecoder/{file}.fasta.transdecoder.pep",
    output:
        hits="results/{folder}/eggnog/{file}_eggnog.emapper.hits",
        annotations="results/{folder}/eggnog/{file}_eggnog.emapper.annotations",
        seed_orthologs="results/{folder}/eggnog/{file}_eggnog.emapper.seed_orthologs",
    log:
        "logs/{folder}/eggnog/{file}_annotation.log",
    shadow:
        "minimal"
    threads: config["threads"]["eggnog"]
    params:
        EGGNOG_DIR=config.get("EGGNOG_DIR", "~/.eggnog/"),
    conda:
        "../envs/eggnogmapper.yaml"
    shell:
        """
        emapper.py -i {input.pep} -o {wildcards.file}_eggnog --cpu {threads} -m \
        diamond --itype proteins --data_dir {params.EGGNOG_DIR} &> {log}

        mv {wildcards.file}_eggnog.emapper.hits {output.hits}
        mv {wildcards.file}_eggnog.emapper.annotations {output.annotations}
        mv {wildcards.file}_eggnog.emapper.seed_orthologs {output.seed_orthologs}

        """


rule eggnog_database:
    input:
        eggnog="results/{folder}/eggnog/{file}_eggnog.emapper.annotations",
        database="results/{folder}/database.db",
    output:
        touch("results/{folder}/eggnog/{file}.done"),
    log:
        "logs/{folder}/eggnog/{file}_database.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/eggnog_database.py"
