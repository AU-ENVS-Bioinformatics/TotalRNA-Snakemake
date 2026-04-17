fasta_file="trinity.Trinity"

rule transdecoder:
    input:
        fasta=f"results/trinity/{fasta_file}.fasta",
    output:
        pep=f"results/eggnog/{fasta_file}.fasta.transdecoder.pep",
    log:
        f"logs/eggnog/transdecoder_{fasta_file}.log",
    benchmark:
        f"results/benchmarks/transdecoder_{fasta_file}.txt",
    shadow:
        "minimal"
    conda:
        "../envs/transdecoder.yaml"
    shell:
        """
        TransDecoder.LongOrfs -t {input.fasta} &> {log}
        TransDecoder.Predict -t {input.fasta} &>> {log}
        
        cp {fasta_file}.fasta.transdecoder.pep {output.pep}
        """


rule annotation_eggnog:
    input:
        pep=f"results/eggnog/{fasta_file}.fasta.transdecoder.pep",
    output:
        hits=f"results/eggnog/{fasta_file}_eggnog.emapper.hits",
        annotations=f"results/eggnog/{fasta_file}_eggnog.emapper.annotations",
        seed_orthologs=f"results/eggnog/{fasta_file}_eggnog.emapper.seed_orthologs",
    log:
        f"logs/eggnog/{fasta_file}_annotation.log",
    shadow:
        "minimal"
    threads: config["threads"]["eggnog"]
    benchmark:
        f"results/benchmarks/annotation_eggnog_{fasta_file}.txt",
    params:
        EGGNOG_DIR=config.get("EGGNOG_DIR", "~/.eggnog/"),
    conda:
        "../envs/eggnogmapper.yaml"
    shell:
        """
        emapper.py -i {input.pep} -o {fasta_file}_eggnog --cpu {threads} -m \
        diamond --itype proteins --data_dir {params.EGGNOG_DIR} &> {log}

        cp {fasta_file}_eggnog.emapper.hits {output.hits}
        cp {fasta_file}_eggnog.emapper.annotations {output.annotations}
        cp {fasta_file}_eggnog.emapper.seed_orthologs {output.seed_orthologs}
        """


rule eggnog_database:
    input:
        eggnog=f"results/eggnog/{fasta_file}_eggnog.emapper.annotations",
        database="results/mRNA/database.db",
    output:
        KO_file="results/eggnog/deseq2_KO_counts.tsv",
        gene_file="results/eggnog/deseq2_gene_counts.tsv",
        eggnog_output_file="results/eggnog/eggnog_output.tsv",
    log:
        f"logs/eggnog/{fasta_file}_database.log",
    conda:
        "../envs/duckdb.yaml"
    benchmark:
        f"results/benchmarks/eggnog_database_{fasta_file}.txt",
    script:
        "../scripts/eggnog_database.py"
