rule eggnog:
    conda:
        "../envs/eggnog.yaml"
    message:
        "[Eggnog] Annotation of GO, KEEG, COG, genes of the predicted ORF's"
    input:
        pep = f"{RESULTS_DIR}/nonrRNA/predicted/transcripts.fasta.transdecoder.pep"
    output:
        annotations=f"{RESULTS_DIR}/nonrRNA/eggnog/eggnog.emapper.annotations",
        hits=f"{RESULTS_DIR}/nonrRNA/eggnog/eggnog.emapper.hits",
        seed_orthologs=f"{RESULTS_DIR}/nonrRNA/eggnog/eggnog.emapper.seed_orthologs",
    log:
        stdout = f"{RESULTS_DIR}/nonrRNA/eggnog/logs/eggnog.log"
    benchmark:
        f"{RESULTS_DIR}/nonrRNA/eggnog/benchmarks/eggnog.txt"
    params:
        diamond_db=config["databases"]["emapper_diamond_proteindb"],
        options=config["nonrRNA"]["eggnog"]["options"],
        outdir=f"{RESULTS_DIR}/nonrRNA/eggnog/",
        prefix="eggnog"
    threads:
        config["nonrRNA"]["eggnog"].get("threads", 16)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.annotations})

        emapper.py \
            -i {input.pep} \
            --cpu {threads} \
            --data_dir {params.diamond_db} \
            --output_dir {params.outdir} \
            -o {params.prefix} \
            {params.options} \
            > {log.stdout} 2>&1
        """

#functional annotation & orthology assignment
# GO terms
# KEGG
# COG categories
# gene names

#emapper.py -i transcripts.fasta.transdecoder.pep 
#--itype proteins 
#-m diamond 
#--data_dir /data_2/Databases/eggnog_mapper 
#--cpu 32 
#--output_dir ../eggnog/ 
#--sensmode sensitive 
#-o eggnog
