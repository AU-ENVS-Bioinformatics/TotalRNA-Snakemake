rule reformat_rRNA_for_bbduk:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBDuk] Reformatting ribodetected seperated reads for {wildcards.sample} to ensure read order for downstream processing"
    input:
        r1=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r1.fastq.gz",
        r2=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r2.fastq.gz",
        ref=config["databases"]["sortmeRNA_ssu"]
    output:
        interleaved_reads=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}.interleaved.rRNA.fastq.gz",
    log:
        stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbduk_ssu.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbduk_ssu.txt"
    params:
        options="-Xmx40g"     
    threads:
        config["RNA"]["bbduk"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.interleaved_reads})

        reformat.sh in={input.r1} in2={input.r2} interleaved=t threads={threads} out={output.interleaved_reads} {params.options} > {log.stdout} 2>&1
        """

rule bbduk_ssu:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBDuk] Extracting SSU reads for {wildcards.sample}"
    input:
        interleaved_reads=rules.reformat_rRNA_for_bbduk.output.interleaved_reads,
        ref=config["databases"]["sortmeRNA_ssu"]
    output:
        ssu_interleaved=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.interleaved.rRNA.fastq.gz",
    log:
        stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbduk_ssu.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbduk_ssu.txt"
    params:
        options="interleaved=t -Xmx40g k=31 hdist=0"
    threads:
        config["RNA"]["bbduk"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.ssu_interleaved})

        bbduk.sh \
            in={input.interleaved_reads} \
            outm={output.ssu_interleaved} \
            ref={input.ref} \
            threads={threads} \
            {params.options}  \
            > {log.stdout} 2>&1
        """

rule reformat_rRNA_for_SSU:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBDuk] Reformatting ribodetected seperated reads for {wildcards.sample} to ensure read order for downstream processing"
    input:
        interleaved_reads=rules.bbduk_ssu.output.ssu_interleaved,
    output:
        rRNA_ssu_r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.rRNA.R1.fastq.gz",
        rRNA_ssu_r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_SSU.rRNA.R2.fastq.gz",
    log:
        stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbduk_ssu.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbduk_ssu.txt"
    params:
        options="-Xmx40g"     
    threads:
        config["RNA"]["bbduk"]["threads"]
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.rRNA_ssu_r1})

        reformat.sh in={input.interleaved_reads} threads={threads} out={output.rRNA_ssu_r1} out2={output.rRNA_ssu_r2} {params.options} > {log.stdout} 2>&1
        """
        
# reformat.sh in=ANN_10.rRNA.r1.fastq.gz in2=ANN_10.rRNA.r2.fastq.gz interleaved=t threads=6 out=interleaved.fq.gz -Xmx40g
# bbduk.sh in=interleaved.fq.gz ref=/data_2/Databases/SILVA_138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta interleaved=t outm=matched.fq.gz threads=6 k=31 hdist=0 -Xmx40g
# reformat.sh in=matched.fq.gz threads=6 out=ANN_10.rRNA.SSU.r1.fastq.gz out2=ANN_10.rRNA.SSU.r2.fastq.gz -Xmx40g



## BBDUK is somehow not sorted when using threads as the chunks are out of order, and thus the read pairs are not in the same order as the input, which causes problems for downstream processing. To solve this, we first reformat the ribodetected reads to be interleaved, then run bbduk on the interleaved reads, and finally reformat the output back to paired-end format. This ensures that the read pairs are in the same order as the input, and thus can be processed correctly downstream.
#bbduk.sh in1=/data/rasmus/Dev/SnakeResDev/RNA/ANN_10/ribodetector/ANN_10.rRNA.r1.fastq.gz in2=/data/rasmus/Dev/SnakeResDev/RNA/ANN_10/ribodetector/ANN_10.rRNA.r2.fastq.gz -Xmx40g outm1=/data/rasmus/Dev/SnakeResDev/rRNA/ANN_10/filtered/ANN_10_R1.rRNA.fastq.gz outm2=/data/rasmus/Dev/SnakeResDev/rRNA/ANN_10/filtered/ANN_10_R2.rRNA.fastq.gz ref=/data_2/Databases/SILVA_138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta k=31 hdist=0