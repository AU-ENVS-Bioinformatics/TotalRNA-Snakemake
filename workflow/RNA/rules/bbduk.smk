import subprocess

def decide_threads_lines(input):
    """Return 1 thread if <100k reads, else use configured threads."""
    fq = input.r1
    cmd = f"zcat {fq} | wc -l"
    lines = int(subprocess.check_output(cmd, shell=True))
    reads = lines // 4

    if reads < 100000:
        return 1
    else:
        return config["RNA"]["bbduk"]["threads"]

def decide_threads_size(input):
    size_mb = os.path.getsize(input.r1) / 1e6

    if size_mb < 50:   # tune threshold
        return 1
    else:
        return config["RNA"]["bbduk"]["threads"]


if config["RNA"]["method"] == "ribodetector" and config["RNA"]["refinement"] == "bbduk":
    rule bbduk_ssu:
        conda:
            "../envs/bbmap.yaml"
        message:
            "[BBDuk] Extracting SSU reads for {wildcards.sample}"
        input:
            r1=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r1.fastq.gz",
            r2=f"{RESULTS_DIR}/RNA/{{sample}}/ribodetector/{{sample}}.rRNA.r2.fastq.gz",
            ref=config["databases"]["sortmeRNA_ssu"]
        output:
            ssu_r1=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_R1.rRNA.fastq.gz",
            ssu_r2=f"{RESULTS_DIR}/rRNA/{{sample}}/filtered/{{sample}}_R2.rRNA.fastq.gz",
        log:
            stdout=f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbduk_ssu.log"
        benchmark:
            f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbduk_ssu.txt"
        params:
            options=config["RNA"]["bbduk"]["options"]        
        threads:
            decide_threads_size(input)
        shell:
            r"""
            set -euo pipefail

            mkdir -p $(dirname {output.ssu_r1})

            bbduk.sh \
                in1={input.r1} \
                in2={input.r2} \
                outm1={output.ssu_r1} \
                outm2={output.ssu_r2} \
                ref={input.ref} \
                threads={threads} \
                -Xmx40g k=31 hdist=0\
                > {log.stdout} 2>&1
            """

#bbduk.sh in1=/data/rasmus/Dev/SnakeResDev/RNA/ANN_10/ribodetector/ANN_10.rRNA.r1.fastq.gz in2=/data/rasmus/Dev/SnakeResDev/RNA/ANN_10/ribodetector/ANN_10.rRNA.r2.fastq.gz -Xmx40g outm1=/data/rasmus/Dev/SnakeResDev/rRNA/ANN_10/filtered/ANN_10_R1.rRNA.fastq.gz outm2=/data/rasmus/Dev/SnakeResDev/rRNA/ANN_10/filtered/ANN_10_R2.rRNA.fastq.gz ref=/data_2/Databases/SILVA_138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta k=31 hdist=0