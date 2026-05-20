
#bbmap.sh ref=combined.fa in=reads.fq out=mapped.sam ambig=best
rule bbmap:
    conda:
        "../envs/bbmap.yaml"
    message:
        "[BBMap] Running BBMap for {wildcards.sample}"
    input:
        cleaned_r1=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R1.cleaned.fastq.gz",
        cleaned_r2=f"{RESULTS_DIR}/qc/{{sample}}/decontamination/{{sample}}_R2.cleaned.fastq.gz",
        ref=config["databases"]["sortmeRNA_ssu"]
    output:
        aligned=f"{RESULTS_DIR}/RNA/{{sample}}/bbmap/{{sample}}_aligned.sam",
    log:
        stdout = f"{RESULTS_DIR}/RNA/{{sample}}/logs/bbmap.log"
    benchmark:
        f"{RESULTS_DIR}/RNA/{{sample}}/benchmarks/bbmap.txt"
    params:
        options=config["RNA"]["bbmap"]["options"]
    threads:
        config["RNA"]["bbmap"]["threads"]
    shell:
        r"""
        set -euo pipefail
        
        bbmap.sh \
            in1={input.cleaned_r1} \
            in2={input.cleaned_r2} \
            ref={input.ref} \
            out={output.aligned} \
            threads={threads} \
            {params.options} \
            > {log.stdout} 2>&1
        """


java -ea   -Xmx384453m -Xms384453m -cp /data/rasmus/conda-envs/RAHCS_env/opt/bbmap-39.81-1/bbtools.jar align2.BBMap build=1 overwrite=true fastareadlen=500 in1=A263.rRNA.r1.fastq.gz in2=A263.rRNA.r2.fastq.gz ref=/data_2/Databases/SILVA_138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta out=mapped.sam threads=16 maxindel=20 strictmaxindel=t minid=0.90 pairlen=1000 rescuedist=1200 minhits=2
Executing align2.BBMap [build=1, overwrite=true, fastareadlen=500, in1=A263.rRNA.r1.fastq.gz, in2=A263.rRNA.r2.fastq.gz, ref=/data_2/Databases/SILVA_138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta, out=mapped.sam, threads=16, maxindel=20, strictmaxindel=t, minid=0.90, pairlen=1000, rescuedist=1200, minhits=2]
Version 39.81
