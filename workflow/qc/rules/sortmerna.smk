# change it to specify config bechaviour as opposed to hardcoded SSU and LSU in this script, so now 
# i define it as stages, such that i can change order, use other databases, add more databases like archea or eukaryotic. 
# I also added a workdir for sortmerna, which is needed for large databases to avoid memory issues.

rule sortmerna:
    conda:
        "../envs/sortmerna.yaml"
    message:
        "[SortMeRNA] seperate reads according to the databases from stages {wildcards.stage} for {wildcards.sample}"
    input:
        r1=sortmerna_r1,
        r2=sortmerna_r2,
        db=sortmerna_db,
        db_idx=sortmerna_db_idx
    output:
        aligned_fwd="{outdir}/{sample}/QC/sortmerna/{stage}/{sample}_{stage}.aligned_fwd.fq.gz",
        aligned_rev="{outdir}/{sample}/QC/sortmerna/{stage}/{sample}_{stage}.aligned_rev.fq.gz",
        nonaligned_fwd="{outdir}/{sample}/QC/sortmerna/{stage}/{sample}_{stage}.nonaligned_fwd.fq.gz",
        nonaligned_rev="{outdir}/{sample}/QC/sortmerna/{stage}/{sample}_{stage}.nonaligned_rev.fq.gz",
    log:
        stdout = "{outdir}/{sample}/logs/sortmerna_{stage}.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/sortmerna_{stage}.txt"
    params:
        workdir="{outdir}/{sample}/QC/sortmerna/workdir/{stage}",
        options=config["qc"]["sortmerna"]["options"],
        aligned_prefix="{outdir}/{sample}/QC/sortmerna/{stage}/{sample}_{stage}.aligned",
        nonaligned_prefix="{outdir}/{sample}/QC/sortmerna/{stage}/{sample}_{stage}.nonaligned"
    threads:
        config["qc"]["sortmerna"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir}

        sortmerna \
            --ref {input.db} \
            --idx-dir {input.db_idx} \
            --workdir {params.workdir} \
            --threads {threads} \
            {params.options} \
            --reads {input.r1} \
            --reads {input.r2} \
            --aligned {params.aligned_prefix} \
            --other {params.nonaligned_prefix} \
            > {log.stdout} 2>&1

        rm -rf {params.workdir}/kvdb {params.workdir}/readb
        """
#/usr/bin/time -v sortmerna --ref /data_2/Databases/SILVA_138/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta --ref /data_2/Databases/SILVA_138/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta --idx-dir /data_2/Databases/sortmerna_idx/LSU/ --idx-dir /data_2/Databases/sortmerna_idx/SSU/ --workdir sortmeRNA  --fastx --paired_in --out2 --threads 10 --aligned sortmeRNA/aligned.fastq --other sortmeRNA/notaligned.fastq --reads decontamination/ANN_10_R1.cleaned.fastq.gz --reads decontamination/ANN_10_R2.cleaned.fastq.gz --otu_map ANN_10.otu