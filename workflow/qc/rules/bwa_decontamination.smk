rule bwa_contaminants:
    conda:
        "../envs/bwa.yaml"
    message:
        "[BWA] decontaminate {wildcards.sample}"
    input:
        trim_r1="{outdir}/{sample}/QC/trimmed/{sample}_R1.fastq.gz",
        trim_r2="{outdir}/{sample}/QC/trimmed/{sample}_R2.fastq.gz"
    output:
        aligned="{outdir}/{sample}/QC/decontamination/{sample}_aligned_contaminants.sam",
    log:
        stdout="{outdir}/{sample}/logs/aligned_contaminants.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/{sample}_bwa_aligned_contaminants.txt"
    params:
        genome=config["databases"]["genome_for_decontamination"],
    threads: config["qc"]["bwa_contaminants"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.aligned})

        bwa-mem2 mem -t {threads} {params.genome} {input.trim_r1} {input.trim_r2} -o {output.aligned} > {log.stdout} 2>&1
        """

rule contamination_source_ids:
    conda:
        "../envs/bwa.yaml"
    message:
        "[BWA] decontaminate {wildcards.sample}"
    input:
        aligned="{outdir}/{sample}/QC/decontamination/{sample}_aligned_contaminants.sam",
    output:
        contaminants_id="{outdir}/{sample}/QC/decontamination/{sample}_contaminants_id.txt",
    log:
        stdout="{outdir}/{sample}/logs/contamination_source_ids.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/contamination_source_ids.txt"
    params:
        filters=config["qc"]["contamination_source_ids"]["options"],
        script=f"{workflow.basedir}/qc/scripts/dehosting.sh",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.contaminants_id})
        
        {params.script} \
          -i {input.aligned} \
          -o {output.contaminants_id} \
          {params.filters} > {log.stdout} 2>&1
        """

#Consider that the decontamination is currently limitied to paired-end reads, as the contaminant read ids are extracted from the aligned paired-end reads, and then used to filter both R1 and R2. 
#If we want to decontaminate single-end reads, we would need to modify the workflow to handle single-end data separately, and ensure that the read ids are correctly extracted and used for filtering.

rule decontamination:
    conda:
        "../envs/seqkit.yaml"
    message:
        "[Seqkit] decontaminate {wildcards.sample}"
    input:
        trim_r1="{outdir}/{sample}/QC/trimmed/{sample}_R1.fastq.gz",
        trim_r2="{outdir}/{sample}/QC/trimmed/{sample}_R2.fastq.gz",
        read_ids="{outdir}/{sample}/QC/decontamination/{sample}_contaminants_id.txt"
    output:
        cleaned_r1="{outdir}/{sample}/QC/decontamination/{sample}_R1.cleaned.fastq.gz",
        cleaned_r2="{outdir}/{sample}/QC/decontamination/{sample}_R2.cleaned.fastq.gz",
    log:
        stdout="{outdir}/{sample}/logs/decontamination.log"
    benchmark:
        "{outdir}/{sample}/benchmarks/{sample}_decontamination.txt"
    threads: config["qc"]["decontamination"]["threads"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cleaned_r1})

        seqkit grep -f {input.read_ids} -v --threads {threads} {input.trim_r1} -o {output.cleaned_r1} 2> {log.stdout}

        seqkit grep -f {input.read_ids} -v --threads {threads} {input.trim_r2} -o {output.cleaned_r2} 2>> {log.stdout}
        """