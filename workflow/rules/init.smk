"""
Data Preparation Module
=======================
Handles initial data processing steps including:
- Database preparation for SortMeRNA
- Quality trimming with Trim Galore
- rRNA/tRNA filtering with SortMeRNA
"""

import os
import re
from dataclasses import dataclass
from collections import defaultdict


@dataclass(frozen=True)
class Sample:
    name: str
    r1: str
    r2: str


def build_samples(raw_dir, pattern=None, normalize_fn=None):
    files = os.listdir(raw_dir)
    tmp = defaultdict(dict)

    if pattern is None:
        pattern = re.compile(r"(.+?)[-_]?(?:R|read)?([12])\.fastq\.gz$")
    if normalize_fn is None:
        normalize_fn = lambda x: x.replace("-", "_")

    for f in files:
        match = pattern.match(f)
        if not match:
            continue
        raw_sample = match.group(1)
        read = match.group(2)

        clean_sample = normalize_fn(raw_sample)

        tmp[clean_sample][f"R{read}"] = os.path.join(raw_dir, f)

    samples = {}

    for name, reads in tmp.items():
        if "R1" not in reads or "R2" not in reads:
            raise ValueError(f"Sample {name} is missing R1 or R2")

        samples[name] = Sample(name, reads["R1"], reads["R2"])

    return samples


samples_dict = build_samples("reads/")
unique_samples = sorted(samples_dict.keys())
unique_samples = unique_samples[:2]  # Limit to first 10 samples for testing


onstart:
    print("#### Total-RNA workflow")
    print("Checking for required software...")
    shell("type usearch || {{ echo 'usearch not found'; }}")
    if len(unique_samples) == 0:
        print(
            "There are no samples name detected! Please check your samples file and config files."
        )
    print(f"Detected {len(unique_samples)} samples: {unique_samples}")
    print("Reading samples and config files...")


rule trim_files:
    input:
        expand(
            "results/trim_galore/{sample}_R1.fq.gz",
            sample=unique_samples,
        ),
        expand(
            "results/trim_galore/{sample}_R2.fq.gz",
            sample=unique_samples,
        ),


rule sortmerna_rRNA:
    input:
        expand(
            "results/sortmerna/SSU/{sample}_{dir}.fq.gz",
            sample=unique_samples,
            dir=["fwd", "rev"],
        ),


rule sortmerna_mRNA:
    input:
        expand(
            "results/sortmerna/not_LSU/{sample}_{dir}.fq.gz",
            sample=unique_samples,
            dir=["fwd", "rev"],
        ),


include: "../rules/trim_galore.smk"
include: "../rules/sortmerna.smk"
