import os
import pandas as pd
import yaml
from typing import Dict, Tuple, List, Optional

### ANALYSIS SPECIFIC FUNCTIONALITIES

def fastqc_input_path(
    sample: str,
    stage: str,
    readtag: str,
    sample_to_illumina: Dict[str, List[str]],
    illumina_trim_read_path: str,
) -> str:
    """
    Resolve which FASTQ FastQC should run on, given:
      stage: "raw" or "trimmed"
      readtag: "1","2","1.trim","2.trim"
    """
    if stage == "raw":
        if readtag not in ("1", "2"):
            raise ValueError(f"Invalid readtag for raw stage: {readtag}")
        idx = 0 if readtag == "1" else 1
        return sample_to_illumina[sample][idx]

    if stage == "trimmed":
        if readtag not in ("1.trim", "2.trim"):
            raise ValueError(f"Invalid readtag for trimmed stage: {readtag}")
        return f"{illumina_trim_read_path}/{sample}_{readtag}.fastq.gz"

    raise ValueError(f"Unknown stage: {stage}")

### DEFINE AND LIST THE EXPECTED RESULTS
def read_results_catalogue(results_catalogue_path: str) -> Dict[str, str]:
    with open(results_catalogue_path, "r") as f:
        return yaml.safe_load(f) or {}

def list_results(
    samples: List[str],
    results_catalogue: Dict,
    context: Optional[Dict[str, str]] = None,
) -> List[str]:
    """
    Supports two catalogue formats:

    Old:
      fastp:
        - "{sample}_1.trim.fastq.gz"

    New:
      fastp:
        base: "{illumina_trim}"
        files:
          - "{sample}_1.trim.fastq.gz"
    """
    context = context or {}
    out: List[str] = []

    for step, spec in results_catalogue.items():

        # --- normalize to (base_template, patterns_list) ---
        if isinstance(spec, dict):
            base_t = spec.get("base", "") or ""
            patterns = spec.get("files", []) or []
        else:
            base_t = ""
            patterns = spec or []

        if isinstance(patterns, str):
            patterns = [patterns]

        for pat in patterns:
            base_has_sample = "{sample}" in base_t
            pat_has_sample = "{sample}" in pat

            if base_has_sample or pat_has_sample:
                for s in samples:
                    base = base_t.format(sample=s, **context) if base_t else ""
                    full = os.path.join(base, pat) if base else pat
                    out.append(full.format(sample=s, **context))
            else:
                base = base_t.format(**context) if base_t else ""
                full = os.path.join(base, pat) if base else pat
                out.append(full.format(**context))

    # de-duplicate while keeping order
    seen = set()
    out_unique = []
    for x in out:
        if x not in seen:
            out_unique.append(x)
            seen.add(x)

    return out_unique

def sample_read_map(
    samplesheet: pd.DataFrame,
    sample_id: str = "sample_id",
    sample_name: str = "sample_name",
    illumina_reads: str = "illumina_reads",
    illumina_read_config_path: str = "",
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:

    samples: Dict[str, str] = {}
    sample_to_illumina: Dict[str, List[str]] = {}

    for _, row in samplesheet.iterrows():
        sample = str(row[sample_id]).strip()

        # Split the *value* in the illumina_reads column (do not overwrite the column-name variable)
        reads = [r.strip() for r in str(row[illumina_reads]).split(",")]

        if len(reads) != 2:
            raise ValueError(f"Expected 2 read files for sample '{sample}', got: {reads}")

        illumina_full = [os.path.join(illumina_read_config_path, r) for r in reads]

        samples[sample] = sample
        sample_to_illumina[sample] = illumina_full

    return samples, sample_to_illumina
