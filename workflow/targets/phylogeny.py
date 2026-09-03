from typing import List
from pathlib import Path


def phylogeny_outputs(
    results_dir: Path,
    samples: List[str],
    config: dict,
) -> List[str]:
    """
    Return final expected SSU phylogeny outputs.
    """

    return [
        f"{results_dir}/Phylogeny/tree/cross_sample_SSU.treefile",
        f"{results_dir}/Phylogeny/tree/cross_sample_SSU.iqtree",
        f"{results_dir}/Phylogeny/tree/cross_sample_SSU.model.gz",
    ]