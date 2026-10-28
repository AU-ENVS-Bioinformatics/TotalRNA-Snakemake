"""
Helpers for viral sequence analysis.

Validates dependencies required by viral postprocessing, which reuses
TransDecoder ORFs, eggNOG annotations, and Salmon transcript abundance.

Expected completion outputs remain in targets/viral_sequence_analysis.py.
"""

def validate_viral_config(config: dict) -> None:
    """Validate dependencies required by viral postprocessing."""
    viral_config = config["viral_sequence_analysis"]

    if not viral_config.get("postprocessing", True):
        return

    functional_config = config["functional_profiling"]

    if functional_config["module"] not in {
        "assembly",
        "both",
        "all",
    }:
        raise ValueError(
            "Viral postprocessing requires assembly-based "
            "functional profiling."
        )

    if functional_config["ORF_predictor"] != "transdecoder":
        raise ValueError(
            "Viral postprocessing requires "
            "functional_profiling.ORF_predictor: transdecoder."
        )