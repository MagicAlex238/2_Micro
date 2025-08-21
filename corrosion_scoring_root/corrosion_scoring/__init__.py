# corrosion_scoring/__init__.py

from .term_processor import TermProcessor
from .scoring_system import (
    calculate_overall_scores,
    calculate_corrosion_relevance_score,
    assign_mechanism_from_pathway,
    assign_corrosion_mechanisms,
)
from .utils_ec import normalize_ec_id, strip_all_ec_tokens
from .name_utils import enhanced_clean_protein_name

__all__ = [
    "TermProcessor",
    "calculate_overall_scores",
    "calculate_corrosion_relevance_score",
    "assign_mechanism_from_pathway",
    "assign_corrosion_mechanisms",
    "normalize_ec_id",
    "strip_all_ec_tokens",
    "enhanced_clean_protein_name",
]