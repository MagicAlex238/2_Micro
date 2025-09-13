"""
Refactored Corrosion Scoring System v2.0

This module provides similar scoring system to the first version where the text mining, scoring, and synergy detection
where incorporated between each on the architecture.
Unfortunately the refractoring took several months and was not suscessfull in retrieving the output that version 1 suscessfully retrieved

"""
# ========== __init__.py ==========
# corrosion_scoring/__init__.py

from .term_processor import TermProcessor
from .scoring_system import (
    calculate_overall_scores,
    calculate_corrosion_relevance_score,
    assign_mechanism_from_pathway,
    assign_corrosion_mechanisms,
    infer_mechanisms_from_pathway_category,
    consolidate_metal_terms,
    validate_against_pathways,
)

from .utils_ec import normalize_ec_id, strip_all_ec_tokens
from .name_utils import enhanced_clean_protein_name, clean_protein_name

# Re-export global dictionaries so they are available as `cs.metal_mapping`, etc.
from .global_terms import (
        metal_terms_dict,
        corrosion_synergies_dict,
        functional_categories_dict,
        metal_mapping, pathway_dict, mechanisms_dict, operational_environmental_factors_dict # Not for scoring
    )

__all__ = [
    "TermProcessor",
    "calculate_overall_scores",
    "calculate_corrosion_relevance_score",
    "assign_mechanism_from_pathway",
    "assign_corrosion_mechanisms",
    "infer_mechanisms_from_pathway_category",
    "consolidate_metal_terms",
    "validate_against_pathways",
    "normalize_ec_id",
    "strip_all_ec_tokens",
    "enhanced_clean_protein_name",
    "clean_protein_name",
    # re-exported globals as names
    "metal_terms_dict",
    "corrosion_synergies_dict",
    "functional_categories_dict",
    "metal_mapping",
    "pathway_dict",
    "mechanisms_dict",
    "operational_environmental_factors_dict",
]