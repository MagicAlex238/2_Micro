# corrosion_scoring/__init__.py

from .term_processor import TermProcessor
from .scoring_system import (
    calculate_overall_scores,
    calculate_corrosion_relevance_score,
    assign_mechanism_from_pathway,
    assign_corrosion_mechanisms,
    infer_mechanisms_from_pathway_label,
    consolidate_metal_terms,
    validate_against_pathways
)

from .utils_ec import normalize_ec_id, strip_all_ec_tokens
from .name_utils import enhanced_clean_protein_name, clean_protein_name

# Re-export global dictionaries so they are available as `cs.metal_mapping`, etc.
from .global_terms import (
        metal_terms_dict, # For retrieval and scoring
        corrosion_synergies_dict, # For retrieval and scoring
        functional_categories_dict, # For retrieval and scoring
        pathway_dict, mechanisms_dict, operational_environmental_factors_dict, # for retrieval but Not for scoring just for poppulating
        metal_mapping # no for retrieval
    )

__all__ = [
    "TermProcessor",
    "calculate_overall_scores",
    "calculate_corrosion_relevance_score",
    "assign_mechanism_from_pathway",
    "assign_corrosion_mechanisms",
    "infer_mechanisms_from_pathway_label", 
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