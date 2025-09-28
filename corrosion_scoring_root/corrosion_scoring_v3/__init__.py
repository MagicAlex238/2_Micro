"""
Refactored Corrosion Scoring System v3.0

This module provides a clean separation between text mining, scoring, and synergy detection
with improved error handling and maintainable architecture.:)
"""
from .text_miner import TextMiner
from .score_calculator import ScoreCalculator
from .synergy_detector import SynergyDetector
from .config import ScoringConfig

from .exceptions import ScoringError, TextMiningError, SynergyDetectionError
from .global_terms import (
    functional_categories_dict, metal_terms_dict, corrosion_synergies_dict, # scoring dictionaries
    metal_mapping, # cannonical mapping
    mechanisms_dict, pathway_dict, operational_environmental_factors_dict) # additional dictionaries)
from .validators import ValidationError
from .term_processor import TermProcessor


# Convenience function for backward compatibility
def calculate_overall_scores(text: str, processors: dict, config: ScoringConfig = None, brenda_metals: list = None):
    """
    Backward-compatible wrapper for the new scoring system.
    
    Args:
        text: Text to analyze
        processors: Dictionary containing processors
        config: Scoring configuration
        brenda_metals: List of metals from BRENDA
        
    Returns:
        Dictionary containing all scores and analysis results
    """
    # Simplified wrapper: accept processors as a dict (or None). Merge any provided
    # processors into TextMiner's internal default processors so population-only
    # processors (pathway, mechanisms, operational) are preserved.
    if config is None:
        config = ScoringConfig()

    # Initialize components
    text_miner = TextMiner(config)

    if processors is not None:
        if not isinstance(processors, dict):
            raise TypeError("processors must be a dict mapping names to TermProcessor instances or None")
        # Merge: keep TextMiner's defaults but override with provided processors
        merged = dict(text_miner.processors)
        merged.update(processors)
        text_miner.processors = merged
    score_calculator = ScoreCalculator(config)
    synergy_detector = SynergyDetector(config)

    # Extract features
    features = text_miner.extract_all_features(text, brenda_metals)

    # Calculate scores
    scores = score_calculator.calculate_scores(features)

    # Detect synergies
    synergies = synergy_detector.detect_synergies(features)

    # Combine results
    return {**features, **scores, **synergies}


def get_default_processors():
    """Return the canonical processors dict used by example notebooks.

    This factory creates TermProcessor instances from the package's global
    dictionaries and returns a mapping similar to what notebooks expect.
    """
    return {
        'fc_processor': TermProcessor(functional_categories_dict),
        'metal_processor': TermProcessor(metal_terms_dict),
        'synergy_processor': TermProcessor(corrosion_synergies_dict),
    }

__all__ = [
    'TextMiner',
    'ScoreCalculator', 
    'SynergyDetector',
    'ScoringConfig',
    'ScoringError',
    'TextMiningError',
    'SynergyDetectionError',
    'calculate_overall_scores',
    # re-exported globals as names
    'functional_categories_dict',
    'metal_terms_dict',
    'corrosion_synergies_dict',
    'metal_mapping',
    'pathway_dict',
    'mechanisms_dict',
    'operational_environmental_factors_dict',
    'ValidationError',
    'TermProcessor'
]
