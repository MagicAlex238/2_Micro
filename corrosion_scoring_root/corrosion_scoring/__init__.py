"""
Refactored Corrosion Scoring System v2.0

This module provides a clean separation between text mining, scoring, and synergy detection
with improved error handling and maintainable architecture.
"""

from .text_miner import TextMiner
from .score_calculator import ScoreCalculator
from .synergy_detector import SynergyDetector
from .config import ScoringConfig
from .exceptions import ScoringError, TextMiningError, SynergyDetectionError

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
    if config is None:
        config = ScoringConfig()
    
    # Initialize components
    text_miner = TextMiner(processors)
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

__all__ = [
    'TextMiner',
    'ScoreCalculator', 
    'SynergyDetector',
    'ScoringConfig',
    'ScoringError',
    'TextMiningError',
    'SynergyDetectionError',
    'calculate_overall_scores'
]
