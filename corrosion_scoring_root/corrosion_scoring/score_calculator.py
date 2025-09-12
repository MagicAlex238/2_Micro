"""
ScoreCalculator class for computing scores from extracted features.
"""

import math
from typing import Dict, Any, List, Tuple
from .config import ScoringConfig
from .exceptions import ScoringError
from .validators import ValidationError

class ScoreCalculator:
    """
    Handles scoring calculations for corrosion relevance assessment.
    
    This class takes extracted features and calculates various scores
    without performing any text mining operations.
    """
    
    def __init__(self, config: ScoringConfig = None):
        """
        Initialize the ScoreCalculator.
        
        Args:
            config: Configuration settings for scoring
        """
        self.config = config or ScoringConfig()
        
        # Import global dictionaries safely
        try:
            from .global_terms import (
                functional_categories_dict,
                corrosion_synergies_dict
            )
            self.functional_categories_dict = functional_categories_dict
            self.corrosion_synergies_dict = corrosion_synergies_dict
        except ImportError as e:
            raise ScoringError(f"Failed to import global terms for scoring: {e}") from e
    
    def calculate_scores(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate all scores from extracted features.
        
        Args:
            features: Dictionary containing extracted features
            
        Returns:
            Dictionary containing calculated scores
            
        Raises:
            ScoringError: If score calculation fails
        """
        try:
            scores = {}
            
            # Calculate functional category scores
            fc_scores = self._calculate_functional_scores(features)
            scores.update(fc_scores)
            
            # Calculate metal scores
            metal_scores = self._calculate_metal_scores(features)
            scores.update(metal_scores)
            
            # Calculate legacy synergy scores (for fallback)
            legacy_synergy_scores = self._calculate_legacy_synergy_scores(features)
            scores.update(legacy_synergy_scores)
            
            # Calculate overall scores
            overall_scores = self._calculate_overall_scores(scores)
            scores.update(overall_scores)
            
            # Calculate final corrosion relevance score
            scores['corrosion_relevance_score'] = self._calculate_corrosion_relevance_score(scores)
            
            return scores
            
        except Exception as e:
            raise ScoringError(f"Score calculation failed: {e}") from e
    
    def _calculate_functional_scores(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate functional category scores."""
        try:
            fc_matches = features.get('functional_matches_detailed', {})
            functional_score = 0.0
            functional_categories_scored = []
            
            for category, terms in fc_matches.items():
                if not terms:
                    continue
                
                # Calculate score for this category
                category_hits = len(set(terms))
                base_score = math.log(category_hits + 1)
                
                # Get weight from global dictionary
                weight = 1.0  # Default weight
                if category in self.functional_categories_dict:
                    cat_data = self.functional_categories_dict[category]
                    if isinstance(cat_data, dict) and 'score' in cat_data:
                        weight = float(cat_data['score'])
                
                weighted_score = base_score * weight
                functional_score += weighted_score
                
                functional_categories_scored.append({
                    'category': category,
                    'score': weighted_score,
                    'terms': sorted(set(terms)),
                    'hit_count': category_hits,
                    'weight': weight
                })
            
            return {
                'functional_score': float(functional_score),
                'functional_categories': functional_categories_scored
            }
            
        except Exception as e:
            raise ScoringError(f"Functional score calculation failed: {e}") from e
    
    def _calculate_metal_scores(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate metal-related scores."""
        try:
            metal_matches = features.get('metal_matches_detailed', {})
            consolidated_metals = features.get('consolidated_metals', [])
            
            # Count unique metal terms detected
            total_metal_terms = sum(len(set(terms)) for terms in metal_matches.values())
            metal_score = float(total_metal_terms)
            
            return {
                'metal_score': metal_score,
                'metal_count': len(consolidated_metals),
                'metal_categories_detected': list(metal_matches.keys())
            }
            
        except Exception as e:
            raise ScoringError(f"Metal score calculation failed: {e}") from e
    
    def _calculate_legacy_synergy_scores(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate legacy keyword-based synergy scores for fallback."""
        try:
            # This is a simplified version for legacy compatibility
            # The main synergy detection is handled by SynergyDetector
            
            return {
                'legacy_synergy_score': 0.0,
                'legacy_synergy_groups': []
            }
            
        except Exception as e:
            raise ScoringError(f"Legacy synergy score calculation failed: {e}") from e
    
    def _calculate_overall_scores(self, scores: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate weighted overall scores."""
        try:
            functional_score = scores.get('functional_score', 0.0)
            metal_score = scores.get('metal_score', 0.0)
            
            overall_functional_score = functional_score * self.config.functional_score_weight
            overall_metal_score = metal_score * self.config.metal_score_weight
            
            return {
                'overall_functional_score': float(overall_functional_score),
                'overall_metal_score': float(overall_metal_score)
            }
            
        except Exception as e:
            raise ScoringError(f"Overall score calculation failed: {e}") from e
    
    def _calculate_corrosion_relevance_score(self, scores: Dict[str, Any]) -> float:
        """Calculate the final corrosion relevance score."""
        try:
            overall_functional = scores.get('overall_functional_score', 0.0)
            overall_metal = scores.get('overall_metal_score', 0.0)
            overall_synergy = scores.get('overall_synergy_score', 0.0)
            
            return float(overall_functional + overall_metal + overall_synergy)
            
        except Exception as e:
            raise ScoringError(f"Corrosion relevance score calculation failed: {e}") from e