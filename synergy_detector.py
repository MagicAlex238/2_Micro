"""
SynergyDetector class for detecting corrosion synergies.
"""

from typing import Dict, Any, List, Set, Tuple, Optional
from .config import ScoringConfig
from .exceptions import SynergyDetectionError

class SynergyDetector:
    """
    Handles detection and scoring of corrosion synergies.
    
    This class analyzes functional category co-occurrence patterns
    to identify synergistic corrosion mechanisms.
    """
    
    def __init__(self, config: ScoringConfig = None):
        """
        Initialize the SynergyDetector.
        
        Args:
            config: Configuration settings for synergy detection
        """
        self.config = config or ScoringConfig()
        
        # Functional categories to check for co-occurrence
        self.subcategories_fc = [
            "o2_consumption", "nitrogen_metabolism", "iron_metabolism",
            "sulfur_metabolism", "h2_consumption", "direct_eet",
            "carbon_metabolism", "indirect_eet", "organic_acid_metabolism",
            "metal_binding_chelation", "biofilm_formation", "manganese_processes",
            "methanogenesis", "fumarate_formation", "halogen_related",
            "phosphorus_metabolism"
        ]
    
    def detect_synergies(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """
        Detect synergies from extracted features.
        
        Args:
            features: Dictionary containing extracted features
            
        Returns:
            Dictionary containing synergy detection results
            
        Raises:
            SynergyDetectionError: If synergy detection fails
        """
        try:
            # Get functional category matches
            fc_matches = features.get('functional_matches_detailed', {})
            
            # Detect functional category synergies
            fc_synergy_results = self._detect_functional_category_synergies(fc_matches)
            
            # Calculate overall synergy score with weighting
            overall_synergy_score = fc_synergy_results['synergy_score'] * self.config.synergy_score_weight
            
            return {
                **fc_synergy_results,
                'overall_synergy_score': float(overall_synergy_score)
            }
            
        except Exception as e:
            raise SynergyDetectionError(f"Synergy detection failed: {e}") from e
    
    def _detect_functional_category_synergies(self, fc_matches: Dict[str, List[str]]) -> Dict[str, Any]:
        """
        Detect synergies based on co-occurrence of terms from different functional categories.
        
        Args:
            fc_matches: Dictionary of functional category matches
            
        Returns:
            Dictionary containing synergy detection results
        """
        try:
            # Step 1: Filter relevant categories and collect their terms
            detected_categories = {}
            all_found_terms = set()
            
            for category in self.subcategories_fc:
                if category in fc_matches and fc_matches[category]:
                    detected_categories[category] = list(fc_matches[category])
                    all_found_terms.update(fc_matches[category])
            
            # Step 2: Initialize synergy results
            synergy_results = {
                'fc_cooccurrence_synergy_hit': False,
                'synergy_score': 0.0,
                'synergy_child_terms_found': [],
                'synergy_categories_involved': [],
                'synergy_description': '',
                'synergy_type': 'none'
            }
            
            # Step 3: Check minimum requirements
            if len(detected_categories) < self.config.min_synergy_categories:
                return synergy_results
            
            # Step 4: Check for priority synergies
            max_synergy_score = 0.0
            best_synergy = None
            involved_categories = []
            
            for synergy_pair, synergy_info in self.config.priority_synergies.items():
                cat1, cat2 = synergy_pair
                if cat1 in detected_categories and cat2 in detected_categories:
                    # Calculate combined term count from both categories
                    combined_terms = set(detected_categories[cat1] + detected_categories[cat2])
                    
                    # Require minimum terms for high-confidence synergy
                    if len(combined_terms) >= self.config.min_synergy_terms:
                        current_score = synergy_info['score']
                        if current_score > max_synergy_score:
                            max_synergy_score = current_score
                            best_synergy = synergy_info
                            involved_categories = [cat1, cat2]
            
            # Step 5: Check for general multi-category synergy if no priority synergy found
            if max_synergy_score == 0.0 and len(detected_categories) >= self.config.min_synergy_categories:
                if len(all_found_terms) >= self.config.min_general_synergy_terms:
                    max_synergy_score = 1.5  # Lower score for general synergy
                    involved_categories = list(detected_categories.keys())
                    best_synergy = {
                        'description': f'Multi-pathway Synergy ({len(detected_categories)} categories)'
                    }
            
            # Step 6: Populate results if synergy detected
            if max_synergy_score > 0.0:
                synergy_terms = set()
                for cat in involved_categories:
                    if cat in detected_categories:
                        synergy_terms.update(detected_categories[cat])
                
                synergy_results.update({
                    'fc_cooccurrence_synergy_hit': True,
                    'synergy_score': max_synergy_score,
                    'synergy_child_terms_found': sorted(list(synergy_terms)),
                    'synergy_categories_involved': involved_categories,
                    'synergy_description': best_synergy['description'],
                    'synergy_type': 'functional_category_cooccurrence'
                })
            
            return synergy_results
            
        except Exception as e:
            raise SynergyDetectionError(f"Functional category synergy detection failed: {e}") from e
    
    def get_synergy_explanation(self, synergy_categories: List[str]) -> str:
        """
        Get a detailed explanation of detected synergies.
        
        Args:
            synergy_categories: List of involved synergy categories
            
        Returns:
            Detailed explanation string
        """
        if not synergy_categories or len(synergy_categories) < 2:
            return "No significant synergies detected."
        
        # Check if it's a known priority synergy
        if len(synergy_categories) == 2:
            synergy_pair = tuple(sorted(synergy_categories))
            if synergy_pair in self.config.priority_synergies:
                return self.config.priority_synergies[synergy_pair]['description']
        
        # General explanation for multi-category synergies
        return f"Multi-pathway synergy detected involving {len(synergy_categories)} functional categories: {', '.join(synergy_categories)}"