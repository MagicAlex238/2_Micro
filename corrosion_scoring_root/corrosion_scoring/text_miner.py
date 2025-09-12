"""
TextMiner class for extracting features from text without scoring.
"""

import math
import re
from typing import Dict, List, Set, Tuple, Any, Optional
from collections import defaultdict

from .validators import validate_text, validate_processors, validate_metals_list
from .exceptions import TextMiningError
from .config import ScoringConfig

class TextMiner:
    """
    Handles text mining and feature extraction for corrosion analysis.
    
    This class is responsible for extracting relevant terms and categories
    from text using various processors, but does NOT perform scoring.
    """
    
    def __init__(self, processors: Dict[str, Any], config: ScoringConfig = None):
        """
        Initialize the TextMiner.
        
        Args:
            processors: Dictionary containing required processors
            config: Configuration settings
            
        Raises:
            TextMiningError: If processors are invalid
        """
        try:
            validate_processors(processors)
        except Exception as e:
            raise TextMiningError(f"Invalid processors: {e}") from e
        
        self.processors = processors
        self.config = config or ScoringConfig()
        
        # Import global dictionaries safely
        try:
            from .global_terms import (
                metal_terms_dict,
                corrosion_synergies_dict,
                functional_categories_dict,
                metal_mapping,
                pathway_dict,
                mechanisms_dict,
                operational_environmental_factors_dict
            )
            self.metal_mapping = metal_mapping
            self.mechanisms_dict = mechanisms_dict
            self.pathway_dict = pathway_dict
            self.operational_dict = operational_environmental_factors_dict
        except ImportError as e:
            raise TextMiningError(f"Failed to import global terms: {e}") from e
    
    def extract_all_features(self, text: str, brenda_metals: List[str] = None) -> Dict[str, Any]:
        """
        Extract all features from text.
        
        Args:
            text: Text to analyze
            brenda_metals: Optional list of metals from BRENDA
            
        Returns:
            Dictionary containing extracted features
            
        Raises:
            TextMiningError: If feature extraction fails
        """
        try:
            # Validate and sanitize input
            clean_text = validate_text(text, self.config.max_text_length)
            clean_metals = validate_metals_list(brenda_metals)
            
            features = {}
            
            # Extract functional categories
            fc_features = self._extract_functional_categories(clean_text)
            features.update(fc_features)
            
            # Extract metals
            metal_features = self._extract_metals(clean_text, clean_metals)
            features.update(metal_features)
            
            # Extract mechanisms (for population only, not scoring)
            mechanism_features = self._extract_mechanisms(clean_text)
            features.update(mechanism_features)
            
            # Extract pathways (for population only, not scoring)
            pathway_features = self._extract_pathways(clean_text)
            features.update(pathway_features)
            
            # Extract operational factors (for population only, not scoring)
            operational_features = self._extract_operational_factors(clean_text)
            features.update(operational_features)
            
            return features
            
        except Exception as e:
            raise TextMiningError(f"Feature extraction failed: {e}") from e
    
    def _extract_functional_categories(self, text: str) -> Dict[str, Any]:
        """Extract functional category features."""
        try:
            fc_processor = self.processors['fc_processor']
            fc_matches = fc_processor.find_all_matches(text)
            
            fc_categories = []
            fc_child_terms = []
            fc_matches_detailed = {}
            
            for category, terms in fc_matches.items():
                if terms:
                    unique_terms = sorted(set(terms))
                    fc_categories.append(category)
                    fc_child_terms.extend(unique_terms)
                    fc_matches_detailed[category] = unique_terms
            
            return {
                'functional_categories_detected': fc_categories,
                'functional_child_terms': sorted(set(fc_child_terms)),
                'functional_matches_detailed': fc_matches_detailed
            }
            
        except Exception as e:
            raise TextMiningError(f"Functional category extraction failed: {e}") from e
    
    def _extract_metals(self, text: str, brenda_metals: List[str]) -> Dict[str, Any]:
        """Extract metal features."""
        try:
            metal_processor = self.processors['metal_processor']
            
            # Extract metals from text
            metal_matches = metal_processor.find_all_matches(text)
            detected_metal_categories = list(metal_matches.keys())
            
            # Consolidate metals from BRENDA and detected categories
            consolidated_metals = self._consolidate_metal_terms(brenda_metals, detected_metal_categories)
            
            return {
                'detected_metal_categories': detected_metal_categories,
                'brenda_metals': brenda_metals,
                'consolidated_metals': consolidated_metals,
                'metal_matches_detailed': {k: sorted(set(v)) for k, v in metal_matches.items()}
            }
            
        except Exception as e:
            raise TextMiningError(f"Metal extraction failed: {e}") from e
    
    def _extract_mechanisms(self, text: str) -> Dict[str, Any]:
        """Extract corrosion mechanisms for population (not scoring)."""
        try:
            from .term_processor import TermProcessor
            mechanisms_processor = TermProcessor(self.mechanisms_dict)
            
            mechanism_matches = mechanisms_processor.find_all_matches(text)
            detected_mechanisms = list(mechanism_matches.keys())
            
            return {
                'corrosion_mechanisms': detected_mechanisms,
                'mechanism_matches_detailed': {k: sorted(set(v)) for k, v in mechanism_matches.items()}
            }
            
        except Exception as e:
            raise TextMiningError(f"Mechanism extraction failed: {e}") from e
    
    def _extract_pathways(self, text: str) -> Dict[str, Any]:
        """Extract pathway information for population (not scoring)."""
        try:
            from .term_processor import TermProcessor
            pathway_processor = TermProcessor(self.pathway_dict)
            
            pathway_matches = pathway_processor.find_all_matches(text)
            detected_pathways = list(pathway_matches.keys())
            
            return {
                'detected_pathways': detected_pathways,
                'pathway_matches_detailed': {k: sorted(set(v)) for k, v in pathway_matches.items()}
            }
            
        except Exception as e:
            raise TextMiningError(f"Pathway extraction failed: {e}") from e
    
    def _extract_operational_factors(self, text: str) -> Dict[str, Any]:
        """Extract operational environmental factors for population (not scoring)."""
        try:
            from .term_processor import TermProcessor
            operational_processor = TermProcessor(self.operational_dict)
            
            operational_matches = operational_processor.find_all_matches(text)
            detected_operational = list(operational_matches.keys())
            
            return {
                'operational_environmental_factors': detected_operational,
                'operational_matches_detailed': {k: sorted(set(v)) for k, v in operational_matches.items()}
            }
            
        except Exception as e:
            raise TextMiningError(f"Operational factor extraction failed: {e}") from e
    
    def _consolidate_metal_terms(self, brenda_metals: List[str], detected_categories: List[str]) -> List[str]:
        """
        Consolidate metal names from BRENDA and detected categories into standardized symbols.
        
        Args:
            brenda_metals: Raw metal terms from BRENDA data
            detected_categories: Metal category names detected by TermProcessor
            
        Returns:
            List of consolidated, unique, standardized metal symbols
        """
        consolidated = set()
        
        # Process BRENDA metals (raw terms)
        for metal in brenda_metals:
            metal_raw = str(metal)
            metal_norm = metal_raw.strip().lower()
            
            if any(term in metal_norm for term in ['not detected', 'none']):
                continue
                
            # Strip brackets and the word ion
            metal_norm = re.sub(r'[\[\]\(\)]', '', metal_norm)
            metal_norm = re.sub(r'\bions?\b', '', metal_norm).strip()
            
            if 'not detected' in metal_norm:
                continue
                
            raw_tokens = re.findall(r'[a-z0-9\+\-]+', metal_norm)
            tokens = {re.sub(r'^(fe|cu|zn|ni|co|mn|cr|al|mg|ca|ba|sr|pb|as|hg)\d+\+?$', r'\1', t)
                     for t in raw_tokens}
            
            # Handle iron-sulfur cluster notation
            if re.search(r'\b\d+\s*fe\s*[-–]\s*\d+\s*s\b', metal_norm):
                tokens.add('fe')
            
            # Map tokens to standard symbols
            for token in tokens:
                # Check mapping first
                for key, symbol in self.metal_mapping.items():
                    if key.lower() == token:
                        consolidated.add(symbol)
                        break
                else:
                    # Check if it's a direct metal symbol
                    standard_metals = {
                        'fe', 'cu', 'ni', 'zn', 'co', 'mn', 'cr', 'al', 'mg', 'ca',
                        'ba', 'sr', 'pb', 'as', 'hg', 'mo', 'w', 'v', 'ti', 'sn',
                        'sb', 'cd', 'ag', 'au', 'pt'
                    }
                    if token in standard_metals:
                        consolidated.add(token.upper())
        
        # Process detected categories (category names like 'iron', 'copper')
        for category in detected_categories:
            if category in self.metal_mapping:
                consolidated.add(self.metal_mapping[category])
            else:
                consolidated.add(category)  # Fallback to category name
        
        return sorted(consolidated)