# ========== __init__.py ==========
"""
Refactored Corrosion Scoring System v3.0

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


# ========== config.py ==========
"""
Configuration settings for the corrosion scoring system.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional

@dataclass
class ScoringConfig:
    """Configuration class for scoring parameters."""
    
    # Scoring weights
    metal_score_weight: float = 0.5
    functional_score_weight: float = 1.5
    synergy_score_weight: float = 2.0
    
    # Text processing limits
    max_text_length: int = 10000
    max_regex_keywords: int = 1000
    
    # Synergy detection parameters
    min_synergy_categories: int = 2
    min_synergy_terms: int = 2
    min_general_synergy_terms: int = 3
    
    # Priority categories for functional categories
    functional_priority_order: List[str] = None
    
    # Synergy priority mappings with scores
    priority_synergies: Dict[tuple, Dict[str, float]] = None
    
    def __post_init__(self):
        """Set default values for complex fields."""
        if self.functional_priority_order is None:
            self.functional_priority_order = [
                'iron_metabolism', 'sulfur_metabolism', 'organic_acid_metabolism',
                'biofilm_formation', 'o2_consumption', 'metal_binding_chelation',
                'nitrogen_metabolism', 'manganese_processes', 'h2_consumption',
                'methanogenesis', 'carbon_metabolism', 'fumarate_formation',
                'phosphorus_metabolism', 'direct_eet', 'indirect_eet', 'halogen_related'
            ]
        
        if self.priority_synergies is None:
            self.priority_synergies = {
                ('organic_acid_metabolism', 'metal_binding_chelation'): {
                    'score': 3.0, 
                    'description': 'TOC-chelation Synergy (TOC-chelate)'
                },
                ('iron_metabolism', 'organic_acid_metabolism'): {
                    'score': 2.8, 
                    'description': 'Iron-Organic Acid Synergy (acid-enhanced Fe corrosion)'
                },
                ('biofilm_formation', 'metal_binding_chelation'): {
                    'score': 2.7,
                    'description': 'biofilm-chelate Synergy (biofilm-chelate-corrosion)'
                },
                ('o2_consumption', 'iron_metabolism'): {
                    'score': 2.5,
                    'description': 'Oxygen-Iron Synergy (aerobic Fe corrosion)'
                },
                ('sulfur_metabolism', 'iron_metabolism'): {
                    'score': 2.4,
                    'description': 'Sulfur-iron Synergy (SRB-mediated corrosion)'
                },
                ('sulfur_metabolism', 'h2_consumption'): {
                    'score': 2.3,
                    'description': 'Sulfur-Hydrogen Synergy (SRB-mediated corrosion)'
                },
                ('biofilm_formation', 'iron_metabolism'): {
                    'score': 2.2,
                    'description': 'Biofilm-Iron Synergy (biofilm-enhanced Fe corrosion)'
                },
                ('nitrogen_metabolism', 'iron_metabolism'): {
                    'score': 2.0,
                    'description': 'Nitrogen-Iron Synergy (nitrate-enhanced Fe corrosion)'
                }
            }

# ========== exceptions.py ==========
"""
Custom exceptions for the corrosion scoring system.
"""

class ScoringError(Exception):
    """Base exception for scoring system errors."""
    pass

class TextMiningError(ScoringError):
    """Exception raised during text mining operations."""
    pass

class SynergyDetectionError(ScoringError):
    """Exception raised during synergy detection."""
    pass

class ProcessorError(ScoringError):
    """Exception raised when processors are invalid or missing."""
    pass

class ValidationError(ScoringError):
    """Exception raised when input validation fails."""
    pass

# ========== global_terms.py ========== too large


# ========== score_calculator.py ==========
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
            from ..corrosion_scoring_v3.global_terms import (
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

# ========== synergy_detector.py =========
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

# ========== term_processor.py ==========
"""
Fixed TermProcessor class with proper error handling and validation.
"""

import re
from collections import OrderedDict, defaultdict
from typing import Dict, List, Set, Optional, Pattern
from ..corrosion_scoring_v3.exceptions import ProcessorError

class TermProcessor:
    """
    Normalizes and matches corrosion terms with priority handling.
    Fixed version with proper error handling and validation.
    """

    def __init__(self, taxonomy: dict, synonyms: dict = None):
        """
        Initialize the TermProcessor.
        
        Args:
            taxonomy: Dictionary containing term taxonomy
            synonyms: Optional dictionary of synonyms
            
        Raises:
            ProcessorError: If initialization fails
        """
        if not isinstance(taxonomy, dict):
            raise ProcessorError("Taxonomy must be a dictionary")
        
        self.synonyms = synonyms or {
            "oxidoreductase": ["oxid reduction enzyme"],
            "iron sulfur": ["fe s", "fe s cluster", "iron-sulfur", "fe-s"],
            "sulfide": ["sulphide"],
            "sulfite": ["sulphite"],
            "sulfur": ["sulphur"],
        }
        
        # Priority order for resolving conflicts
        self.priority_order: List[str] = [
            'iron_metabolism', 'sulfur_metabolism', 'organic_acid_metabolism',
            'biofilm_formation', 'o2_consumption', 'metal_binding_chelation',
            'nitrogen_metabolism', 'manganese_processes', 'h2_consumption',
            'methanogenesis', 'carbon_metabolism', 'fumarate_formation',
            'phosphorus_metabolism', 'direct_eet', 'indirect_eet', 'halogen_related'
        ]

        try:
            # Build category -> set(normalized terms) preserving priority
            self.category_to_terms: Dict[str, Set[str]] = self._build_category_to_terms(taxonomy)

            # Build keyword (normalized term) -> category, preferring priority categories
            self.keyword_to_category_map: Dict[str, str] = self._build_keyword_map_with_priority()

            # Precompile a regex over normalized keywords
            self.master_regex: Optional[Pattern] = self._compile_master_pattern()
            
        except Exception as e:
            raise ProcessorError(f"Failed to initialize TermProcessor: {e}") from e

    def _build_category_to_terms(self, taxonomy: dict) -> Dict[str, Set[str]]:
        """Create category -> set of normalized terms, honoring priority_order."""
        cat_to_terms: OrderedDict[str, Set[str]] = OrderedDict()

        # Start with priority categories in order
        for cat in self.priority_order:
            if cat in taxonomy:
                try:
                    terms = self._extract_terms_from_category(taxonomy[cat])
                    expanded_terms = self._expand_terms_with_synonyms(terms)
                    cat_to_terms[cat] = {self._normalize_term(t) for t in expanded_terms if self._normalize_term(t)}
                except Exception as e:
                    raise ProcessorError(f"Error processing category '{cat}': {e}") from e

        # Add remaining categories not in priority list
        for cat, data in taxonomy.items():
            if cat not in cat_to_terms:
                try:
                    terms = self._extract_terms_from_category(data)
                    expanded_terms = self._expand_terms_with_synonyms(terms)
                    cat_to_terms[cat] = {self._normalize_term(t) for t in expanded_terms if self._normalize_term(t)}
                except Exception as e:
                    raise ProcessorError(f"Error processing category '{cat}': {e}") from e

        return cat_to_terms
    
    def _extract_terms_from_category(self, data) -> List[str]:
        """Extract terms from category data, handling different formats."""
        if isinstance(data, dict) and 'terms' in data:
            terms = data['terms']
        else:
            terms = data
        
        if not isinstance(terms, (list, tuple)):
            raise ProcessorError(f"Terms must be a list or tuple, got {type(terms)}")
        
        return list(terms)
    
    def _expand_terms_with_synonyms(self, terms: List[str]) -> List[str]:
        """Expand terms with synonyms."""
        expanded = []
        for term in terms:
            expanded.append(term)
            base = self._normalize_term(term)
            if base in self.synonyms:
                expanded.extend(self.synonyms[base])
        return expanded

    def _build_keyword_map_with_priority(self) -> Dict[str, str]:
        """Reverse map terms -> category, giving precedence to earlier (priority) categories."""
        kw_to_cat: Dict[str, str] = {}
        for cat, terms in self.category_to_terms.items():
            for term in terms:
                # First category wins; this respects priority ordering
                if term not in kw_to_cat:
                    kw_to_cat[term] = cat
        return kw_to_cat

    def _compile_master_pattern(self) -> Optional[Pattern]:
        """Compile a single, case-insensitive regex for all normalized keywords."""
        if not self.keyword_to_category_map:
            return None

        # Sort by length desc to prefer longer matches
        kws = sorted(self.keyword_to_category_map.keys(), key=len, reverse=True)
        
        # Limit pattern size to avoid catastrophic backtracking
        max_keywords = 1000
        if len(kws) > max_keywords:
            kws = kws[:max_keywords]

        try:
            # Use word boundaries for better matching
            pattern = r'\b(' + '|'.join(re.escape(k) for k in kws) + r')\b'
            return re.compile(pattern, re.IGNORECASE)
        except re.error as e:
            raise ProcessorError(f"Failed to compile regex pattern: {e}") from e

    def _normalize_term(self, term: str) -> str:
        """Core normalizer used internally by the processor."""
        if not isinstance(term, str):
            return ""

        t = term.lower()

        # Conservative substitutions (avoid corrupting 'oxidase', 'reductase', etc.)
        t = re.sub(r'\bsulphur\b', 'sulfur', t)
        t = re.sub(r'\bsulphate\b', 'sulfate', t)
        t = re.sub(r'\bsulphide\b', 'sulfide', t)
        t = re.sub(r'\bmetallo(?=enzyme|protein)\b', 'metal', t)
        t = re.sub(r'\bcorrosion\w*\b', 'corrosion', t)

        # Replace non-word/underscore runs with a single space
        t = re.sub(r'[\W_]+', ' ', t)

        # Collapse multiple spaces and trim
        t = re.sub(r'\s+', ' ', t).strip()

        return t

    def normalize_text(self, text: str) -> str:
        """Public wrapper for external callers."""
        return self._normalize_term(text)

    def find_all_matches(self, text: str) -> Dict[str, List[str]]:
        """
        Find all keywords in text and group them by functional category.
        
        Args:
            text: Text to search
            
        Returns:
            Dictionary mapping categories to lists of found terms
            
        Raises:
            ProcessorError: If matching fails
        """
        if self.master_regex is None:
            return {}
        
        if not isinstance(text, str):
            return {}

        try:
            normalized_text = self._normalize_term(text)
            if not normalized_text:
                return {}

            found_keywords = set(self.master_regex.findall(normalized_text))

            categorized_matches = defaultdict(list)
            for keyword in found_keywords:
                normalized_keyword = self._normalize_term(keyword)
                category = self.keyword_to_category_map.get(normalized_keyword)
                if category:
                    categorized_matches[category].append(normalized_keyword)

            return dict(categorized_matches)
            
        except Exception as e:
            raise ProcessorError(f"Error finding matches in text: {e}") from e

    def matches_normalized(self, term: str, text: str) -> bool:
        """Check if the normalized term appears in the normalized text."""
        try:
            norm_term = self._normalize_term(term)
            if not norm_term:
                return False
            norm_text = self._normalize_term(text)
            return bool(re.search(rf'\b{re.escape(norm_term)}\b', norm_text))
        except Exception:
            return False

    def find_first_category(self, term: str) -> Optional[str]:
        """Return the category for a single token/term, honoring priority."""
        try:
            norm_term = self._normalize_term(term)
            if not norm_term:
                return None

            # Direct lookup (fast path)
            category = self.keyword_to_category_map.get(norm_term)
            if category:
                return category

            # Fallback: check membership across categories in priority order
            for category in self.priority_order:
                if norm_term in self.category_to_terms.get(category, set()):
                    return category

            # Last resort: scan remaining categories
            for category, terms in self.category_to_terms.items():
                if category not in self.priority_order and norm_term in terms:
                    return category

            return None
            
        except Exception:
            return None

# ========== text_miner.py ==========
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
            from ..corrosion_scoring_v3.global_terms import (
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
            from ..corrosion_scoring_v3.term_processor import TermProcessor
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
            from ..corrosion_scoring_v3.term_processor import TermProcessor
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
            from ..corrosion_scoring_v3.term_processor import TermProcessor
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

# ========== validators.py ==========
"""
Input validation utilities for the scoring system.
"""

import re
from typing import Any, Dict, List, Optional, Union
from .exceptions import ValidationError

def validate_text(text: Any, max_length: int = 10000, allow_empty: bool = False) -> str:
    """
    Validate and sanitize input text.
    
    Args:
        text: Input text to validate
        max_length: Maximum allowed text length
        allow_empty: Whether to allow empty text
        
    Returns:
        Sanitized text string
        
    Raises:
        ValidationError: If text is invalid
    """
    if text is None:
        if allow_empty:
            return ""
        raise ValidationError("Text cannot be None")
    
    if not isinstance(text, str):
        try:
            text = str(text)
        except Exception as e:
            raise ValidationError(f"Cannot convert text to string: {e}") from e
    
    if not allow_empty and not text.strip():
        raise ValidationError("Text cannot be empty")
    
    if len(text) > max_length:
        text = text[:max_length]
    
    # Sanitize text - remove potentially problematic characters for regex
    text = re.sub(r'[\\^$.*+?{}[\]|()\x00-\x1f]', ' ', text)
    text = re.sub(r'\s+', ' ', text).strip()
    
    return text

def validate_processor(processor: Any, processor_name: str) -> None:
    """
    Validate that a processor has required methods.
    
    Args:
        processor: Processor object to validate
        processor_name: Name of the processor for error messages
        
    Raises:
        ValidationError: If processor is invalid
    """
    if processor is None:
        raise ValidationError(f"{processor_name} processor cannot be None")
    
    required_methods = ['find_all_matches', 'matches_normalized']
    for method in required_methods:
        if not hasattr(processor, method):
            raise ValidationError(f"{processor_name} processor missing required method: {method}")

def validate_processors(processors: Dict[str, Any]) -> None:
    """
    Validate all required processors are present and valid.
    
    Args:
        processors: Dictionary of processors
        
    Raises:
        ValidationError: If processors are invalid
    """
    required_processors = ['fc_processor', 'metal_processor', 'synergy_processor']
    
    if not isinstance(processors, dict):
        raise ValidationError("Processors must be provided as a dictionary")
    
    for proc_name in required_processors:
        if proc_name not in processors:
            raise ValidationError(f"Missing required processor: {proc_name}")
        validate_processor(processors[proc_name], proc_name)

def validate_metals_list(metals: Any) -> List[str]:
    """
    Validate and clean metals list.
    
    Args:
        metals: Input metals list
        
    Returns:
        Cleaned list of metal strings
    """
    if metals is None:
        return []
    
    if not isinstance(metals, (list, tuple)):
        metals = [metals]
    
    cleaned_metals = []
    for metal in metals:
        if metal is not None:
            metal_str = str(metal).strip()
            if metal_str and metal_str.lower() not in ['not detected', 'none', '']:
                cleaned_metals.append(metal_str)
    
    return cleaned_metals