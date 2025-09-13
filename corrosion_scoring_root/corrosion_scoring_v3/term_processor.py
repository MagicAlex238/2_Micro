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