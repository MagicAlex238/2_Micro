import re
from collections import OrderedDict, defaultdict
from typing import Dict, List, Set, Optional

'''class TermProcessor:
    """Normalizes and matches corrosion terms with priority handling"""
    
    def __init__(self, taxonomy: dict):# to garanty priority found as subcategories on fc.global terms
        self.priority_order = ['iron_metabolism', 'sulfur_metabolism', 'organic_acid_metabolism',
        'biofilm_formation', 'o2_consumption', 'metal_binding_chelation',
        'nitrogen_metabolism', 'manganese_processes', 'h2_consumption',
        'methanogenesis', 'carbon_metabolism',  'fumarate_formation',
        'phosphorus_metabolism']
        
        # This map will link a specific keyword back to its category. e.g., {'iron': 'iron_metabolism', 'sulfide': 'sulfur_metabolism'}
        self.keyword_to_category_map = self._build_keyword_map(taxonomy)
        
        # This pre-compiles one single, highly efficient regex pattern from all keywords.
        self.master_regex = self._compile_master_pattern()
    
    def _build_keyword_map(self, taxonomy: dict) -> dict:
        """Creates a lookup map from a term to its parent category."""
        keyword_map = {}
        for category, data in taxonomy.items():
            # Handles both dictionary and list formats for terms
            terms = data.get('terms', data) if isinstance(data, dict) else data
            for term in terms:
                # We normalize the term here for robust matching
                keyword_map[self._normalize_term(term)] = category
        return keyword_map
        
    def _compile_master_pattern(self):
        """Compiles all keywords into a single, case-insensitive regex pattern."""
        if not self.keyword_to_category_map:
            return None
        # The `\b` ensures we match whole words only.
        pattern = r'\b(' + '|'.join(re.escape(k) for k in self.keyword_to_category_map.keys()) + r')\b'
        return re.compile(pattern, re.IGNORECASE)

    def find_all_matches(self, text: str) -> dict:
        """
        Finds all keywords in a text and groups them by their functional category.
        Returns a dictionary like: {'iron_metabolism': ['iron'], 'sulfur_metabolism': ['sulfide']}
        """
        if self.master_regex is None or not isinstance(text, str):
            return {}
            
        # Normalize the input text once
        normalized_text = self._normalize_term(text)
        
        # Find all unique keywords in the text
        found_keywords = set(self.master_regex.findall(normalized_text))
        
        # Group the found keywords by their category using the map
        categorized_matches = defaultdict(list)
        for keyword in found_keywords:
            category = self.keyword_to_category_map.get(keyword)
            if category:
                categorized_matches[category].append(keyword)
                
        return dict(categorized_matches)

    def _normalize_term(self, term: str) -> str:
        """Enhanced normalization with corrosion-specific substitutions to catch missmatches"""
        substitutions = {
            'reduc': 'reduction',
            'oxid': 'oxidation', 
            'ferri': 'iron',
            'ferro': 'iron',
            'sulph': 'sulf',
            'metallo': 'metal',
            'corrosi': 'corrosion'
        }
        
        base_term = re.sub(r'[\W_]+', '', term).lower()
        for pattern, replacement in substitutions.items():
            base_term = re.sub(f'^{pattern}', replacement, base_term)
        return base_term

    def matches_normalized(self, term: str, text: str) -> bool:
        """Check if normalized term appears in text"""
        norm_term = self._normalize_term(term)
        norm_text = self._normalize_term(text)
        return norm_term in norm_text
    
    def _create_priority_dict(self, taxonomy: dict) -> OrderedDict:
        """Create search priority structure"""
        priority_dict = OrderedDict()
        
        for cat in self.priority_order:
            if cat in taxonomy:
                if isinstance(taxonomy[cat], dict) and 'terms' in taxonomy[cat]:
                    # Handle functional_categories format
                    priority_dict[cat] = {self._normalize_term(t) for t in taxonomy[cat]['terms']}
                elif isinstance(taxonomy[cat], list):
                    # Handle simple list format
                    priority_dict[cat] = {self._normalize_term(t) for t in taxonomy[cat]}
        
        # Add any remaining categories not in priority order
        for cat, terms in taxonomy.items():
            if cat not in priority_dict:
                if isinstance(terms, dict) and 'terms' in terms:
                    priority_dict[cat] = {self._normalize_term(t) for t in terms['terms']}
                elif isinstance(terms, list):
                    priority_dict[cat] = {self._normalize_term(t) for t in terms}
        
        return priority_dict
    
    def find_first_category(self, term: str) -> str:
        """Returns first matching category based on priority"""
        norm_term = self._normalize_term(term)
        for category in self.priority_order:
            if norm_term in [self._normalize_term(t) for t in self.keyword_to_category_map.get(category, {}).get('terms', [])]:
                return category

        # Fallback for non-priority categories
        for keyword, category in self.keyword_to_category_map.items():
            if norm_term == keyword:
                return category
        return None
'''
class TermProcessor:
    """Normalizes and matches corrosion terms with priority handling."""

    def __init__(self, taxonomy: dict):
        # Priority includes all categories you referenced elsewhere (synergy and FC lists)
        self.priority_order: List[str] = [
            'iron_metabolism', 'sulfur_metabolism', 'organic_acid_metabolism',
            'biofilm_formation', 'o2_consumption', 'metal_binding_chelation',
            'nitrogen_metabolism', 'manganese_processes', 'h2_consumption',
            'methanogenesis', 'carbon_metabolism', 'fumarate_formation',
            'phosphorus_metabolism', 'direct_eet', 'indirect_eet', 'halogen_related'
        ]

        # Build category -> set(normalized terms) preserving priority
        self.category_to_terms: Dict[str, Set[str]] = self._build_category_to_terms(taxonomy)

        # Build keyword (normalized term) -> category, preferring priority categories
        self.keyword_to_category_map: Dict[str, str] = self._build_keyword_map_with_priority()

        # Precompile a regex over normalized keywords; longest first to prefer longer phrases
        self.master_regex = self._compile_master_pattern()

    def _build_category_to_terms(self, taxonomy: dict) -> Dict[str, Set[str]]:
        """Create category -> set of normalized terms, honoring priority_order."""
        cat_to_terms: "OrderedDict[str, Set[str]]" = OrderedDict()

        # Start with priority categories in order
        for cat in self.priority_order:
            if cat in taxonomy:
                terms = taxonomy[cat]['terms'] if isinstance(taxonomy[cat], dict) and 'terms' in taxonomy[cat] else taxonomy[cat]
                cat_to_terms[cat] = {self._normalize_term(t) for t in terms}

        # Add remaining categories not in priority list
        for cat, data in taxonomy.items():
            if cat not in cat_to_terms:
                terms = data['terms'] if isinstance(data, dict) and 'terms' in data else data
                cat_to_terms[cat] = {self._normalize_term(t) for t in terms}

        return cat_to_terms

    def _build_keyword_map_with_priority(self) -> Dict[str, str]:
        """Reverse map terms -> category, giving precedence to earlier (priority) categories."""
        kw_to_cat: Dict[str, str] = {}
        for cat, terms in self.category_to_terms.items():
            for term in terms:
                # First category wins; this respects priority ordering already applied
                kw_to_cat.setdefault(term, cat)
        return kw_to_cat

    def _compile_master_pattern(self) -> Optional[re.Pattern]:
        """Compile a single, case-insensitive regex for all normalized keywords."""
        if not self.keyword_to_category_map:
            return None

        # Sort by length desc to reduce premature short-term matches
        kws = sorted(self.keyword_to_category_map.keys(), key=len, reverse=True)

        # NOTE: Text is normalized to lowercase with spaces; keywords are too.
        # We can safely use word boundaries.
        pattern = r'\b(' + '|'.join(re.escape(k) for k in kws) + r')\b'
        return re.compile(pattern, re.IGNORECASE)

    def _normalize_term(self, term: str) -> str:
        """Normalization that keeps token boundaries, aiding regex word-boundary matching."""
        if not isinstance(term, str):
            return ""

        # Lowercase
        t = term.lower()

        # Corrosion-specific substitutions (apply anywhere, not only at start)
        substitutions = {
            'reduc': 'reduction',
            'oxid': 'oxidation',
            'ferri': 'iron',
            'ferro': 'iron',
            'sulph': 'sulf',
            'metallo': 'metal',
            'corrosi': 'corrosion',
        }
        for pattern, replacement in substitutions.items():
            t = re.sub(pattern, replacement, t)

        # Replace non-word/underscore runs with a single space to preserve boundaries
        t = re.sub(r'[\W_]+', ' ', t)

        # Collapse multiple spaces and strip
        t = re.sub(r'\s+', ' ', t).strip()

        return t

    def find_all_matches(self, text: str) -> dict:
        """
        Finds all keywords in a text and groups them by their functional category.
        Returns: {'iron_metabolism': ['iron'], 'sulfur_metabolism': ['sulfide'], ...}
        """
        if self.master_regex is None or not isinstance(text, str):
            return {}

        normalized_text = self._normalize_term(text)
        if not normalized_text:
            return {}

        found_keywords = set(self.master_regex.findall(normalized_text))  # already normalized words/phrases

        categorized_matches = defaultdict(list)
        for keyword in found_keywords:
            # Normalize match to be safe, though it should already be normalized casing-wise
            k = self._normalize_term(keyword)
            cat = self.keyword_to_category_map.get(k)
            if cat:
                categorized_matches[cat].append(k)

        return dict(categorized_matches)

    def matches_normalized(self, term: str, text: str) -> bool:
        """Check if the normalized term substring appears in the normalized text."""
        norm_term = self._normalize_term(term)
        if not norm_term:
            return False
        norm_text = self._normalize_term(text)
        return norm_term in norm_text

    def find_first_category(self, term: str) -> Optional[str]:
        """Return the category for a single token/term, honoring priority where ambiguous."""
        norm_term = self._normalize_term(term)
        if not norm_term:
            return None

        # Direct lookup (fast path)
        cat = self.keyword_to_category_map.get(norm_term)
        if cat:
            return cat

        # Fallback: check membership across categories in priority order
        for category in self.priority_order:
            if norm_term in self.category_to_terms.get(category, set()):
                return category

        # Last resort: scan remaining categories
        for category, terms in self.category_to_terms.items():
            if category in self.priority_order:
                continue
            if norm_term in terms:
                return category

        return None