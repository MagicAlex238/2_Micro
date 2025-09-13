import re
from collections import OrderedDict, defaultdict
from typing import Dict, List, Set, Optional
class TermProcessor:
    """Normalizes and matches corrosion terms with priority handling."""

    def __init__(self, taxonomy: dict, synonyms: dict = None):
        self.synonyms = synonyms or {
            "oxidoreductase": ["oxid reduction enzyme"],
            "iron sulfur": ["fe s", "fe s cluster", "iron-sulfur", "fe-s"],
            "sulfide": ["sulphide"],
            "sulfite": ["sulphite"],
            "sulfur": ["sulphur"],
        }
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
                expanded = []
                for t in terms:
                    expanded.append(t)
                    base = self._normalize_term(t)
                    for syn in self.synonyms.get(base, []):
                        expanded.append(syn)
                cat_to_terms[cat] = {nt for t in expanded if (nt := self._normalize_term(t))}

        # Add remaining categories not in priority list
        for cat, data in taxonomy.items():
            if cat not in cat_to_terms:
                terms = data['terms'] if isinstance(data, dict) and 'terms' in data else data
                cat_to_terms[cat] = {nt for t in terms if (nt := self._normalize_term(t))}

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
        # limit pattern size to avoid backtracking
        if len(kws) >1000:
            kws = kws[:1000]
        # NOTE: Text is normalized to lowercase with spaces; keywords are too. safely use word boundaries.
        pattern = r'\b(' + '|'.join(re.escape(k) for k in kws) + r')\b'
        return re.compile(pattern, re.IGNORECASE)

    def _normalize_term(self, term: str) -> str:
        """Core normalizer used internally by the processor."""
        if not isinstance(term, str):
            return ""

        t = term.lower()

        # Conservative substitutions (avoid corrupting 'oxidase', 'reductase', etc.)
        t = re.sub(r'\bsulphur\b', 'sulfur', t)     # Complete word only
        t = re.sub(r'\bsulphate\b', 'sulfate', t)   # Complete word only  
        t = re.sub(r'\bsulphide\b', 'sulfide', t)   # Complete word only
        t = re.sub(r'\bmetallo(?=enzyme|protein)\b', 'metal', t)  # Only before enzyme/protein
        t = re.sub(r'\bcorrosion\w*\b', 'corrosion', t)  # Normalize corrosion variants

        # Replace non-word/underscore runs with a single space (keeps alnum boundaries)
        t = re.sub(r'[\W_]+', ' ', t)

        # Collapse multiple spaces and trim
        t = re.sub(r'\s+', ' ', t).strip()

        return t

    def normalize_text(self, text: str) -> str:
        """Public wrapper for external callers; do not call the private method outside."""
        return self._normalize_term(text)

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
        return re.search(rf'\b{re.escape(norm_term)}\b', norm_text) is not None

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