import re
from collections import OrderedDict

class TermProcessor:
    """Normalizes and matches corrosion terms with priority handling"""
    
    def __init__(self, taxonomy: dict):
        self.priority_order = ['iron_metabolism', 'sulfur_metabolism', 'organic_acid_metabolism',
                                 'direct_eet', 'biofilm_formation', 'o2_consumption', 'metal binding / chelation', 'nitrogen_metabolism', 'manganese_processes', 'h2_consumption', 'halogen_related', 'methanogenesis', 'carbon_metabolism', 'indirect_eet', 'fumarate_formation', 'phosphorus_metabolism']
    self.normalized_taxonomy = self._create_priority_dict(taxonomy)

    def _normalize_term(self, term: str) -> str:
        """Enhanced normalization with corrosion-specific substitutions"""
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
        for category, terms in self.normalized_taxonomy.items():
            if norm_term in terms:
                return category
        return None