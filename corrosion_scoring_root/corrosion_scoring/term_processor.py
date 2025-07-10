import re
from collections import OrderedDict

class TermProcessor:
    """Normalizes and matches corrosion terms with priority handling"""
    
    def __init__(self, taxonomy: dict):
        self.priority_order = [
            'iron', 'nickel', 'Mo', 'V5+', 'Cr3+',  # Metal priorities
            'chloride', 'sulfide', 'sulfate',  # Anion priorities
            'corrosion_mechanisms', 'pathway_categories'  # Process priorities
        ]
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
        return OrderedDict(
            (cat, {self._normalize_term(t) for t in terms})
            for cat in self.priority_order if cat in taxonomy
        )

    def find_first_category(self, term: str) -> str:
        """Returns first matching category based on priority"""
        norm_term = self._normalize_term(term)
        for category, terms in self.normalized_taxonomy.items():
            if norm_term in terms:
                return category
        return None
