

from taxonomy import metal_terms, corrosion_mechanisms
from term_processor import TermProcessor

class CorrosionScorer:
    def __init__(self):
        self.term_processor = TermProcessor({
            **metal_terms,
            **corrosion_mechanisms
        })
    
    def score_abundances(self, abundance_data: dict) -> dict:
        """Process abundance data with priority-based scoring"""
        scores = defaultdict(float)
        for term, value in abundance_data.items():
            category = self.term_processor.find_first_category(term)
            if category:
                scores[category] += value
        return scores
