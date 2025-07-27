

def calculate_overall_scores(text, brenda_metals=None, pathways=None):
    """Calculate all the overall scores for a given text.
    Args:
        text: Text to analyze (combined enzyme names, class, pathways, reactions)
        brenda_metals: Metals from BRENDA database
        pathways: Pathway information for additional mechanism detection
    Returns:
        Dictionary containing all overall scores and matched categories
    """
    if brenda_metals is None:
        brenda_metals = []

    results = {}

    # Score metals
    metal_score, metal_matches = score_keyword_matches(text, metal_terms)
    for metal in brenda_metals:
        if metal not in metal_matches:
            metal_matches[metal] = 1.0

    results["metals_involved"] = list(metal_matches.keys())
    results["metal_scores"] = metal_matches
    results["overall_metal_score"] = float(metal_score)

    # Score corrosion mechanisms
    corrosion_score, corrosion_matches = score_keyword_matches(text, corrosion_mechanisms)   
    if pathways is None:
        pathways = []
    if pathways:
        pathway_mechanisms = assign_mechanism_from_pathway(pathways)
        for mechanism in pathway_mechanisms:
            if mechanism not in corrosion_matches:
                corrosion_matches[mechanism] = 1.0  

    results["overall_corrosion_score"] = float(corrosion_score)

    # Score synergies
    synergy_score, synergy_matches = score_keyword_matches(text, corrosion_synergies)
    results["corrosion_synergies"] = list(synergy_matches.keys())
    results["corrosion_synergy_scores"] = synergy_matches
    results["overall_synergy_score"] = float(synergy_score)

    # Score functional categories - ENHANCED to be primary scoring method
    functional_terms = {
        cat: details["terms"] for cat, details in functional_categories.items()
    }
    func_score, func_matches = score_keyword_matches(text, functional_terms)
    weighted_func_matches = {}
    for cat, match_score in func_matches.items():
        original_weight = functional_categories[cat]["score"]
        weighted_func_matches[cat] = match_score * original_weight

    results["functional_categories"] = [
        {"category": cat, "score": score}
        for cat, score in weighted_func_matches.items()
    ]
    results["overall_functional_score"] = float(sum(weighted_func_matches.values()))

    return results

def check_metal_organic_synergy(metals, enzyme_names, pathways):
    """Check for synergistic effects between metals and organic compounds.
    Args:
        metals: List of metals
        enzyme_names: List of enzyme names
        pathways: List of pathways
    Returns:
        Synergy score
    """
    name_text = " ".join(enzyme_names or []).lower()
    pathway_text = " ".join(pathways or []).lower()

    # Organic acid terms from pathway_categories and organic_categories
    organic_terms = []
    for category in cs.functional_categories

    if any(
        metal in metals for metal in cs.metal_terms):
        return 1.0  

    return 0.0


def calculate_corrosion_relevance_score(
    metal_score,
    synergy_score=0,
    functional_score=0,
):
    """Calculate final corrosion relevance score and category.
    UPDATED: functional_score now carries more weight as primary scoring method
    Args:
        metal_score: Metal involvement score
        synergy_score: Synergy score
        functional_score: Functional category score (PRIMARY)
    Returns:
        Corrosion relevance score and category
    """
    # Apply weights - UPDATED to emphasize functional_score
    weighted_metal_score = metal_score * METAL_SCORE_WEIGHT
    weighted_synergy_score = synergy_score * SYNERGY_SCORE_WEIGHT
    weighted_functional_score = functional_score * (FUNCTIONAL_SCORE_WEIGHT * 1.5)  # Increased weight as primary method

    # Calculate final score
    corrosion_relevance_score = float(
        weighted_metal_score
        + weighted_synergy_score
        + weighted_functional_score
    )

    # Determine category
    if corrosion_relevance_score >= HIGH_RELEVANCE_THRESHOLD:
        corrosion_relevance = "high"
    elif corrosion_relevance_score >= MEDIUM_RELEVANCE_THRESHOLD:
        corrosion_relevance = "medium"
    else:
        corrosion_relevance = "low"

    return corrosion_relevance_score, corrosion_relevance


#=================================================================================================================
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
