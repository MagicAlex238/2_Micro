#!/usr/bin/env python
# coding: utf-8

"""
Corrosion Relevance Scoring System: functions for corrosion relevance evaluation.
Updated to use unified functional_categories scoring approach.
"""
import math
import sys
import os

try:
    # Try relative import (for package installation)
    from .global_terms import (
        metal_terms,
        corrosion_synergies,
        functional_categories,
        metal_mapping,
    )
except ImportError:
    print("Critical error")

# Scoring weights - keeping existing structure
METAL_SCORE_WEIGHT = 0.5 
FUNCTIONAL_SCORE_WEIGHT = 1.5
SYNERGY_SCORE_WEIGHT = 2.0 

# Classification thresholds
HIGH_RELEVANCE_THRESHOLD = 5.0
MEDIUM_RELEVANCE_THRESHOLD = 2.0

def consolidate_metal_terms(brenda_metals, text_detected_metals):
    """Consolidate metal names from BRENDA and text mining into standardized symbols.
    Args:
        brenda_metals: Metals obtained from BRENDA data
        text_detected_metals: Metals detected from text mining
    Returns:
        Consolidated list of unique, standardized metal symbols
    """
    consolidated = set()
    all_metals = (brenda_metals or []) + (text_detected_metals or [])

    for metal in all_metals:
        metal_norm = metal.strip().lower()
        # Check if the normalized term matches any key in the standard mapping
        for key, symbol in metal_mapping.items():
            if key in metal_norm:
                consolidated.add(symbol)
                break
        else:
            # If no mapping is found, add the original
            consolidated.add(metal.strip())
    return list(consolidated)

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

###=========================================================
"""
Corrosion Relevance Scoring System: functions for corrosion relevance evaluation.
Updated to use unified functional_categories scoring approach.
"""
import math
import sys
import os

try:
    # Try relative import (for package installation)
    from .global_terms import (
        metal_terms,
        corrosion_synergies,
        functional_categories,
        metal_mapping,
    )
except ImportError:
    print("Critical error")

# Scoring weights - updated structure
METAL_SCORE_WEIGHT = 0.5
FUNCTIONAL_SCORE_WEIGHT = 1.5
SYNERGY_SCORE_WEIGHT = 2.0

# Classification thresholds
HIGH_RELEVANCE_THRESHOLD = 5.0
MEDIUM_RELEVANCE_THRESHOLD = 2.0

def consolidate_metal_terms(brenda_metals, text_detected_metals):
    """Consolidate metal names from BRENDA and text mining into standardized symbols.
    
    Args:
        brenda_metals: Metals obtained from BRENDA data
        text_detected_metals: Metals detected from text mining
        
    Returns:
        Consolidated list of unique, standardized metal symbols
    """
    consolidated = set()
    all_metals = (brenda_metals or []) + (text_detected_metals or [])

    for metal in all_metals:
        metal_norm = metal.strip().lower()
        # Check if the normalized term matches any key in the standard mapping
        for key, symbol in metal_mapping.items():
            if key in metal_norm:
                consolidated.add(symbol)
                break
        else:
            # If no mapping is found, add the original
            consolidated.add(metal.strip())
    return list(consolidated)

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
    text_lower = text.lower()

    # Score metals - simplified since all data has metals
    metal_score = 0.0
    detected_metals = []
    for metal in metal_terms:
        if metal.lower() in text_lower:
            detected_metals.append(metal)
            metal_score += 1.0
    
    # Consolidate metals
    consolidated_metals = consolidate_metal_terms(brenda_metals, detected_metals)
    results["consolidated_metals"] = consolidated_metals
    results["overall_metal_score"] = float(metal_score)

    # Score synergies
    synergy_score = 0.0
    synergy_matches = {}
    for synergy_group, terms in corrosion_synergies.items():
        group_hits = 0
        for term in terms:
            if term.lower() in text_lower:
                group_hits += 1
        if group_hits > 0:
            synergy_matches[synergy_group] = math.log(group_hits + 1)
            synergy_score += synergy_matches[synergy_group]
    
    results["corrosion_synergies"] = list(synergy_matches.keys())
    results["corrosion_synergy_scores"] = synergy_matches
    results["overall_synergy_score"] = float(synergy_score)

    # Score functional categories - PRIMARY scoring method
    functional_score = 0.0
    func_matches = {}
    for cat, details in functional_categories.items():
        category_hits = 0
        for term in details["terms"]:
            if term.lower() in text_lower:
                category_hits += 1
        
        if category_hits > 0:
            base_score = math.log(category_hits + 1)
            weighted_score = base_score * details["score"]
            func_matches[cat] = weighted_score
            functional_score += weighted_score

    results["functional_categories"] = [
        {"category": cat, "score": score}
        for cat, score in func_matches.items()
    ]
    results["overall_functional_score"] = float(functional_score)

    return results

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
    # Apply weights - functional_score is now primary
    weighted_metal_score = metal_score * METAL_SCORE_WEIGHT
    weighted_synergy_score = synergy_score * SYNERGY_SCORE_WEIGHT
    weighted_functional_score = functional_score * FUNCTIONAL_SCORE_WEIGHT

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