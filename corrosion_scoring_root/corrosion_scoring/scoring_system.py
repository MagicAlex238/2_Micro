#!/usr/bin/env python
# coding: utf-8
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
        metal_mapping
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
    total_individual_synergy_terms_found = 0 

    for synergy_group, terms in corrosion_synergies.items():
        group_hits = 0
        for term in terms:
            if term.lower() in text_lower:
                group_hits += 1
                total_individual_synergy_terms_found += 1

        if group_hits > 0:
            synergy_matches[synergy_group] = math.log(group_hits + 1)

    if total_individual_synergy_terms_found >= 2:
        for score_val in synergy_matches.values():
            synergy_score += score_val # Sum up scores if threshold met
        results["corrosion_synergies"] = list(synergy_matches.keys())
    else:
        # If less than 2 individual terms found, reset score and matches
        synergy_score = 0.0
        synergy_matches = {}
        results["corrosion_synergies"] = []
    # contextual matching of synergies across different categories0

    inferred_synergies_found = []
    inferred_synergy_score = 0.0

    # Rule 1: Metal + Specific Functional Category + Specific Mechanism
    if 'Fe' in consolidated_metals and \
    'Oxidoreductase Activity' in functional_categories and \
    'Microbial Influenced Corrosion' in mechanisms_found:
    inferred_synergies_found.append('Fe-MIC-Redox Synergy')
    inferred_synergy_score += 5.0 # Example score

    # Rule 2: Presence of 'sulfide' compound + Metal + Reductase FC
    if any('sulfide' in c for c in corrosion_relevant_compounds) and \
    ('Fe' in consolidated_metals or 'Cu' in consolidated_metals) and \
    'Reductase Activity' in functional_categories: # Assuming 'Reductase Activity' is an FC
    inferred_synergies_found.append('Sulfidic Metal Biocorrosion')
    inferred_synergy_score += 7.0

    # Add these to your results dictionary:
    results['inferred_synergies'] = inferred_synergies_found
    results['overall_inferred_synergy_score'] = inferred_synergy_score
    
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