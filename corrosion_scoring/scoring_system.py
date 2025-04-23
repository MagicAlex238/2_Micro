#!/usr/bin/env python
# coding: utf-8

"""
Corrosion Relevance Scoring System: functions for corrosion relevance evaluation.
"""
import math
import sys
sys.path.append('/home/beatriz/MIC/2_Micro/corrosion_scoring') # /kaggle/input/corrosion-scoring
from corrosion_scoring.global_terms import (
    metal_terms, 
    corrosion_mechanisms, 
    pathway_categories, 
    organic_categories,
    corrosion_synergies, 
    functional_categories, 
    corrosion_keyword_groups, 
    metal_mapping
)

# Scoring weights
METAL_SCORE_WEIGHT = 1.5
CORROSION_MECHANISM_WEIGHT = 2.0
ORGANIC_PROCESS_WEIGHT = 1.0
KEYWORD_SCORE_WEIGHT = 0.5
FUNCTIONAL_SCORE_WEIGHT = 0.7
SYNERGY_SCORE_WEIGHT = 0.6

# Classification thresholds
HIGH_RELEVANCE_THRESHOLD = 5.0
MEDIUM_RELEVANCE_THRESHOLD = 2.0

def score_keyword_matches(text, term_dict):
    """Score text content against a dictionary of keywords with diminishing returns.
    Args:
        text: The text to analyze for keyword matches
        term_dict: Dictionary mapping categories to lists of related terms
    Returns:
        Overall score and dictionary of category scores
    """
    score = 0
    matches = {}
    
    if not text:
        return 0, {}
        
    text = text.lower()
    
    for category, terms in term_dict.items():
        category_matches = 0
        for term in terms:
            if term.lower() in text:
                category_matches += 1
        
        if category_matches > 0:
            # Score is proportional to match percentage but with diminishing returns
            match_ratio = category_matches / len(terms)
            category_score = math.sqrt(match_ratio) * 2  # Square root for diminishing returns
            score += category_score
            matches[category] = category_score
    
    return score, matches

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

def calculate_overall_scores(text, brenda_metals=None):
    """Calculate all the overall scores for a given text.
    Args:
        text: Text to analyze (combined enzyme names, class, pathways, reactions)
        brenda_metals: Metals from BRENDA database
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
    
    results['metals_involved'] = list(metal_matches.keys())
    results['metal_scores'] = metal_matches
    results['overall_metal_score'] = float(metal_score)
    
    # Score corrosion mechanisms
    corrosion_score, corrosion_matches = score_keyword_matches(text, corrosion_mechanisms)
    results['corrosion_mechanisms'] = list(corrosion_matches.keys())
    results['corrosion_mechanism_scores'] = corrosion_matches
    results['overall_corrosion_score'] = float(corrosion_score)
    
    # Score synergies
    synergy_score, synergy_matches = score_keyword_matches(text, corrosion_synergies)
    results['corrosion_synergies'] = list(synergy_matches.keys())
    results['corrosion_synergy_scores'] = synergy_matches
    results['overall_synergy_score'] = float(synergy_score)
    
    # Score functional categories
    functional_terms = {cat: details['terms'] for cat, details in functional_categories.items()}
    func_score, func_matches = score_keyword_matches(text, functional_terms)
    weighted_func_matches = {}
    for cat, match_score in func_matches.items():
        original_weight = functional_categories[cat]['score']
        weighted_func_matches[cat] = match_score * original_weight
    
    results['functional_categories'] = [{"category": cat, "score": score} 
                                     for cat, score in weighted_func_matches.items()]
    results['overall_functional_score'] = float(sum(weighted_func_matches.values()))
    
    # Score organic processes
    organic_score, organic_matches = score_keyword_matches(text, organic_categories)
    results['organic_processes'] = list(organic_matches.keys())
    results['organic_process_scores'] = organic_matches
    results['overall_organic_process_score'] = float(organic_score)
    
    # Score keyword groups
    keyword_score, keyword_matches = score_keyword_matches(text, corrosion_keyword_groups)
    results['corrosion_keyword_groups'] = list(keyword_matches.keys())
    results['corrosion_keyword_scores'] = keyword_matches
    results['overall_keyword_score'] = float(keyword_score)
    
    return results

def calculate_pathway_score(pathways, enzyme_names, enzyme_class):
    """Calculate pathway score based on corrosion relevance.
    Args:
        pathways: List of pathway names/descriptions
        enzyme_names: List of enzyme names
        enzyme_class: Enzyme class description
    Returns:
        Pathway score and dictionary of pathway category scores
    """
    pathway_score = 0
    pathway_category_scores = {}
    
    # Build text for analysis
    all_text = ' '.join(enzyme_names or []) + ' ' + (enzyme_class or '') + ' ' + ' '.join(pathways or [])
    all_text = all_text.lower()
    
    # Score pathway categories
    for category, terms in pathway_categories.items():
        if any(term.lower() in all_text for term in terms):
            pathway_category_scores[category] = 1.0
            pathway_score += 1.0
    
    # Use corrosion-related terms from corrosion_keyword_groups
    corrosion_relevant_groups = ['iron_sulfur_redox', 'ocre', 'acid_production', 
                               'electron_transfer', 'biofilm', 'sulfide']
    
    for pathway in pathways or []:
        pathway_lower = pathway.lower()
        for group in corrosion_relevant_groups:
            if group in corrosion_keyword_groups:
                for term in corrosion_keyword_groups[group]:
                    if term.lower() in pathway_lower:
                        pathway_score += 1.0
                        break
    
    # Get organic acid terms from pathway_categories
    organic_terms = []
    if 'organic_acid_metabolism' in pathway_categories:
        organic_terms = pathway_categories['organic_acid_metabolism']
    
    name_text = ' '.join(enzyme_names or []).lower()
    if any(term.lower() in name_text for term in organic_terms):
        pathway_score += 0.5  # Smaller weight for organic terms alone
    
    return pathway_score, pathway_category_scores

def check_metal_organic_synergy(metals, enzyme_names):
    """Check for synergistic effects between metals and organic compounds.
    Args:
        metals: List of metals
        enzyme_names: List of enzyme names
    Returns:
        Synergy score
    """
    name_text = ' '.join(enzyme_names or []).lower()
    
    # Get organic acid terms from pathway_categories
    organic_terms = []
    if 'organic_acid_metabolism' in pathway_categories:
        organic_terms = pathway_categories['organic_acid_metabolism']
    
    if (any(metal in metals for metal in ['iron', 'Fe', 'manganese', 'Mn', 'copper', 'Cu']) and
        any(term.lower() in name_text for term in organic_terms)):
        return 1.0  # Weight for metal-organic combination
    
    return 0.0

def calculate_corrosion_relevance_score(metal_score, mech_score, pathway_score, process_score, 
                                       keyword_score, synergy_score=0, functional_score=0):
    """Calculate final corrosion relevance score and category.
    Args:
        metal_score: Metal involvement score
        mech_score: Corrosion mechanism score
        pathway_score: Pathway relevance score
        process_score: Organic process score
        keyword_score: Keyword group score
        synergy_score: Synergy score
        functional_score: Functional category score
    Returns:
        Corrosion relevance score and category
    """
    # Apply weights
    weighted_metal_score = metal_score * METAL_SCORE_WEIGHT
    weighted_mech_score = mech_score * CORROSION_MECHANISM_WEIGHT
    weighted_process_score = process_score * ORGANIC_PROCESS_WEIGHT
    weighted_keyword_score = keyword_score * KEYWORD_SCORE_WEIGHT
    weighted_synergy_score = synergy_score * SYNERGY_SCORE_WEIGHT
    weighted_functional_score = functional_score * FUNCTIONAL_SCORE_WEIGHT
    
    # Calculate final score
    corrosion_relevance_score = float(
        weighted_metal_score + 
        weighted_mech_score + 
        pathway_score + 
        weighted_process_score + 
        weighted_keyword_score + 
        weighted_synergy_score + 
        weighted_functional_score
    )
    
    # Determine category
    if corrosion_relevance_score >= HIGH_RELEVANCE_THRESHOLD:
        corrosion_relevance = 'high'
    elif corrosion_relevance_score >= MEDIUM_RELEVANCE_THRESHOLD:
        corrosion_relevance = 'medium'
    else:
        corrosion_relevance = 'low'
    
    return corrosion_relevance_score, corrosion_relevance