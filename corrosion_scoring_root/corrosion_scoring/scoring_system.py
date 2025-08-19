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
import re
import math
from collections import defaultdict
from .term_processor import TermProcessor 

try:
    # Try relative import (for package installation)
    from .global_terms import (
        metal_terms,
        corrosion_synergies,
        functional_categories,
        metal_mapping, pathway_categories, corrosion_mechanisms # Not for scoring
    )
except ImportError:
    print("Critical error")#

# Scoring weights - 
METAL_SCORE_WEIGHT = 0.5
FUNCTIONAL_SCORE_WEIGHT = 1.5
SYNERGY_SCORE_WEIGHT = 2.0

# Classification thresholds
HIGH_RELEVANCE_THRESHOLD = 5.0
MEDIUM_RELEVANCE_THRESHOLD = 2.0

def score_keyword_matches(text, keyword_list, processor = None):
    """Score keyword matches in text.
    
    Args:
        text: Text to analyze
        keyword_list: List of keywords to search for
        
    Returns:
        Tuple of (score, matched_terms)
    """
    text_lower = text.lower()
    score = 0.0
    matches = []
    for keyword in keyword_list:
    if processor.matches_normalized(keyword, text_lower):
        matches.append(keyword)
        score += 1.0

    return score, matches

#============================================================================================================
def assign_corrosion_mechanisms(text_to_analyze: str) -> list[str]:
    """
    Identifies and extracts corrosion mechanism terms from the given text
    using the module's _mechanism_processor instance.
    This function is for POPULATION, NOT SCORING.
    """
    mechanism_processor = TermProcessor(corrosion_mechanisms)
    # Use the processor to find all matching terms and their categories in one fast pass
    found_mechanisms = set() # Use a set to ensure unique mechanisms

    # Tokenize the input text into individual words/phrases to process
    # This regex pulls out sequences of word characters (letters, numbers, underscore)
    terms_to_process = re.findall(r'\b\w+\b', text_to_analyze.lower())

    for term in terms_to_process:
        # Use the mechanism_processor to find the highest priority category match
        category = mechanism_processor.find_first_category(term)
        if category:
            found_mechanisms.add(category)

    return list(found_mechanisms)

def infer_mechanisms_from_pathway_category(pathway_category: str) -> list[str]:
    """
    Maps pathway categories to likely corrosion mechanisms.
    This creates the bridge between pathways and mechanisms.
    """
    pathway_to_mechanism_map = {'oxygen_metabolism': ['O2_consumption'], 
        'nitrogen_metabolism':  ['nitrogen_metabolism'],
        'iron_sulfur_redox': ['iron_metabolism', 'sulfur_metabolism', 'ocre_formation'],
        'manganese_processes': ['manganese_metabolism'],
        'sulfur_metabolism':['sulfur_metabolism'],
        'hydrogen_metabolism': ['h2_consumption'],
        'organic_acid_metabolism': ['acid_production'],
        'metal_organic_interaction': ['metal_chelation', 'microbe_metal_synergy'],
        'biofilm_formation': ['biofilm_formation'],
        'carbon_metabolism': ['carbon_metabolism'],
        'halogen_related': ['chloride_attack'],
        'methanogenesis': ['methanogenesis']  
    }
    
    return pathway_to_mechanism_map.get(pathway_category, [])

def assign_mechanism_from_pathway(pathway_text: str) -> list[str]:
    """
    Extracts corrosion mechanisms from pathway text using both pathway and mechanism processors.
    This function looks for direct mechanism terms AND infers mechanisms from pathway names.
    """
    mechanism_processor = TermProcessor(corrosion_mechanisms)
    pathway_processor = TermProcessor(pathway_categories)
    found_mechanisms = set()
    
    if not pathway_text:
        return []
    
    # Tokenize the pathway text
    terms_to_process = re.findall(r'\b\w+\b', pathway_text.lower())
    
    # Method 1: Direct mechanism detection
    for term in terms_to_process:
        mechanism = mechanism_processor.find_first_category(term)
        if mechanism:
            found_mechanisms.add(mechanism)
    
    # Method 2: Pathway-to-mechanism inference using pathway processor
    for term in terms_to_process:
        pathway_category = pathway_processor.find_first_category(term)
        if pathway_category:
            # Map pathway categories to mechanisms 
            inferred_mechanisms = infer_mechanisms_from_pathway_category(pathway_category)
            found_mechanisms.update(inferred_mechanisms)
    
    # Method 3: Pattern-based inference for common pathway types
    pathway_lower = pathway_text.lower()
    
    # Sulfur-related pathways
    if any(term in pathway_lower for term in ['sulfur', 'sulfate', 'thiosulfate', 'sulfide']):
        found_mechanisms.add('sulfur_metabolism')
    
    # Iron-related pathways
    if any(term in pathway_lower for term in ['iron', 'ferric', 'ferrous', 'heme']):
        found_mechanisms.add('iron_metabolism')
    
    # Organic acid pathways
    if any(term in pathway_lower for term in ['acetate', 'lactate', 'citrate', 'organic acid']):
        found_mechanisms.add('acid_production')
    
    # Biofilm-related pathways
    if any(term in pathway_lower for term in ['biofilm', 'exopolysaccharide', 'quorum']):
        found_mechanisms.add('biofilm_formation')
    
    return list(found_mechanisms)
#==============================================================================================================

def consolidate_metal_terms(brenda_metals, text_detected_metals= None):
    """
    Consolidates metal names from BRENDA and text mining into standardized symbols.
    Parameters:  brenda_metals (list of str): Metals obtained from BRENDA data.
        text_detected_metals (list of str): Metals detected from text mining.
    Returns: list: Consolidated list of unique, standardized metal symbols.
    """   
    consolidated = set()
    all_metals = (brenda_metals or []) + (text_detected_metals or [])
    
    for metal in all_metals:
        metal_norm = metal.strip().lower() #processor.matches_normalized(term, text_lower)
        found_mapping = False # Flag to check if a mapping was found
        # Check if the normalized term matches any key in the standard mapping
        for key, symbol in metal_mapping.items():
            if key in metal_norm:
                consolidated.add(symbol)
                found_mapping = True
                #break # Found a mapping, move to the next metal
        if not found_mapping:
            # If no mapping is found after checking all keys, add the original
            consolidated.add(metal.strip())
    return list(consolidated)

def calculate_overall_scores(text, fc_processor, metal_processor, synergy_processor, brenda_metals=None, pathways=None):
    """Calculate all the overall scores for a given text.
    
    Args: text: Text to analyze (combined enzyme names, class, pathways, reactions)
        brenda_metals: Metals from BRENDA database
        pathways: Pathway information for additional mechanism detection
    Returns: Dictionary containing all overall scores and matched categories
    """
    if brenda_metals is None:
        brenda_metals = []

    results = {}
    text_lower = text.lower()
    # Convert pathways to a single string if it's a list for assign_mechanism_from_pathway
    if pathways:
        pathways_for_mechanism = ' '.join(pathways) if isinstance(pathways, list) else (pathways or "")
        # Call assign_mechanism_from_pathway which uses the module's global TermProcessor instance
        results['corrosion_mechanisms'] = assign_mechanism_from_pathway(pathways_for_mechanism)
    else:
        results['corrosion_mechanisms'] = [] # Default if no pathways

    # Score metals with all terms
    all_metal_keywords = [term for sublist in metal_terms.values() for term in sublist]
    metal_score, detected_metals = score_keyword_matches(text, all_metal_keywords, metal_processor)
    
    # Consolidate metals
    consolidated_metals = consolidate_metal_terms(brenda_metals, detected_metals)
    results["consolidated_metals"] = consolidated_metals
    results["overall_metal_score"] = float(metal_score)

    #========================================================================
    # Score functional categories - PRIMARY scoring method
    functional_score = 0.0
    func_matches = {} #dict storages weighted score
    detected_fc_names = set()  # to take into account for coocurrence
    for cat, details in functional_categories.items():
        category_hits = 0
        for term in details["terms"]:
            if fc_processor.matches_normalized(term, text_lower):
                category_hits += 1
        
        if category_hits > 0:
            base_score = math.log(category_hits + 1) #log prevents inflation, promotes balance
            weighted_score = base_score * details["score"]
            func_matches[cat] = weighted_score
            functional_score += weighted_score
            detected_fc_names.add(cat)

    results["functional_categories"] = [
        {"category": cat, "score": score}
        for cat, score in func_matches.items()
    ]
    results["overall_functional_score"] = float(functional_score)
    #=====================================================================

    # The functional categories to check for co-occurrence 
    subcategories_fc = ["o2_consumption","nitrogen_metabolism", "iron_metabolism","sulfur_metabolism","h2_consumption","direct_eet",
    "carbon_metabolism","indirect_eet","organic_acid_metabolism","metal_binding_chelation","biofilm_formation","manganese_processes",
    "methanogenesis","fumarate_formation","halogen_related","phosphorus_metabolism"
    ]
    
    def detect_functional_category_synergies(text_lower, functional_categories, subcategories_fc, fc_processor):
        """ Detect synergies based on co-occurrence of terms from different functional categories.
        Args:text_lower: Lowercase text to analyze
            functional_categories: Dictionary of functional categories and their terms
            subcategories_fc: List of functional category names to check
        Returns: Dict containing synergy detection results"""
        
        # Step 1: Detect which functional categories have hits and collect their terms
        detected_categories = {}
        all_found_terms = set()
        
        for category in subcategories_fc:
            if category in functional_categories:
                terms = functional_categories[category]['terms']
                found_terms = [term for term in terms if fc_processor.matches_normalized(term, text_lower)]
                
                if found_terms:
                    detected_categories[category] = found_terms
                    all_found_terms.update(found_terms)
        
        # Step 2: Check for meaningful synergies
        synergy_results = {
            'fc_cooccurrence_synergy_hit': False,
            'synergy_score': 0.0,
            'synergy_child_terms_found': [],
            'synergy_categories_involved': [],
            'synergy_description': ''
        }
        
        # Require at least 2 different functional categories
        if len(detected_categories) < 2:
            return synergy_results
        
        # Define high-priority synergy combinations with their scores
        priority_synergies = {
            ('organic_acid_metabolism', 'metal_binding_chelation'): {
                'score': 3.0,
                'description': 'TOC-chelation Synergy (TOC-chelate)'
            },
            ('iron_metabolism', 'organic_acid_metabolism'): {
                'score': 2.8,
                'description': 'Iron-Organic Acid Synergy (acid-enhanced Fe corrosion)'
            },
            ('biofilm_formation', 'metal_binding_chelation'): {
                'score': 2.7,
                'description': 'biofilm-chelate Synergy (biofilm-chelate-corrosion)'
            },
            ('o2_consumption', 'iron_metabolism'): {
                'score': 2.5,
                'description': 'Oxygen-Iron Synergy (aerobic Fe corrosion)'
            },
            ('sulfur_metabolism', 'iron_metabolism'): {
                'score': 2.4,
                'description': 'Sulfur-iron Synergy (SRB-mediated corrosion)'
            },
            ('sulfur_metabolism', 'h2_consumption'): {
                'score': 2.3,
                'description': 'Sulfur-Hydrogen Synergy (SRB-mediated corrosion)'
            },
            ('biofilm_formation', 'iron_metabolism'): {
                'score': 2.2,
                'description': 'Biofilm-Iron Synergy (biofilm-enhanced Fe corrosion)'
            },
            ('nitrogen_metabolism', 'iron_metabolism'): {
                'score': 2.0,
                'description': 'Nitrogen-Iron Synergy (nitrate-enhanced Fe corrosion)'
            }
        }
        
        # Step 3: Check for priority synergies first
        max_synergy_score = 0.0
        best_synergy = None
        involved_categories = []
        
        for synergy_pair, synergy_info in priority_synergies.items():
            cat1, cat2 = synergy_pair
            if cat1 in detected_categories and cat2 in detected_categories:
                # Calculate combined term count from both categories
                combined_terms = set(detected_categories[cat1] + detected_categories[cat2])
                
                # Require at least 3 terms total for high-confidence synergy
                if len(combined_terms) >= 3:
                    current_score = synergy_info['score']
                    if current_score > max_synergy_score:
                        max_synergy_score = current_score
                        best_synergy = synergy_info
                        involved_categories = [cat1, cat2]
        
        # Step 4: If no priority synergy found, check for general multi-category synergy
        if max_synergy_score == 0.0 and len(detected_categories) >= 2:
            # General synergy: any 2+ categories with sufficient terms
            if len(all_found_terms) >= 4:  # Require more terms for general synergy
                max_synergy_score = 1.5  # Lower score for general synergy
                involved_categories = list(detected_categories.keys())
                best_synergy = {
                    'description': f'Multi-pathway Synergy ({len(detected_categories)} categories)'
                }
        
        # Step 5: Populate results if synergy detected
        if max_synergy_score > 0.0:
            # Collect all terms from involved categories
            synergy_terms = set()
            for cat in involved_categories:
                if cat in detected_categories:
                    synergy_terms.update(detected_categories[cat])
            
            synergy_results.update({
                'fc_cooccurrence_synergy_hit': True,
                'synergy_score': max_synergy_score,
                'synergy_child_terms_found': sorted(list(synergy_terms)),
                'synergy_categories_involved': involved_categories,
                'synergy_description': best_synergy['description']
            })
        
        return synergy_results

    #=================
    synergy_results = detect_functional_category_synergies(
    text_lower, functional_categories, subcategories_fc, synergy_processor
    )
    
    # Legacy keyword synergy detection (as fallback)
    keyword_synergy_score = 0.0
    keyword_synergy_groups_found = set()
    
    for synergy_group, terms in corrosion_synergies.items():
        group_hits = 0
        for term in terms:
            if synergy_processor.matches_normalized(term, text_lower):
                group_hits += 1
        if group_hits > 0:
            keyword_synergy_score += math.log(group_hits + 1)
            keyword_synergy_groups_found.add(synergy_group)

    # Choose best synergy result
    if synergy_results['fc_cooccurrence_synergy_hit']:
        # Use functional category synergy
        results["corrosion_synergies"] = synergy_results['synergy_child_terms_found']
        results["overall_synergy_score"] = float(synergy_results['synergy_score'])
        results["synergy_type"] = "functional_category_cooccurrence"
        results["synergy_description"] = synergy_results['synergy_description']
        results["synergy_categories"] = synergy_results['synergy_categories_involved']
    elif keyword_synergy_score > 0:
        # Fall back to keyword synergy
        results["corrosion_synergies"] = sorted(list(keyword_synergy_groups_found))
        results["overall_synergy_score"] = float(keyword_synergy_score)
        results["synergy_type"] = "keyword_based"
        results["synergy_description"] = "Keyword-based synergy detection"
    else:
        # No synergy detected
        results["corrosion_synergies"] = []
        results["overall_synergy_score"] = 0.0
        results["synergy_type"] = "none"
        results["synergy_description"] = "No synergies detected"
   
    return results

def calculate_corrosion_relevance_score(
    metal_score,
    synergy_score=0,
    functional_score=0,
    ):
    """Calculate final corrosion relevance score and category.
    UPDATED: functional_score now carries more weight as primary scoring method
    Args:functional_score: Functional category score (PRIMARY)
        metal_score: Metal involvement score
        synergy_score: Synergy score  
    Returns:  Corrosion relevance score and category
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
#=================================================================================
# function for fc in brenda
# This new function uses the processor but mimics the old function's output
def calculate_scores_with_processor(text, fc_processor, original_fc_dict):
    """
    Calculates granular corrosion scores using the pre-compiled TermProcessor
    but returns the same data structure as the original calculate_overall_scores.
    """
    # 1. Use the processor to find all matching terms and their categories in one fast pass
    categorized_matches = fc_processor.find_all_matches(text)
    
    # 2. Now, replicate the original scoring logic using these matches
    functional_score = 0.0
    func_matches_output = [] # This will be the list of dicts for the output
    
    for category, found_terms in categorized_matches.items():
        # Count the number of unique terms found for this category
        num_hits = len(found_terms)
        
        # Get the base score for this category from the original dictionary
        base_score = original_fc_dict[category]['score']
        
        # Apply the same log scaling as original function
        weighted_score = math.log(num_hits + 1) * base_score
        
        functional_score += weighted_score
        func_matches_output.append({"category": category, "score": weighted_score})
        
    # 3. Return the results in the EXACT same format as before
    return {
        "functional_categories": func_matches_output,
        "overall_functional_score": float(functional_score)
    }

#=====================================================================================================================0
def validate_against_pathways(record, pathways_data):
    """
    Validates detected pathways, mechanisms, and functional categories against ipath ground truth.
    Returns validation metrics and suggestions.
    """
    mechanism_processor = TermProcessor(corrosion_mechanisms)
    pathway_processor = TermProcessor(pathway_categories)
    ec_number = record['ec_number']
    validation_results = {
        'pathway_validation': {},
        'mechanism_validation': {},
        'functional_category_validation': {},
        'overall_confidence': 0.0
    }
    
    if ec_number not in pathways_data:
        validation_results['overall_confidence'] = 0.5  # No ground truth available
        return validation_results
    
    pathways = pathways_data[ec_number]
    detected_pathways = record.get('pathway_ko', [])
    detected_mechanisms = record.get('corrosion_mechanisms', [])
    detected_fc = [fc['category'] for fc in record.get('functional_categories', [])]
    
    # Pathway validation
    if pathways and detected_pathways:
        # Normalize both for comparison
        norm_pathways = {pathway_processor._normalize_term(p) for p in pathways}
        norm_detected = {pathway_processor._normalize_term(p) for p in detected_pathways}
        
        overlap = norm_pathways.intersection(norm_detected)
        pathway_precision = len(overlap) / len(norm_detected) if norm_detected else 0
        pathway_recall = len(overlap) / len(norm_pathways) if norm_pathways else 0
        
        validation_results['pathway_validation'] = {
            'precision': pathway_precision,
            'recall': pathway_recall,
            'f1_score': 2 * (pathway_precision * pathway_recall) / (pathway_precision + pathway_recall) if (pathway_precision + pathway_recall) > 0 else 0,
            'overlap_terms': list(overlap)
        }
    
    # Mechanism validation (infer expected mechanisms from ipath)
    if pathways:
        expected_mechanisms = set()
        for pathway in pathways:
            expected_mechanisms.update(assign_mechanism_from_pathway(pathway))
        
        if expected_mechanisms and detected_mechanisms:
            norm_expected = {mechanism_processor._normalize_term(m) for m in expected_mechanisms}
            norm_detected = {mechanism_processor._normalize_term(m) for m in detected_mechanisms}
            
            overlap = norm_expected.intersection(norm_detected)
            mech_precision = len(overlap) / len(norm_detected) if norm_detected else 0
            mech_recall = len(overlap) / len(norm_expected) if norm_expected else 0
            
            validation_results['mechanism_validation'] = {
                'precision': mech_precision,
                'recall': mech_recall,
                'f1_score': 2 * (mech_precision * mech_recall) / (mech_precision + mech_recall) if (mech_precision + mech_recall) > 0 else 0,
                'expected_mechanisms': list(expected_mechanisms),
                'overlap_terms': list(overlap)
            }
    
    # Calculate overall confidence
    pathway_f1 = validation_results['pathway_validation'].get('f1_score', 0)
    mechanism_f1 = validation_results['mechanism_validation'].get('f1_score', 0)
    validation_results['overall_confidence'] = (pathway_f1 + mechanism_f1) / 2
    
    return validation_results