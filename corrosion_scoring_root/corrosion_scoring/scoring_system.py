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
        metal_terms_dict,
        corrosion_synergies_dict,
        functional_categories_dict,
        metal_mapping, pathway_dict, mechanisms_dict, operational_environmental_factors_dict# Not for scoring
    )
except ImportError as e:
    print("Critical error: Failed to import global_terms")
    raise ImportError("Failed to import global_terms") from e

# Scoring weights - 
METAL_SCORE_WEIGHT = 0.5
FUNCTIONAL_SCORE_WEIGHT = 1.5
SYNERGY_SCORE_WEIGHT = 2.0

#============================================================================================================
def score_keyword_matches(text, processor):
    """Score keyword matches using TermProcessor.
    Args: text: Text to analyze
    processor: TermProcessor instance
    Returns: Tuple of (score, matched_categories)
    """
    if processor is None:
        raise ValueError("processor is None in score_keyword_matches")

    categorized_matches = processor.find_all_matches(text) or {}
    total_score = 0.0
    matched_categories: list[str] = []

    for category, found_terms in categorized_matches.items():
        if found_terms:
            total_score += len(found_terms)
            # Always record the category if it appears
            if category not in matched_categories:
                matched_categories.append(category)

    return total_score, matched_categories

#============================================================================================================
def assign_corrosion_mechanisms(text_to_analyze: str) -> list[str]:
    """
    Identifies and extracts corrosion mechanism terms from the given text
    using the module's _mechanism_processor instance.
    This function is for POPULATION, NOT SCORING.
    """
    mechanism_processor = TermProcessor(corrosion_mechanisms)
    # Use the processor to find all matching terms and their categories in one fast pass
    '''    found_mechanisms = set() # Use a set to ensure unique mechanisms

    # Tokenize the input text into individual words/phrases to process
    # This regex pulls out sequences of word characters (letters, numbers, underscore)
    terms_to_process = re.findall(r'\b\w+\b', text_to_analyze.lower())

    for term in terms_to_process:
        # Use the mechanism_processor to find the highest priority category match
        category = mechanism_processor.find_first_category(term)
        if category:
            found_mechanisms.add(category)

    return list(found_mechanisms)'''
    if not isinstance(text_to_analyze, str) or not text_to_analyze.strip():
        return []
    # Full‑text matching so multi-word phrases are detected
    mech_matches = mechanism_processor.find_all_matches(text_to_analyze)
    return list(mech_matches.keys())

def infer_mechanisms_from_pathway_category(pathway_category: str) -> list[str]:
    """
    Maps pathway categories to likely corrosion mechanisms.
    This creates the bridge between pathways and mechanisms.
    """
    pathway_to_mechanism_map = {'oxygen_metabolism': ['o2_consumption'], 
        'nitrogen_metabolism':  ['nitrogen_metabolism'],
        'iron_sulfur_redox': ['iron_metabolism', 'sulfur_metabolism', 'ochre_formation'],
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
    
    inferred = pathway_to_mechanism_map.get(pathway_dict, [])
    valid_mech_keys = set(TermProcessor(mechanisms_dict).category_to_terms.keys())
    return [m for m in inferred if m in valid_mech_keys]
#=================================================================================

def assign_mechanism_from_pathway(pathway_text: str) -> list[str]:
    """
    Extracts corrosion mechanisms from pathway text using both pathway and mechanism processors.
    This function looks for direct mechanism terms AND infers mechanisms from pathway names.
    """
    mechanism_processor = TermProcessor(mechanisms_dict)
    pathway_processor = TermProcessor(pathway_dict)
      
    found_mechanisms: set[str] = set()
    if not isinstance(pathway_text, str) or not pathway_text.strip():
         return []
    
    # Tokenize the pathway text
    terms_to_process = re.findall(r'\b\w+\b', pathway_text.lower())
   
    # Method 1: Direct mechanism detection (full-text)
    mech_matches = mechanism_processor.find_all_matches(pathway_text)
    found_mechanisms.update(mech_matches.keys())

    # Method 2: Pathway-to-mechanism inference (full-text)
    if not found_mechanisms:
        path_matches = pathway_processor.find_all_matches(pathway_text)
        for pathway_category in path_matches.keys():
            inferred = infer_mechanisms_from_pathway_category(pathway_category)
            found_mechanisms.update(inferred)
    
    return list(found_mechanisms)
#==============================================================================================================

def consolidate_metal_terms(brenda_metals, detected_metal_categories=None):
    """
    Consolidates metal names from BRENDA and detected categories into standardized symbols.
        Args:    brenda_metals: Raw metal terms from BRENDA data
        detected_metal_categories: Metal category names detected by TermProcessor
    Returns:  list: Consolidated list of unique, standardized metal symbols.
    """   
    consolidated = set()
    
    # Process BRENDA metals (raw terms)
    for metal in (brenda_metals or []):
        metal_raw =str(metal)
        metal_norm = metal_raw.strip().lower()
        # strip brackets aand the word ion
        metal_norm = re.sub(r'[\[\]\(\)]', '', metal_norm)
        metal_norm = re.sub(r'\bions?\b', '', metal_norm).strip()
        raw_tokens = re.findall(r'[a-z0-9\+\-]+', metal_norm)
        tokens = set(re.sub(r'^(fe|cu|zn|ni|co|mn|cr|al|mg|ca|ba|sr|pb|as|hg)\d+\+?$', r'\1', t)
                      for t in raw_tokens)
        # Map all matched tokens (no early break)
        matched = False
        for key, symbol in metal_mapping.items():
            if key.lower() in tokens:
                consolidated.add(symbol)
                matched = True
        if not matched:
            continue

    # Process detected categories (category names like 'iron', 'copper')
    for category in (detected_metal_categories or []):
        # Map category names directly to symbols
        if category in metal_mapping:
            consolidated.add(metal_mapping[category])
        else:
            consolidated.add(category)  # Fallback to category name
    
    return list(consolidated)
#========================================================
def sanitize_input_text(text: str, max_length: int = 10000) -> str:
    """Sanitize input text to prevent regex attacks and excessive processing."""
    if not isinstance(text, str):
        return ""
    
    # Limit length to prevent DoS
    if len(text) > max_length:
        text = text[:max_length]
    
    # Remove potentially problematic characters for regex
    text = re.sub(r'[\\^$.*+?{}[\]|()\x00-\x1f]', ' ', text)
    
    # Collapse whitespace
    text = re.sub(r'\s+', ' ', text).strip()
    
    return text

#======================================================================================================================
def calculate_overall_scores(text, fc_processor, metal_processor, synergy_processor, brenda_metals=None, pathways=None):
    """Calculate all the overall scores for a given text.
    Args: text: Text to analyze (combined enzyme names, class, pathways, reactions)
        brenda_metals: Metals from BRENDA database
        pathways: Pathway information for additional mechanism detection
    Returns: Dictionary containing all overall scores and matched categories
    """
    if not isinstance(text, str):
            raise ValueError("text must be a string")
    text = sanitize_input_text(text)
    if not text.strip():
        raise ValueError("text cannot be empty")
    if fc_processor is None or metal_processor is None or synergy_processor is None:
        raise ValueError("All processors must be provided")
    
    if brenda_metals is None:
        brenda_metals = []
    if any(p is None for p in (fc_processor, metal_processor, synergy_processor)):
        raise ValueError("fc_processor, metal_processor, and synergy_processor must be provided")
    text = '' if text is None else (text if isinstance(text, str) else str(text))

    results = {}
 
    # Convert pathways to a single string if it's a list for assign_mechanism_from_pathway
    if pathways:
        pathways_for_mechanism = ' '.join(pathways) if isinstance(pathways, list) else (pathways or "")
        # Call assign_mechanism_from_pathway which uses the module's global TermProcessor instance
        results['corrosion_mechanisms'] = assign_mechanism_from_pathway(pathways_for_mechanism)
    else:
        results['corrosion_mechanisms'] = [] # Default if no pathways

    # Score metals using processor
    metal_score, detected_metals = score_keyword_matches(text, processor= metal_processor)
    
    # Consolidate metals
    consolidated_metals = consolidate_metal_terms(brenda_metals, detected_metals)
    results["consolidated_metals"] = consolidated_metals
    results["metal_score"] = float(metal_score) if metal_score is not None else 0.0

    #========================================================================
    # Score functional categories - PRIMARY scoring method
    functional_score = 0.0
    func_matches = {} #dict storages weighted score
    func_found_terms = {} # actual terms found per category
    detected_fc_names = set()  # to take into account for coocurrence
    for cat, details in functional_categories_dict.items():
        category_hits = 0
        found_terms_for_cat = []
        for term in details["terms"]:
            if fc_processor.matches_normalized(term, text):
                category_hits += 1
                found_terms_for_cat.append(term)
        
        if category_hits > 0:
            base_score = math.log(category_hits + 1) #log prevents inflation, promotes balance
            weighted_score = base_score * details["score"]
            func_matches[cat] = weighted_score
            functional_score += weighted_score
            detected_fc_names.add(cat)
            func_found_terms[cat] = found_terms_for_cat

    results["functional_categories"] = [
        {"category": cat, "score": score}
        for cat, score in func_matches.items()
    ]
    results["functional_score"] = float(functional_score) if functional_score is not None else 0.0

    #=====================================================================

    # The functional categories to check for co-occurrence 
    subcategories_fc = ["o2_consumption","nitrogen_metabolism", "iron_metabolism","sulfur_metabolism","h2_consumption","direct_eet",
    "carbon_metabolism","indirect_eet","organic_acid_metabolism","metal_binding_chelation","biofilm_formation","manganese_processes",
    "methanogenesis","fumarate_formation","halogen_related","phosphorus_metabolism"
    ]
    
    def detect_functional_category_synergies(text, functional_categories, subcategories_fc, fc_processor):
        """ Detect synergies based on co-occurrence of terms from different functional categories.
        Args:text: Lowercase text to analyze
            functional_categories: Dictionary of functional categories and their terms
            subcategories_fc: List of functional category names to check
        Returns: Dict containing synergy detection results"""
        
        # Step 1: Detect which functional categories have hits and collect their terms
        detected_categories = {} # coming from functional categories only dict of category -> list of terms
        all_found_terms = set() # coming from functional categories only set of all terms found
        
        for category in subcategories_fc:
             if category in func_found_terms and func_found_terms[category]:
                detected_categories[category] = list(func_found_terms[category])
                all_found_terms.update(func_found_terms[category])

        # Step 2: Check for meaningful synergies based on detected categories
        synergy_results = {
            'fc_cooccurrence_synergy_hit': False,
            'synergy_score': 0.0,
            'synergy_child_terms_found': [], # list of terms contributing to synergy found inside fc
            'synergy_categories_involved': [], # list of terms in the upper category of the child term
            'synergy_description': '' # description of the synergy detected
        }
        
        # Require at least 2 different functional categories
        if len(detected_categories) < 2:
            return synergy_results
        #Synergy scores based on corrosion failure analysis literature:
        # Organic acid + metal chelation: highest damage potential (score: 3.0) # For reference and justification please see main manuscript
        # Defining high-priority synergy combinations with their scores based on functional categories
        priority_synergies = { # Based on failure analysis of cooling and heating operational water systems
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
        
        # Step 3: Check for priority synergies first and select the highest scoring one
        max_synergy_score = 0.0
        best_synergy = None
        involved_categories = []
        
        for synergy_pair, synergy_info in priority_synergies.items():
            cat1, cat2 = synergy_pair
            if cat1 in detected_categories and cat2 in detected_categories:
                # Calculate combined term count from both categories
                combined_terms = set(detected_categories[cat1] + detected_categories[cat2])
                
                # Require at least 2 terms total for high-confidence synergy
                if len(combined_terms) >= 2:
                    current_score = synergy_info['score']
                    if current_score > max_synergy_score:
                        max_synergy_score = current_score
                        best_synergy = synergy_info
                        involved_categories = [cat1, cat2]
        
        # Step 4: If no priority synergy found, check for general multi-category synergy
        if max_synergy_score == 0.0 and len(detected_categories) >= 2:
            # General synergy: any 2+ categories with sufficient terms
            if len(all_found_terms) >= 3:  # Require more terms for general synergy
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
    
    synergy_results = detect_functional_category_synergies(
        func_found_terms, subcategories_fc)
    
    # Legacy keyword synergy detection (as fallback) corresponds to the global_terms dictionary corrosion_synergies
    keyword_synergy_score = 0.0
    keyword_synergy_groups_found = set()
    
    for synergy_group, obj in corrosion_synergies_dict.items():
        term_list = obj.get('terms', []) if isinstance(obj, dict) else obj
        group_hits = 0
        terms = obj.get('terms', [])
        for term in terms:
            if synergy_processor.matches_normalized(term, text):
                group_hits += 1
        if group_hits > 0:
            keyword_synergy_score += math.log(group_hits + 1)
            keyword_synergy_groups_found.add(synergy_group)

    # Choose best synergy result
    if synergy_results['fc_cooccurrence_synergy_hit']:
        # Use functional category synergy
        # in calculate_overall_scores(), functional-category synergy branch
        results["corrosion_synergies"] = synergy_results['synergy_categories_involved']   # e.g., ['organic_acid_metabolism','metal_binding_chelation']
        results["synergy_terms"] = sorted(synergy_results['synergy_child_terms_found'])   # keep children separately

        results["synergy_score"] = float(synergy_results['synergy_score'])
        results["synergy_type"] = "functional_category_cooccurrence"
        results["synergy_description"] = synergy_results['synergy_description']
        results["synergy_categories"] = synergy_results['synergy_categories_involved']
    elif keyword_synergy_score > 0:
        # Fall back to keyword synergy
        results["corrosion_synergies"] = sorted(list(keyword_synergy_groups_found))
        results["synergy_terms"] = []
        results["synergy_score"] = float(keyword_synergy_score)
        results["synergy_type"] = "keyword_based"
        results["synergy_description"] = "Keyword-based synergy detection"
        results["synergy_categories"] = sorted(list(keyword_synergy_groups_found))
    else:
        # No synergy detected
        results["corrosion_synergies"] = []
        results["synergy_terms"] = [] 
        results["synergy_score"] = 0.0
        results["synergy_type"] = "none"
        results["synergy_description"] = "No synergies detected"
        results["synergy_categories"] = []

    # save raw scores 
    results["metal_score"] = float(metal_score) if metal_score is not None else 0.0
    results["functional_score"] = float(functional_score) if functional_score is not None else 0.0
    results["synergy_score"] = float(results.get["synergy_score", 0.0])  
   
    # calculate overall scores
    results["overall_metal_score"] = float(metal_score * METAL_SCORE_WEIGHT)
    results["overall_functional_score"] = float(functional_score * FUNCTIONAL_SCORE_WEIGHT)
    results["overall_synergy_score"] = float(results["synergy_score"] * SYNERGY_SCORE_WEIGHT)

    return results

#======================================================================================================

def calculate_corrosion_relevance_score(overall_metal_score, overall_synergy_score, overall_functional_score):
    
    """Calculate final corrosion relevance score (float).
   
    Returns:  Corrosion relevance score 
    """
    corrosion_relevance_score =   float(overall_metal_score + overall_synergy_score + overall_functional_score)
    
    return corrosion_relevance_score #, corrosion_relevance

#=====================================================================================================================0
def validate_against_pathways(record, pathways_data):
    """
    Validates detected pathways, mechanisms, and functional categories against ipath ground truth.
    Returns validation metrics and suggestions.
    """
    mechanism_processor = TermProcessor(mechanisms_dict)
    pathway_processor = TermProcessor(pathway_dict)
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
    if isinstance(detected_pathways, str):
        detected_pathways = [detected_pathways]
    detected_mechanisms = record.get('corrosion_mechanisms', [])
    #detected_fc = [fc['category'] for fc in record.get('functional_categories', [])]
    
    # Pathway validation
    if pathways and detected_pathways:
        # Normalize both for comparison
        norm_pathways = {pathway_processor.normalize_text(p) for p in pathways}
        norm_detected = {pathway_processor.normalize_text(p) for p in detected_pathways}
        
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
            norm_expected = {mechanism_processor.normalize_text(m) for m in expected_mechanisms}
            norm_detected = {mechanism_processor.normalize_text(m) for m in detected_mechanisms}
            
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