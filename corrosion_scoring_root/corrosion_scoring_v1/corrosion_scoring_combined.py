# ========== __init__.py ==========
# corrosion_scoring/__init__.py

from .term_processor import TermProcessor
from .scoring_system import (
    calculate_overall_scores,
    calculate_corrosion_relevance_score,
    assign_mechanism_from_pathway,
    assign_corrosion_mechanisms,
    infer_mechanisms_from_pathway_label,
    consolidate_metal_terms,
    validate_against_pathways
)

from .utils_ec import normalize_ec_id, strip_all_ec_tokens
from .name_utils import enhanced_clean_protein_name, clean_protein_name

# Re-export global dictionaries so they are available as `cs.metal_mapping`, etc.
from .global_terms import (
        metal_terms_dict, # For retrieval and scoring
        corrosion_synergies_dict, # For retrieval and scoring
        functional_categories_dict, # For retrieval and scoring
        pathway_dict, mechanisms_dict, operational_environmental_factors_dict, # for retrieval but Not for scoring just for poppulating
        metal_mapping # no for retrieval
    )

__all__ = [
    "TermProcessor",
    "calculate_overall_scores",
    "calculate_corrosion_relevance_score",
    "assign_mechanism_from_pathway",
    "assign_corrosion_mechanisms",
    "infer_mechanisms_from_pathway_label", 
    "consolidate_metal_terms",
    "validate_against_pathways",
    "normalize_ec_id",
    "strip_all_ec_tokens",
    "enhanced_clean_protein_name",
    "clean_protein_name",
    # re-exported globals as names
    "metal_terms_dict",
    "corrosion_synergies_dict",
    "functional_categories_dict",
    "metal_mapping",
    "pathway_dict",
    "mechanisms_dict",
    "operational_environmental_factors_dict",
]

# ========== name_utils.py ==========
import re
import unicodedata

from .utils_ec import normalize_ec_id, strip_all_ec_tokens

GARBAGE_PREFIXES = (
    'transferred to', 'transferred to and', 'deleted', 'obsolete',
    'reclassified as', 'renamed to', 'see '
)

def enhanced_clean_protein_name(name: str) -> str:
    """Canonicalize DE/protein names for grouping/aggregation."""
    if name is None:
        return "uncharacterized protein"
    if not isinstance(name, str):
        name = str(name)

    # Unicode normalize (e.g., fancy dashes)
    name = unicodedata.normalize("NFKD", name)

    # Remove EC tokens (keep one in case we need a fallback)
    ec_token = normalize_ec_id(name)
    name = strip_all_ec_tokens(name)

    # Remove bracketed qualifiers
    name = re.sub(r'\([^)]*\)', '', name)
    name = re.sub(r'\[[^\]]*\]', '', name)

    # Normalize separators and punctuation
    name = re.sub(r'\s*[-/]\s*', '-', name)
    name = re.sub(r'[,:;|]+', ' ', name)
    name = re.sub(r'\s+', ' ', name).strip()

    # Lowercase
    name = name.lower()

    # Standardize common patterns
    replacements = [
        (r'\bacyl carrier protein\b', 'acp'),
        (r'\balcohol\s+dehydrogenase\b', 'alcohol-dehydrogenase'),
        (r'\bglutathione\s+dehydrogenase\b', 'glutathione-dehydrogenase'),
        (r'\b(?:l-)?threonine\s+dehydrogenase\b', 'threonine-dehydrogenase'),
        (r'\b3\-oxoacyl\-acp\s+reductase\b', '3-oxoacyl-acp-reductase'),
    ]
    for pat, repl in replacements:
        name = re.sub(pat, repl, name)
        
    name = re.compile(r'\b([a-z]{3,})\s+(dehydrogenase|reductase|synthase|synthetase)\b', re.IGNORECASE) \
         .sub(lambda m: m.group(0) if '-' in m.group(1) else f"{m.group(1)}-{m.group(2)}", name)

    # Light suffix cleanup
    name = re.sub(r'\s+(domain|fragment|precursor)$', '', name)

    # Remove immediate duplicate tokens only
    tokens = name.split()
    dedup = []
    prev = None
    for t in tokens:
        if t != prev:
            dedup.append(t)
        prev = t
    name = ' '.join(dedup)

    # Final cleanup
    name = re.sub(r'-{2,}', '-', name).strip(' -')

    if not name and ec_token:
        name = f"ec {ec_token}"
    if not name:
        name = "uncharacterized protein"
    return name
#===================================

custom_dict = {
    "bifunctional glutamine-synthetase adenylyltransferase-adenylyl-removing enzyme": "bifunctional glutamine-synthetase adenylyltransferase",
    "glutamine--fructose-6-phosphate aminotransferase": "glutamine--fructose-6-phosphate aminotransferase",
    "udp-3-o-acyl-n-acetylglucosamine deacetylase": "udp-3-o-acyl-n-acetylglucosamine deacetylase",
    "multifunctional oxoglutarate decarboxylase": "oxoglutarate decarboxylase",
    "aldo-keto-reductase 2 family": "aldo-keto-reductase", 
    "coenzyme a biosynthesis bifunctional protein coabc (dna-pantothenate metabolism flavoprotein) (phosphopantothenoylcysteine-synthetase-decarboxy": "coenzyme a biosynthesis protein ppcs-ppcdc"
}

def clean_protein_name(name: str) -> str:
    """Simplify the names of the protein without losing the biological meaning as literature recommends."""
    if not isinstance(name, str):
        name = '' if name is None else str(name)
    # normalise name
    name = name.strip()
    name_lower = name.lower()
    #drop garbage prefixes early
    if any(name_lower.startswith(p) for p in GARBAGE_PREFIXES):
        return ""
    # Step 2: Remove uncertainty terms at the beginning
    uncertainty_terms = list(GARBAGE_PREFIXES) + ['probable','putative','possible','uncharacterized','hypothetical']
    pattern_uncertainty = r'^(?:' + '|'.join(uncertainty_terms) + r')\s+'
    name = re.sub(pattern_uncertainty, '', name, flags=re.IGNORECASE)

    # Step 3: Apply custom dictionary based on prefix match
    for original_start, short_name in custom_dict.items():
        if name_lower.startswith(original_start.lower()):
            name = short_name
            name_lower = name.lower()
            break

    # step 5 cut at ~50 chars but finish the current word
    if len(name) > 50:
        tail = name[50:]
        m = re.search(r'[\s\.\-\[\]\(\)]', tail)
        if m:
            end = 50 + (m.start()) 
        else:               
            end = min(50, len(name))  # hard cut at 50 if no separator found
        name = name[:end].strip()

    return name


# ========== scoring_system.py ==========
#!/usr/bin/env python
# coding: utf-8
###=========================================================
#10.09.2025
###_________________________
"""
Corrosion Relevance Scoring System: functions for corrosion relevance evaluation.
Updated to use unified functional_categories as major scoring category, corrosion_synergies and metal_terms
all the other categories (mechanisms, pathways, operational_ are just to poppulate and make plots), this was done so not to 
inflate the results and therefore functional_categories is a mixture of all the other entities.
"""
import math
import sys
import os
import re
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
def assign_corrosion_mechanisms(text_to_analyze: str, mechanisms_processor) -> list[str]:
    """
    Identifies and extracts corrosion mechanism terms from the given text
    using the module's _mechanism_processor instance.
    This function is for POPULATION, NOT SCORING.
    """
    mechanisms_processor = TermProcessor(mechanisms_dict)
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
    mech_matches = mechanisms_processor.find_all_matches(text_to_analyze)
    return list(mech_matches.keys())

def infer_mechanisms_from_pathway_label(label: str) -> list[str]:
    """
    Maps pathway categories to likely corrosion mechanisms.
    This creates the bridge between pathways and mechanisms.
    """
    if not isinstance(label, str):
        raise TypeError("label must be a string pathway label")
    label = label.strip()
    if not label:
        raise ValueError("label must be non-empty")

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
    if label not in pathway_to_mechanism_map:
            raise KeyError(f"unknown pathway label: {label}")

    inferred = pathway_to_mechanism_map[label]
    valid = set(TermProcessor(cs.mechanisms_dict).category_to_terms.keys())
    return [m for m in inferred if m in valid]
#=================================================================================

def assign_mechanism_from_pathway(text_to_analyze: str, mechanisms_processor) -> list[str]:
    """
    Extracts corrosion mechanisms from pathway text using both pathway and mechanism processors.
    This function looks for direct mechanism terms AND infers mechanisms from pathway names.
    """
    if not isinstance(text_to_analyse, str) or not text_to_analyse.strip():
        return []
  
    mech_matches = mechanisms_processor.find_all_matches(text_to_analyse)# {'cat': ['child1', ...]}
    
    return sorted({t for terms in mech_matches.values() for t in terms})  # child terms only
#==============================================================================================================

def score_by_processor(text: str, processor: TermProcessor, weight_of) -> tuple[dict, dict]:
    if not isinstance(text, str) or not text.strip():
        raise ValueError("text must be a non-empty string")
    if processor is None:
        raise ValueError("processor is required")
    matched = processor.find_all_matches(text)  # {'cat': [child_terms]}
    cat_scores = {}
    for cat, terms in matched.items():
        if not isinstance(terms, list):
            raise TypeError(f"terms for category '{cat}' must be a list")
        hits = len(set(terms))
        cat_scores[cat] = math.log(hits + 1) * float(weight_of(cat))
    return cat_scores, {k: sorted(set(v)) for k, v in matched.items()}

#==============================================================================================================

def consolidate_metal_terms(metals, detected_metal_categories=None) -> list[str]:
    """
    Consolidates metal names from BRENDA and detected categories into standardized symbols.
        Args:    brenda_metals: Raw metal terms from BRENDA data
        detected_metal_categories: Metal category names detected by TermProcessor
    Returns:  list: Consolidated list of unique, standardized metal symbols.
    """   
    consolidated = set()
    
    # Process BRENDA metals (raw terms)
    for metal in (metals or []):
        metal_raw =str(metal)
        metal_norm = metal_raw.strip().lower()
        if 'not detected' in metal_norm or 'none' in metal_norm or not metal_norm:
            continue
        # strip brackets aand the word ion
        metal_norm = re.sub(r'[\[\]\(\)]', '', metal_norm)
        metal_norm = re.sub(r'\bions?\b', '', metal_norm).strip()
        if 'not detected' in metal_norm:
          continue
        raw_tokens = re.findall(r'[a-z0-9\+\-]+', metal_norm)
        tokens = set(re.sub(r'^(fe|cu|zn|ni|co|mn|cr|al|mg|ca|ba|sr|pb|as|hg)\d+\+?$', r'\1', t)
                      for t in raw_tokens)
        if re.search(r'\b\d+\s*fe\s*[-–]\s*\d+\s*s\b', metal_norm):
            tokens.add('fe')
        # Map all matched tokens (no early break)
        matched = False
        for key, symbol in metal_mapping.items():

            if key.lower() in tokens:
                consolidated.add(symbol)
                matched = True
        if not matched:
            for t in tokens:
                if t in {'fe','cu','ni','zn','co','mn','cr','al','mg','ca','ba','sr','pb','as','hg','mo','w','v','ti','sn','sb','cd','ag','au','pt'}:
                    consolidated.add(t.upper())
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
    
    return sorted(consolidated)
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
# --- RATIONALE (FC matching single-pass) ---------------------------------------
# Functional-category matches is done ONCE via fc_processor.find_all_matches(text),
# then score per category from that result. This enforces SINGLE OWNERSHIP of each
# child term (TermProcessor already resolves collisions via priority), preventing
# double-counting the same normalized token across different categories when
# synonyms overlap. Then:
#   1) unique() the child terms per category,
#   2) score with log(hits + 1) * category_weight,
#   3) record the actual CHILD TERMS per category for auditability.
# Looping over every category’s term list and calling matches_normalized()
# could credit the same token to multiple categories. 
def calculate_overall_scores(text, fc_processor, metal_processor, synergy_processor, brenda_metals=None):
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
    results = {}
 
    # Score metals using processor
    metal_score, detected_metals = score_keyword_matches(text, processor= metal_processor)
 
    #========================================================================
    # Score functional categories - PRIMARY scoring method
    # Precompute deduplicated matches once
    fc_matched = fc_processor.find_all_matches(text) 
    functional_score = 0.0
    func_matches = {} #dict storages weighted score
    func_found_terms = {} # actual terms found per category
    detected_fc_names = set()  # to take into account for coocurrence
    for cat, terms in fc_matched.items():
        if not terms:
            continue

        # category_hits = unique child-terms owned by this category
        category_hits = len(set(terms))
        base_score = math.log(category_hits + 1)
        weight = functional_categories_dict.get(cat, {}).get("score", 1.0)
        weighted_score = base_score * weight

        func_matches[cat] = weighted_score
        func_found_terms[cat] = sorted(set(terms))
        detected_fc_names.add(cat)
        functional_score += weighted_score
        
    results["functional_categories"] = [
        {"category": cat, "score": score, "terms": func_found_terms.get(cat, [])}
    for cat, score in func_matches.items()
    ]
    results["functional_score"] = float(functional_score) if functional_score is not None else 0.0
    # flatten list
    results["functional_child_terms"] = sorted({t for terms in func_found_terms.values() for t in terms})

    #=====================================================================

    # The functional categories to check for co-occurrence 
    subcategories_fc = ["o2_consumption","nitrogen_metabolism", "iron_metabolism","sulfur_metabolism","h2_consumption","direct_eet",
    "carbon_metabolism","indirect_eet","organic_acid_metabolism","metal_binding_chelation","biofilm_formation","manganese_processes",
    "methanogenesis","fumarate_formation","halogen_related","phosphorus_metabolism"
    ]
    
    def detect_functional_category_synergies(text, functional_categories_dict, subcategories_fc, fc_processor):
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
        text, functional_categories_dict, subcategories_fc, fc_processor)
    
    # Legacy keyword synergy detection (as fallback) corresponds to the global_terms dictionary corrosion_synergies
    keyword_synergy_score = 0.0
    keyword_synergy_groups_found = set()
    
    for synergy_group, obj in corrosion_synergies_dict.items():
        terms = obj.get('terms', []) if isinstance(obj, dict) else obj
        group_hits = sum(1 for term in terms if synergy_processor.matches_normalized(term, text))
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
    results["metal_score"] = float(metal_score or 0.0)
    results["functional_score"] = float(functional_score) if functional_score is not None else 0.0
    results["synergy_score"] = float(results.get("synergy_score", 0.0))  
   
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
    mechanisms_processor = TermProcessor(mechanisms_dict)
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
            norm_expected = {mechanisms_processor.normalize_text(m) for m in expected_mechanisms}
            norm_detected = {mechanisms_processor.normalize_text(m) for m in detected_mechanisms}
            
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

# ========== term_processor.py ==========
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
    #===================================================================

# ========== utils_ec.py ==========
import re
from typing import Optional
#====================== Normalizing EC =====================

# Accept "EC 1.1.1.1", "1.1.1.1", and hyphens "1.1.1.-"
_EC_RE = re.compile(r'\b(?:EC\s*)?((?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-))\b', re.IGNORECASE)

def normalize_ec_id(s: str) -> Optional[str]:
    """Return normalized EC id 'x.x.x.x' (digits or '-') or None if not found."""
    if not isinstance(s, str):
        return None
    m = _EC_RE.search(s.strip())
    return m.group(1) if m else None

def strip_all_ec_tokens(text: str) -> str:
    """Remove all EC tokens from text."""
    if not isinstance(text, str):
        return ""
    return _EC_RE.sub("", text).strip()
