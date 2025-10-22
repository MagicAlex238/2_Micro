"""
SynergyDetector class for detecting corrosion synergies.
"""

import pandas as pd
from typing import Dict, Any, List, Set, Tuple, Optional
from .config import ScoringConfig
from .exceptions import SynergyDetectionError

class SynergyDetector:
    """
    Handles detection and scoring of corrosion synergies.
    
    This class analyzes functional category co-occurrence patterns
    to identify synergistic corrosion mechanisms.
    """
    
    def __init__(self, config: ScoringConfig = None):
        """ Initialize the SynergyDetector.
        Args: config: Configuration settings for synergy detection
        """
        self.config = config or ScoringConfig()
        
        # Functional categories to check for co-occurrence
        self.subcategories_fc = [
            "o2_consumption", "nitrogen_metabolism", "iron_metabolism",
            "sulfur_metabolism", "h2_consumption", "direct_eet",
            "carbon_metabolism", "indirect_eet", "organic_acid_metabolism",
            "metal_binding_chelation", "biofilm_formation", "manganese_processes",
            "methanogenesis", "fumarate_formation", "halogen_related",
            "phosphorus_metabolism"
        ]
    
    def detect_synergies(self, features: Dict[str, Any]) -> Dict[str, Any]:
        """
        Detect synergies from extracted features.        
        Args: features: Dictionary containing extracted features
        Returns: Dictionary containing synergy detection results
        Raises:SynergyDetectionError: If synergy detection fails
        """
        try:
            # Get functional category matches
            fc_matches = features.get('functional_matches_detailed', {})
            
            # Detect functional category synergies
            fc_synergy_results = self._detect_functional_category_synergies(fc_matches)
            
            # Calculate overall synergy score with weighting
            overall_synergy_score = fc_synergy_results['synergy_score'] * self.config.synergy_score_weight
            
            return {
                **fc_synergy_results,
                'overall_synergy_score': float(overall_synergy_score)
            }
            
        except Exception as e:
            raise SynergyDetectionError(f"Synergy detection failed: {e}") from e
    
    def _detect_functional_category_synergies(self, fc_matches: Dict[str, List[str]]) -> Dict[str, Any]:
        """
        Detect synergies based on co-occurrence of terms from different functional categories.
        Args:fc_matches: Dictionary of functional category matches
        Returns: Dictionary containing synergy detection results
        """
        try:
            # Step 1: Filter relevant categories and collect their terms
            detected_categories = {}
            all_found_terms = set()
            
            for category in self.subcategories_fc:
                if category in fc_matches and fc_matches[category]:
                    detected_categories[category] = list(fc_matches[category])
                    all_found_terms.update(fc_matches[category])
            
            # Step 2: Initialize synergy results
            synergy_results = {
                'fc_cooccurrence_synergy_hit': False,
                'synergy_score': 0.0,
                'synergy_child_terms_found': [],
                'synergy_categories_involved': [],
                'synergy_description': '',
                'synergy_type': 'none'
            }
            
            # Step 3: Check minimum requirements
            if len(detected_categories) < self.config.min_synergy_categories:
                return synergy_results
            
            # Step 4: Check for priority synergies
            max_synergy_score = 0.0
            best_synergy = None
            involved_categories = []
            
            for synergy_pair, synergy_info in self.config.priority_synergies.items():
                cat1, cat2 = synergy_pair
                if cat1 in detected_categories and cat2 in detected_categories:
                    # Calculate combined term count from both categories
                    combined_terms = set(detected_categories[cat1] + detected_categories[cat2])
                    
                    # Require minimum terms for high-confidence synergy
                    if len(combined_terms) >= self.config.min_synergy_terms:
                        current_score = synergy_info['score']
                        if current_score > max_synergy_score:
                            max_synergy_score = current_score
                            best_synergy = synergy_info
                            involved_categories = [cat1, cat2]
            
            # Step 5: Check for general multi-category synergy if no priority synergy found
            if max_synergy_score == 0.0 and len(detected_categories) >= self.config.min_synergy_categories:
                if len(all_found_terms) >= self.config.min_general_synergy_terms:
                    max_synergy_score = 1.5  # Lower score for general synergy
                    involved_categories = list(detected_categories.keys())
                    best_synergy = {
                        'description': f'Multi-pathway Synergy ({len(detected_categories)} categories)'
                    }
            
            # Step 6: Populate results if synergy detected
            if max_synergy_score > 0.0:
                synergy_terms = set()
                for cat in involved_categories:
                    if cat in detected_categories:
                        synergy_terms.update(detected_categories[cat])
                
                synergy_results.update({
                    'fc_cooccurrence_synergy_hit': True,
                    'synergy_score': max_synergy_score,
                    'synergy_child_terms_found': sorted(list(synergy_terms)),
                    'synergy_categories_involved': involved_categories,
                    'synergy_description': best_synergy['description'],
                    'synergy_type': 'functional_category_cooccurrence'
                })
            
            return synergy_results
            
        except Exception as e:
            raise SynergyDetectionError(f"Functional category synergy detection failed: {e}") from e
    
    def get_synergy_explanation(self, synergy_categories: List[str]) -> str:
        """
        Get a detailed explanation of detected synergies.
        Args:synergy_categories: List of involved synergy categories
        Returns:Detailed explanation string
        """
        if not synergy_categories or len(synergy_categories) < 2:
            return "No significant synergies detected."
        
        # Check if it's a known priority synergy
        if len(synergy_categories) == 2:
            synergy_pair = tuple(sorted(synergy_categories))
            if synergy_pair in self.config.priority_synergies:
                return self.config.priority_synergies[synergy_pair]['description']
        
        # General explanation for multi-category synergies
        return f"Multi-pathway synergy detected involving {len(synergy_categories)} functional categories: {', '.join(synergy_categories)}"
    
    #===============================================================
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

    
    # Short name mapping for functional categories
    name_map = {
        'iron_metabolism': 'Fe_met',
        'sulfur_metabolism': 'S_met',
        'organic_acid_metabolism': 'organic_acid',
        'metal_binding_chelation': 'chelation',
        'biofilm_formation': 'biofilm',
        'o2_consumption': 'O2',
        'nitrogen_metabolism': 'N',
        'h2_consumption': 'H2',
        'carbon_metabolism': 'C_met',
        'manganese_processes': 'Mn_met',
        'methanogenesis': 'methanogenesis', 
        'fumarate_formation': 'fumarate',
        'phosphorus_metabolism': 'P_met'
    }
    def _detect_row_synergies(self, row: pd.Series) -> Dict[str, Any]:
        """    Collect priority subcategories from 3 columns.        """
        # Priority synergy definitions # Priority lists for each dimension
        priority_functional = ['o2_consumption', 'nitrogen_metabolism', 'h2_consumption','iron_metabolism','sulfur_metabolism', 
                            'organic_acid_metabolism','carbon_metabolism', 'manganese_processes', 'methanogenesis', 'fumarate_formation', 'phosphorus_metabolism',
                                'metal_binding_chelation', 'biofilm_formation']
        
        priority_metals = ['Fe', 'S', 'Cl', 'Mn', 'Ni', 'Cr']
        priority_operational = ['halogen_related', 'microaerobic_conditions','ph_modulation','direct_eet',
                    'indirect_eet','exoelectrogenesis','enzymatic_corrosion','dealloying','galvanic_corrosion','chloride_attack','microbe_metal_synergy','cathodic_depolarization',
                    'passivity_breakdown','concentration_cells','syntrophic_interactions'
                ]
        all_subcategories = []
    
        # 1. Check functional_sub (single string)
        functional_sub = row.get('functional_sub', '')
        if pd.notna(functional_sub) and str(functional_sub).strip():
            func_val = str(functional_sub).strip()
            if func_val in priority_functional:
                all_subcategories.append(func_val)
        
        # 2. Check consolidated_metals (semicolon-separated string)
        consolidated_metals = row.get('consolidated_metals', '')
        if pd.notna(consolidated_metals) and str(consolidated_metals).strip():
            metal_items = [m.strip() for m in str(consolidated_metals).split(';')]
            for metal in metal_items:
                if metal in priority_metals:
                    all_subcategories.append(metal)
        
        # 3. Check operational_sub (single string)
        operational_sub = row.get('operational_sub', '')
        if pd.notna(operational_sub) and str(operational_sub).strip():
            ope_val = str(operational_sub).strip()
            if ope_val in priority_operational:
                all_subcategories.append(ope_val)
        
        # Weighted by importance 
        high_priority = {'Fe', 'S', 'Cl', 'organic_acid_metabolism', 'enzymatic_corrosion', 'microbe_metal_synergy',
                          'iron_metabolism', 'biofilm_formation', 'metal_binding_chelation'}
        total_score = sum(1.5 if item in high_priority else 1.0 
                          for item in all_subcategories)

        return {
            "synergy_combi": all_subcategories if all_subcategories else None,
            "synergy_combi_score": total_score
        }