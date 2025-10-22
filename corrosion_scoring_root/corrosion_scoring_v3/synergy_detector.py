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
    #===============================================================
    # Second option chat
    def detect_synergies_from_dataframe(
        self,
        df: pd.DataFrame,
        feature_columns: List[str]
    ) -> pd.DataFrame:
        """
        Detects corrosion-related synergies based only on priority_synergies.
        Returns df with all original columns + 2 new ones.
        """
        try:
            available_columns = [col for col in feature_columns if col in df.columns]
            synergy_combis = []
            synergy_scores = []

            for _, row in df.iterrows():
                result = self._detect_row_synergies(row, available_columns)
                synergy_combis.append(result["synergy_combi"])
                synergy_scores.append(result["synergy_score"])

            df = df.copy()
            df["synergy_combi"] = synergy_combis
            df["synergy_combi_score"] = synergy_scores
            return df

        except Exception as e:
            raise SynergyDetectionError(f"Dataframe synergy detection failed: {e}") from e


    def _detect_row_synergies(self, row: pd.Series, feature_columns: List[str]) -> Dict[str, Any]:
        """
        Detect synergies for a single row using only textual pattern matching.
        """
        detected_text = " ".join(
            str(row[col]).lower() for col in feature_columns if pd.notna(row[col])
        )

        found_cats = set()
        for cat1, cat2 in self.priority_synergies.keys():
            if cat1.split("_")[0] in detected_text:
                found_cats.add(cat1)
            if cat2.split("_")[0] in detected_text:
                found_cats.add(cat2)

        all_synergies = []
        for (cat1, cat2), info in self.priority_synergies.items():
            if cat1 in found_cats and cat2 in found_cats:
                synergy_name = f"{cat1.split('_')[0]}-{cat2.split('_')[0]}"
                synergy_dict = {synergy_name: [cat1, cat2]}
                all_synergies.append({"dict": synergy_dict, "score": info["score"]})

        if all_synergies:
            all_synergies.sort(key=lambda x: x["score"], reverse=True)
            synergy_combi = {k: v for s in all_synergies for k, v in s["dict"].items()}
            max_score = all_synergies[0]["score"]
        else:
            synergy_combi = None
            max_score = 0.0

        return {"synergy_combi": synergy_combi, "synergy_score": max_score}
   


    ''' #first option claude
    #===========================================================
    def detect_synergies_from_dataframe(self, df: pd.DataFrame,
        corrosion_synergy_dict: Dict[str, List[str]],
        feature_columns: List[str]
        ) -> pd.DataFrame:
        """Detect synergies for each row in a dataframe based on feature columns.
        Args:df: Input dataframe with feature columns
            corrosion_synergy_dict: Dictionary mapping functional categories to terms
            feature_columns: List of column names to check for synergy terms
        Returns: DataFrame with added 'synergy_combi' and 'synergy_combi_score' columns
        Raises: SynergyDetectionError: If synergy detection fails
        """
        try:
            # Filter to only existing columns
            available_columns = [col for col in feature_columns if col in df.columns]
            
            # Process each row and collect synergy info
            synergy_combis = []
            synergy_scores = []
            
            for idx, row in df.iterrows():
                row_synergy = self._detect_row_synergies(
                    row, 
                    available_columns, 
                    corrosion_synergy_dict
                )
                
                synergy_combis.append(row_synergy['synergy_combi'])
                synergy_scores.append(row_synergy['synergy_score'])
            
            # Add columns to dataframe
            return pd.DataFrame({
                'synergy_combi': synergy_combis,
                'synergy_combi_score': synergy_scores
            }, index=df.index)
            
        except Exception as e:
            raise SynergyDetectionError(f"Dataframe synergy detection failed: {e}") from e

    def _detect_row_synergies(self, row: pd.Series,feature_columns: List[str],
        corrosion_synergy_dict: Dict[str, List[str]]
        ) -> Dict[str, Any]:
        """ Detect synergies for a single row based on feature columns.
        Uses ONLY priority_synergies for detection.        
        Args:row: Single row from dataframe
            feature_columns: List of columns to check
            corrosion_synergy_dict: Dictionary mapping categories to terms
        Returns:Dictionary with 'synergy_combi' and 'synergy_score'"""
        # Collect all terms from feature columns
        all_terms = set()
        for col in feature_columns:
            value = row[col]
            if pd.notna(value):
                if isinstance(value, (list, set)):
                    all_terms.update(str(v).lower() for v in value)
                elif isinstance(value, str):
                    # Handle comma-separated or pipe-separated values
                    terms = [t.strip().lower() for t in value.replace('|', ',').split(',')]
                    all_terms.update(terms)
        
        # Map terms to functional categories using corrosion_synergy_dict
        detected_categories = {}
        for category, category_terms in corrosion_synergy_dict.items():
            matched_terms = set()
            for term in category_terms:
                term_lower = str(term).lower()
                if term_lower in all_terms:
                    matched_terms.add(term)
            
            if matched_terms:
                detected_categories[category] = list(matched_terms)
        
        # Name mapping for short synergy names
        name_map = {'iron_metabolism': 'Fe',
            'sulfur_metabolism': 'S',
            'organic_acid_metabolism': 'organic_acid',
            'metal_binding_chelation': 'chelation',
            'biofilm_formation': 'biofilm',
            'o2_consumption': 'O2',
            'nitrogen_metabolism': 'N',
            'h2_consumption': 'H2'
        }
        # CHECK FOR ALL MATCHES IN priority_synergies

        all_synergies = []

        for synergy_pair, synergy_info in self.priority_synergies.items():
            cat1, cat2 = synergy_pair
            if cat1 in detected_categories and cat2 in detected_categories:
                current_score = synergy_info['score']
                short_name1 = name_map.get(cat1, cat1)
                short_name2 = name_map.get(cat2, cat2)
                synergy_name = f"{short_name1}-{short_name2}"
                synergy_terms = detected_categories[cat1] + detected_categories[cat2]
                
                all_synergies.append({
                    'dict': {synergy_name: synergy_terms},
                    'score': current_score
                })

        # Sort by score descending
        all_synergies.sort(key=lambda x: x['score'], reverse=True)

        # Extract just the dicts and max score
        synergy_dicts = [s['dict'] for s in all_synergies]
        max_score = all_synergies[0]['score'] if all_synergies else 0.0

        return {
            'synergy_combi': synergy_dicts if synergy_dicts else None,
            'synergy_score': max_score
        }
        #===========================================================
        # Check ONLY priority_synergies for matches
         for synergy_pair, synergy_info in self.priority_synergies.items():
            cat1, cat2 = synergy_pair
            
            # Check if both categories are present
            if cat1 in detected_categories and cat2 in detected_categories:
                current_score = synergy_info['score']
                
                if current_score > max_score:
                    max_score = current_score
                    
                    # Create short synergy name
                    short_name1 = name_map.get(cat1, cat1)
                    short_name2 = name_map.get(cat2, cat2)
                    best_synergy_name = f"{short_name1}-{short_name2}"
                    
                    # Collect all terms from both categories
                    best_synergy_terms = detected_categories[cat1] + detected_categories[cat2]
        
        # Build result
        if max_score > 0.0:
            synergy_combi = {best_synergy_name: best_synergy_terms}
        else:
            synergy_combi = None
        
        return {
            'synergy_combi': synergy_combi,
            'synergy_score': max_score
        }'''
        #===========================================================

# Custom exception (to keep your code functional)
class SynergyDetectionError(Exception):
    pass