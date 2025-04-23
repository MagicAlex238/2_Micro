#====================================NEW VERSION===============================================================
#==================================================================================
#=================== IMPORTS =======================================================
'''!pip install fuzzywuzzy
!pip install lxml pandas
!pip install pyarrow
!pip install openpyxl
!pip install python-Levenshtein
!pip install -U kaleido
!pip install natsort'''

# Standard library imports
import os
import sys
import ast
import subprocess
import logging
import time
from datetime import datetime
import shutil
from io import StringIO
from pathlib import Path
import re
from IPython import get_ipython
from IPython.display import display

# Data processing and analysis
import pandas as pd
import numpy as np
import openpyxl
import Levenshtein
from scipy.signal import savgol_filter
from community import community_louvain
from joblib import Parallel, delayed

import xml.etree.ElementTree as ET
# Retrieval and requesting 
from lxml import etree
import requests

# Utility libraries
import gzip
import random
from natsort import natsorted
from typing import Dict, List, Tuple, Set, Optional
import pickle
import gc
import joblib
import os
import csv
import json
import pyarrow.parquet as pq

os.environ['DISPLAY'] = ':0'

# large galaxies input and output #large size dir for large files hosted instead in kaggle
large_dir = Path("/home/beatriz/MIC")
output_large = large_dir / "output_large"

ECcontri_Uniprot_path = output_large/ "ECcontri_Uniprot.parquet"
ECcontri_Uniprot = pd.read_parquet(ECcontri_Uniprot_path, engine="pyarrow")

# Load from JSON with timing
json_path = output_large/ "ec_records_flat.json"
print(f"Starting JSON load from {json_path}...")
with open(json_path, 'r') as f:
    ec_records_flat = json.load(f)

#====================================NEW VERSION===============================================================
def enrich_eccontri_data(eccontri_df, ec_records_flat):
    """
        import sys
    #sys.path.append('/kaggle/input/corrosion-scoring') 
    sys.path.append('/home/beatriz/MIC/2_Micro/corrosion_scoring')
    from corrosion_scoring.global_terms import (
    #from global_terms import (
        metal_terms, 
    Enrich the ECcontri_Uniprot dataframe with complete information from ec_records dictionary
    It does search protein names, enzyme names and EC numbers for metadata to enrich the orginal eccontri
    then when it founds an Uncharacterised protein, it replaces this name with the enzyme names corresponding by the EC number
    Parameters:    eccontri_df : pandas DataFrame original ECcontri_Uniprot data with EC numbers in format x.x.x.x
                   ec_records : list of Dictionary where keys are EC numbers (without 'EC:' prefix) and values are metadata dictionaries

    Returns:       enriched_df : pandas DataFrame with original OTUs contribution abundances + additional metadata columns
    """
    import sys
    #sys.path.append('/kaggle/input/corrosion-scoring')  
    sys.path.append('/home/beatriz/MIC/2_Micro/corrosion_scoring')

    try:
        #from global_terms import (
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
        import scoring_system as score_sys
        using_imported_modules = True

        print("Successfully imported scoring modules")
    except ImportError as e:
        print(f"Warning: Could not import scoring modules: {e}")
        using_imported_modules = False
    
    # Make a copy to avoid modifying the original
    enriched_df = eccontri_df.copy()

    start_time = time.time()

    # 1. Create dictionaries for faster lookups
    print("Creating lookup dictionaries...")

    # EC number dictionary flaten
    ec_dict = {record['ec_number']: record for record in ec_records_flat if 'ec_number' in record}
    print(f"Created EC dictionary with {len(ec_dict)} entries")

    # Protein name dictionary (handles both string and list formats)
    protein_name_dict = {}
    for record in ec_records_flat:
        enzyme_names = record.get('enzyme_names', [])  # Default to empty list
        if isinstance(enzyme_names, str):
            # Handle semicolon-separated strings
            names = [n.strip().lower() for n in enzyme_names.split(';') if n.strip()]
        elif isinstance(enzyme_names, list):
            # Handle already-split lists
            names = [n.strip().lower() for n in enzyme_names if n.strip()]
        else:
            names = []
        
        for name in names:
            protein_name_dict[name] = record
    
    print(f"Created protein name dictionary with {len(protein_name_dict)} entries")

    # 2. Create a mapping dictionary to store all matches
    idx_to_metadata = {}

    # 3. Initialise Metadata  columns
    metadata_columns = ['enzyme_names', 'enzyme_class', 'pathways', 'hierarchy', 'metals_consolidated', 'corrosion_mechanisms', 'corrosion_relevance_score', 
                        'corrosion_relevance', 'metal_scores', 'overall_metal_score',  'corrosion_mechanism_scores', 'overall_corrosion_score', 
                        'functional_categories', 'overall_functional_score', 'corrosion_keyword_groups',  
                        'corrosion_keyword_scores', 'overall_keyword_score', 'corrosion_synergies',
                        'pathway_categories', 'organic_processes', 'overall_organic_process_score']
    for col in metadata_columns:
        enriched_df[col] = None

    # a boolean 'has_metal' column
    enriched_df['has_metal'] = False
    #=========================================================================================
    # Only proceed if we have metadata (either from EC or from protein/enzyme name)
    def _apply_metadata(df, idx, metadata, row):
        """Applies metadata to a specific row in the DataFrame."""

        """Applies metadata to a specific row in the DataFrame."""
        if metadata:
            # Protein Name Update
            if pd.isna(df.at[idx, 'protein_name']) or df.at[idx, 'protein_name'].lower() == "uncharacterized protein":
                enzyme_names = metadata.get('enzyme_names')
                if enzyme_names:
                    if isinstance(enzyme_names, list):
                        row['protein_name'] = '; '.join(enzyme_names)
                    else:
                        row['protein_name'] = str(enzyme_names)
            # Basic Metadata
            basic_fields = ['enzyme_names', 'enzyme_class', 'pathways', 'pathway_categories', 'hierarchy', 'corrosion_mechanisms',
                            'functional_categories', 'corrosion_keyword_groups', 'corrosion_synergies', 'organic_processes']
            for field in basic_fields:
                value = metadata.get(field)
                if value:
                    if isinstance(value, list):
                        row[field] = '; '.join(map(str, value))
                    else:
                        row[field] = str(value)
            # Corrosion Relevance
            if 'corrosion_relevance' in metadata:
                row['corrosion_relevance'] = metadata['corrosion_relevance']
            # Metals Consolidated + Metal Flag
            metals = metadata.get('metals_consolidated')
            if metals:
                if isinstance(metals, list):
                    row['metals_consolidated'] = '; '.join(metals)
                else:
                    row['metals_consolidated'] = str(metals)

            # Score Fields
            score_fields = ['overall_metal_score', 'overall_corrosion_score', 'overall_functional_score',
                            'overall_keyword_score', 'overall_organic_process_score', 'overall_synergy_score',
                            'corrosion_relevance_score']
            for field in score_fields:
                value = metadata.get(field)
                try:
                    row[field] = float(value) if value is not None else None
                except (ValueError, TypeError):
                    print(f"Row {idx}: Could not convert {field} to float")
                    row[field] = None
            # Dictionary Fields
            dict_fields = ['metal_scores', 'corrosion_mechanism_scores', 'corrosion_keyword_scores',
                        'organic_process_scores', 'corrosion_synergy_scores']
            for field in dict_fields:
                value = metadata.get(field)
                row[field] = json.dumps(value) if value else None  #Handle null values
        return row  
    #===========================================================================================

    # 4. Helper function and Map-based enrichement to avoid DF iteration
    def _enrich_row(row):
        protein_name = row['protein_name']
        ec_num = row['EC']
        idx = row.name
        metadata = None  # Initialize metadata

        if protein_name and protein_name.lower() != "uncharacterized protein":
            protein_name = protein_name.lower()
            metadata = protein_name_dict.get(protein_name)
            if not metadata:  # Try partial match if no direct match
                for name, record in protein_name_dict.items():
                    if name in protein_name or protein_name in name:
                        metadata = record
                        break
        else:
            if ec_num:  # If protein name is missing, try EC-only matching
                metadata = ec_dict.get(ec_num)

        # If metadata found, apply it; otherwise try to generate it with scoring system if available
        if metadata:
            row = _apply_metadata(eccontri_df, idx, metadata, row)
            eccontri_df.loc[idx] = row  # Update row in eccontri_df
        elif using_imported_modules and protein_name:
            # If no metadata found but scoring modules available, try to generate scores
            # This is useful for proteins not in ec_records but still in our dataframe
            try:
                # Generate text to analyze
                analysis_text = protein_name
                if not pd.isna(row.get('enzyme_class')):
                    analysis_text += ' ' + row['enzyme_class']
                
                # Calculate scores using the scoring system
                score_results = score_sys.calculate_overall_scores(analysis_text)
                
                # Apply the generated scores to the row
                for key, value in score_results.items():
                    if key in row:
                        if isinstance(value, list):
                            row[key] = '; '.join(map(str, value))
                        elif isinstance(value, dict):
                            row[key] = json.dumps(value)
                        else:
                            row[key] = value
                
                # Calculate corrosion relevance
                if all(k in score_results for k in ['overall_metal_score', 'overall_corrosion_score', 'overall_organic_process_score', 'overall_keyword_score']):
                    corrosion_relevance_score, corrosion_relevance = score_sys.calculate_corrosion_relevance_score(
                        score_results.get('overall_metal_score', 0),
                        score_results.get('overall_corrosion_score', 0),
                        score_results.get('overall_organic_process_score', 0),
                        score_results.get('overall_keyword_score', 0),
                        score_results.get('overall_synergy_score', 0),
                        score_results.get('overall_functional_score', 0)
                    )
                    row['corrosion_relevance_score'] = corrosion_relevance_score
                    row['corrosion_relevance'] = corrosion_relevance
                
                # Update row in dataframe
                eccontri_df.loc[idx] = row
            except Exception as e:
                print(f"Error generating scores for row {idx}: {e}")
           
        return row
    # Apply function over each row
    eccontri_df = eccontri_df.apply(_enrich_row, axis=1)

    #==========================================================================================================================
     # 7. Row Processing Setup # Define progress reporting
    total_rows = len(enriched_df)
    log_interval = max(1, min(10000, total_rows // 20))  # Log at most 20 times, minimum every 10000 rows

    print(f"Processing {total_rows} rows with logging every {log_interval} rows")

    # 8. Report and Return
    total_time = time.time() - start_time
    total_rows = len(eccontri_df)

    print(f"Completed enrichment in {total_time:.2f} seconds")
    print(f"Processed {total_rows} rows at {total_rows/total_time:.1f} rows/second")

    metadata_counts = {col: eccontri_df[col].notnull().sum() for col in metadata_columns}
    print("\nMetadata population statistics:")
    for col, count in metadata_counts.items():
        print(f"  {col}: {count} rows ({count/total_rows*100:.1f}%)")

    return eccontri_df
#========================================================================================================================+


if __name__ == "__main__":
    pre_ECcontri_Uniprot_enriched= enrich_eccontri_data(ECcontri_Uniprot, ec_records_flat)
    
    pre_path = output_large / "pre_ECcontri_Uniprot_enriched_v2.parquet"
    pre_ECcontri_Uniprot_enriched.to_parquet(pre_path)

    print(f"Saved {pre_path} with {len(pre_ECcontri_Uniprot_enriched)} rows")
    print(f"File exists: {pre_path.exists()} | Size: {pre_path.stat().st_size} bytes")