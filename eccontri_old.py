#====================================OLD VERSION===================================
#=================== IMPORTS =======================================================
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

# Directory to output large files # eccontris, compilated db
output_large = large_dir / "output_large"

ECcontri_Uniprot_path = output_large/ "ECcontri_Uniprot.parquet"
ECcontri_Uniprot = pd.read_parquet(ECcontri_Uniprot_path)

# Load from JSON with timing
json_path = output_large/ "ec_records_flat.json"
print(f"Starting JSON load from {json_path}...")
start_time = time.time()

with open(json_path, 'r') as f:
    ec_records_flat = json.load(f)

end_time = time.time()
elapsed = end_time - start_time
print(f"Loaded {len(ec_records_flat)} records from JSON in {elapsed:.2f} seconds")

#====================================OLD VERSION===============================================================

def enrich_eccontri_data(eccontri_df, ec_records_flat):
    """
    Enrich the ECcontri_Uniprot dataframe with complete information from ec_records dictionary
    It does search protein names, enzyme names and EC numbers for metadata to enrich the orginal eccontri
    then when it founds an Uncharacterised protein, it replaces this name with the enzyme names corresponding by the EC number
    Parameters:    eccontri_df : pandas DataFrame original ECcontri_Uniprot data with EC numbers in format EC:x.x.x.x
                   ec_records : list of Dictionary where keys are EC numbers (without 'EC:' prefix) and values are metadata dictionaries

    Returns:       enriched_df : pandas DataFrame with additional metadata columns
    """
    # Make a copy to avoid modifying the original
    enriched_df = eccontri_df.copy()

    start_time = time.time()

    # Create dictionaries for faster lookups
    print("Creating lookup dictionaries...")

    # EC number dictionary flaten
    ec_dict = {record['ec_number']: record for record in ec_records_flat if 'ec_number' in record}
    print(f"Created EC dictionary with {len(ec_dict)} entries")

    # Protein name dictionary
    protein_name_dict = {}
    for record in ec_records_flat:
        enzyme_names_str = record.get('enzyme_names', '') # dont understand here why create the dict for proteins directly from enzyme_names?
        if enzyme_names_str:
            # Flattened records will have enzyme_names as semicolon-separated string
            names = [name.strip() for name in enzyme_names_str.split(';')]
            for name in names:
                if name:  # Skip empty names
                    protein_name_dict[name.lower()] = record
    print(f"Created protein name dictionary with {len(protein_name_dict)} entries")

    # Create a mapping dictionary to store all matches
    idx_to_metadata = {}

    # Add all metadata columns
    metadata_columns = ['enzyme_names', 'enzyme_class', 'pathways', 'hierarchy', 'metals_consolidated', 'corrosion_mechanisms', 'corrosion_relevance_score', 
                        'corrosion_relevance', 'metal_scores', 'overall_metal_score',  'corrosion_mechanism_scores', 'overall_corrosion_score', 
                        'functional_categories', 'overall_functional_score', 'corrosion_keyword_groups',  
                        'corrosion_keyword_scores', 'overall_keyword_score', 'corrosion_synergies',
                        'pathway_categories', 'organic_processes', 'overall_organic_process_score']
    for col in metadata_columns:
        enriched_df[col] = None

    # a boolean 'has_metal' column
    enriched_df['has_metal'] = False

     # Define progress reporting
    total_rows = len(enriched_df)
    log_interval = max(1, min(10000, total_rows // 20))  # Log at most 20 times, minimum every 10000 rows

    print(f"Processing {total_rows} rows with logging every {log_interval} rows")

    # try protein name matches
    print("Performing protein name matches...")
    # Get rows without EC+Genus matches
    remaining_indices = set(enriched_df.index) - set(idx_to_metadata.keys())
    mask_remaining = enriched_df.index.isin(remaining_indices)
    mask_valid_protein = enriched_df['protein_name'].notna() & (enriched_df['protein_name'] != "Uncharacterized protein")
    mask_protein_match = mask_remaining & mask_valid_protein

    # This part still needs row-by-row processing for fuzzy matching
    protein_matches = 0
    for idx in enriched_df.index[mask_protein_match]:
        protein_name = enriched_df.loc[idx, 'protein_name'].lower()

        # Direct lookup in protein name dictionary
        if protein_name in protein_name_dict:
            idx_to_metadata[idx] = protein_name_dict[protein_name]
            protein_matches += 1
        else:
            # Try partial matches
            for name, record in protein_name_dict.items():
                if name in protein_name or protein_name in name:
                    idx_to_metadata[idx] = record
                    protein_matches += 1
                    break

    print(f"Found {protein_matches} protein name matches")

    # For any remaining rows, try EC-only matching
    print("Performing EC-only matches...")
    # Get rows without matches so far
    remaining_indices = set(enriched_df.index) - set(idx_to_metadata.keys())
    mask_remaining = enriched_df.index.isin(remaining_indices)
    mask_valid_ec = enriched_df['EC'].notna()
    mask_ec_match = mask_remaining & mask_valid_ec

    ec_only_matches = 0
    for idx in enriched_df.index[mask_ec_match]:
        ec_num = enriched_df.loc[idx, 'EC']
        if ec_num in ec_dict:
            idx_to_metadata[idx] = ec_dict[ec_num]
            ec_only_matches += 1

    print(f"Found {ec_only_matches} EC-only matches")

    # Apply all metadata in one go based on the matches we found
    print("Applying metadata to matched rows...")
    for idx, metadata in idx_to_metadata.items():

        # Only proceed if we have metadata (either from EC or from protein/enzyme name)
        if metadata is not None:

            # If protein_name is missing or uncharacterized, try to fill it with enzyme_names
            if pd.isna(enriched_df.at[idx, 'protein_name']) or enriched_df.at[idx, 'protein_name'].lower() == "Uncharacterized protein":
                if 'enzyme_names' in metadata and metadata['enzyme_names']:
                    # Handle both flattened (string) and non-flattened (list) records
                    if isinstance(metadata['enzyme_names'], list):
                        enriched_df.at[idx, 'protein_name'] = '; '.join(metadata['enzyme_names'])
                    else:
                        enriched_df.at[idx, 'protein_name'] = str(metadata['enzyme_names'])

            # Add basic metadata
            for field in ['enzyme_names', 'enzyme_class', 'pathways','pathway_categories', 'hierarchy', 'corrosion_mechanisms',
                          'functional_categories',  'corrosion_keyword_groups', 'corrosion_synergies', 'organic_processes']:
                if field in metadata and metadata[field]:
                    if isinstance(metadata[field], list):
                        enriched_df.at[idx, field] = '; '.join(str(v) for v in metadata[field])
                    else:
                        enriched_df.at[idx, field] = str(metadata[field])
            
            # add corrosion_relevance category:
            if 'corrosion_relevance' in metadata:
                enriched_df.at[idx, 'corrosion_relevance'] = metadata['corrosion_relevance']

            # add the consolidated metals field:
            if 'metals_consolidated' in metadata and metadata['metals_consolidated']:
                if isinstance(metadata['metals_consolidated'], list):
                    enriched_df.at[idx, 'metals_consolidated'] = '; '.join(metadata['metals_consolidated'])
                    # Set has_metal flag based on whether metals exist
                    enriched_df.at[idx, 'has_metal'] = len(metadata['metals_consolidated']) > 0
                else:
                    enriched_df.at[idx, 'metals_consolidated'] = str(metadata['metals_consolidated'])
                    # Set has_metal flag for string case too
                    enriched_df.at[idx, 'has_metal'] = bool(metadata['metals_consolidated'])

            # Adding the scores for each component , simple score fields, handle in a loop
            score_fields = ['overall_metal_score', 'overall_corrosion_score', 'overall_functional_score', 
                            'overall_keyword_score', 'overall_organic_process_score', 'overall_synergy_score',
                            'corrosion_relevance_score']
            
            for field in score_fields:
                if field in metadata:
                    try:
                        enriched_df.at[idx, field] = float(metadata.get(field, 0))
                    except (ValueError, TypeError):
                        print(f"Row {idx}: Could not convert {field} to float")
                        enriched_df.at[idx, field] = None
            
            # For dictionary score fields, handle in a loop too
            dict_fields = ['metal_scores', 'corrosion_mechanism_scores', 'corrosion_keyword_scores', 
                           'organic_process_scores', 'corrosion_synergy_scores']
            
            for field in dict_fields:
                if field in metadata:
                    enriched_df.at[idx, field] = json.dumps(metadata.get(field, {}))
                               
    # Final report
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Completed enrichment in {total_time:.2f} seconds")
    print(f"Processed {total_rows} rows at {total_rows/total_time:.1f} rows/second")

    # Count non-null values in the metadata columns to see success rate
    metadata_counts = {col: enriched_df[col].notnull().sum() for col in metadata_columns}
    print("\nMetadata population statistics:")
    for col, count in metadata_counts.items():
        print(f"  {col}: {count} rows ({count/total_rows*100:.1f}%)")

    return enriched_df

#========================================================================

if __name__ == "__main__":
    pre_ECcontri_Uniprot_enriched= enrich_eccontri_data(ECcontri_Uniprot, ec_records_flat)
    
    pre_path = output_large / "pre_ECcontri_Uniprot_enriched.parquet"
    pre_ECcontri_Uniprot_enriched.to_parquet(pre_path)

    print(f"Saved {pre_path} with {len(pre_ECcontri_Uniprot_enriched)} rows")
    print(f"File exists: {pre_path.exists()} | Size: {pre_path.stat().st_size} bytes")

    elapsed_time_main = time.time() - start_time_main
    print(f"Processing took {elapsed_time_main:.2f} seconds")
