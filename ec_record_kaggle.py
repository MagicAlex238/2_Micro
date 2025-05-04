 #==================================RUN METABOLISM DB =============================
#============ IMPORTS AND INSTALLS =======================
# Making sure to use same python version for compatibility
'''!sudo apt-get update -y
!sudo apt-get install python3.10
!sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1
!python --version
!pip install fuzzywuzzy
!pip install lxml pandas
!pip install pyarrow
!pip install python-Levenshtein
!pip install openpyxl
!pip install adjustText
!pip install git+https://github.com/MagicAlex238/2_Micro.git'''
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

#from joblib import Parallel, delayed
import xml.etree.ElementTree as ET
# Retrieval and requesting 
from lxml import etree
import requests

# Utility libraries
import gzip
import random
import pickle
import gc
import joblib
import os
import csv
import json
import pyarrow.parquet as pq
# Own Scoring system
import corrosion_scoring as cs

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
# ===========================## 9.1 Setting up Paths and Parsing the Dataframes ==================

#db_dir = Path("/kaggle/input/databases/Databases")
db_dir = Path("/home/beatriz/MIC/Databases")

def setup_paths():
    """Set up paths for database access"""

    # Database paths
    db_paths = {
        'enzyme': db_dir / 'enzyme',
        'enzyme_class': db_dir / 'enzclass.txt',
        'enzyme_brenda' : db_dir/ 'brenda_2024.txt',
        'ko': db_dir / 'ko',
        'ko_hierarchy': db_dir / 'ko_hierarchy.txt',
        'pathway': db_dir / 'pathway',
        'module': db_dir / 'module',
        'reaction': db_dir / 'reaction',
        'compound': db_dir / 'compound',
        'metalpdb': db_dir / 'flat_db_file.xml',
        'ko_pathway': db_dir / 'ec_pathway.list'
    }

    return db_paths

#  Calling the paths
if __name__ == "__main__":
    paths = setup_paths()
    # Print paths to verify
    for db_name, path in paths.items():
        print(f"{db_name}: {path}")
        print(f"Exists: {path.exists()}")

#====================== Clean Names ====================
from fuzzywuzzy import fuzz

def enhanced_clean_protein_name(name):
    """
    Enhanced protein name cleaning that builds on the original function
    with additional normalization steps:

    1. Handles EC numbers appropriately
    2. Removes parenthetical content and brackets
    3. Removes duplicated terms and standardizes spacing
    4. Handles specific protein families consistently
    5. Standardizes capitalizations and hyphenations
    6. Removes common suffixes that don't change the protein identity
    7. Normalizes common protein name patterns
    """
    if pd.isna(name):
        return "Uncharacterized protein"

    # Save EC numbers for later reattachment if needed
    ec_match = re.search(r'(EC\s*[\d\.]+)', name)
    ec_number = ec_match.group(1) if ec_match else None

    # If the name is just an EC number in any format, return it
    if re.match(r'^[\s]*EC\s*[\d\.]+[\s]*$', name):
        return name.strip()

    # Remove EC numbers and content in parentheses
    name = re.sub(r'EC\s∗[\d\.]+EC\s*[\d\.]+', '', name)
    name = re.sub(r'[)]∗[^)]*', '', name)
    name = re.sub(r'[^[^]*\]', '', name)  # Also remove content in square brackets

    # Convert to lowercase for better matching
    name = name.lower()

    # Standardize spacing around hyphens, slashes
    name = re.sub(r'[\s]*[\-\/][\s]*', '-', name)

    # Remove common suffixes that don't change protein identity
    name = re.sub(r'\s+(protein|domain|enzyme|family|subunit|chain|component|type|complex|fragment|precursor)$', '', name)

    # Remove specific protein ID suffixes (like FabG)
    name = re.sub(r'\s+[a-z]+\d+$', '', name)  # e.g., "reductase FabG" -> "reductase"

    # Normalize common protein name patterns
    replacements = {
        # Standardize dehydrogenases
        r'(\w+)\s*dehydrogenase': r'\1-dehydrogenase',

        # Standardize reductases
        r'(\w+)\s*reductase': r'\1-reductase',
        r'\[acyl-carrier-protein\]': 'acp',
        r'\[acyl carrier protein\]': 'acp',
        r'3-oxoacyl-acp reductase': '3-oxoacyl-acp-reductase',

        # Standardize synthetases
        r'(\w+)\s*synthase': r'\1-synthase',
        r'(\w+)\s*synthetase': r'\1-synthetase',

        # Standardize common protein terms
        r'alcohol dehydrogenase': 'alcohol-dehydrogenase',
        r'glutathione dehydrogenase': 'glutathione-dehydrogenase',
        r's-glutathione dehydrogenase': 'glutathione-dehydrogenase',
        r'threonine dehydrogenase': 'threonine-dehydrogenase',
        r'l-threonine dehydrogenase': 'threonine-dehydrogenase'
    }

    for pattern, replacement in replacements.items():
        name = re.sub(pattern, replacement, name)

    # Split into words and remove duplicates while preserving order
    words = name.split()
    seen = set()
    unique_words = []
    for word in words:
        if word not in seen:
            seen.add(word)
            unique_words.append(word)

    # Rejoin words
    name = ' '.join(unique_words)

    # Remove specific redundant patterns
    redundant_patterns = [
        (r'enzyme\s+enzyme', 'enzyme'),
        (r'synthase\s+synthase', 'synthase'),
        (r'transferase\s+transferase', 'transferase'),
        (r'-glucan\s+glucan', 'glucan'),
        (r'protein\s+protein', 'protein'),
        (r'reductase\s+reductase', 'reductase'),
        (r'dehydrogenase\s+dehydrogenase', 'dehydrogenase')
    ]

    for pattern, replacement in redundant_patterns:
        name = re.sub(pattern, replacement, name, flags=re.IGNORECASE)

    # Reattach EC number if it was the main identifier
    if ec_number and len(name.strip()) < 5:  # If the remaining name is very short
        name = f"{name} {ec_number}" if name.strip() else ec_number

    return name.strip()

def normalize_dataset(df, name_col='protein_name', sample_size=20):
    """
    Apply enhanced protein name cleaning to a dataset and
    show before/after examples
    """
    # Create a copy with normalized names
    normalized_df = df.copy()
    normalized_df['normalized_protein'] = normalized_df[name_col].apply(enhanced_clean_protein_name)

    # Display samples of the normalization
    samples = df[[name_col]].drop_duplicates().sample(min(sample_size, df[name_col].nunique()))

    print(f"Original unique protein names: {df[name_col].nunique()}")
    samples['normalized'] = samples[name_col].apply(enhanced_clean_protein_name)

    print("\nSample normalization results:")
    for _, row in samples.iterrows():
        print(f"\nOriginal:   {row[name_col]}")
        print(f"Normalized: {row['normalized']}")

    # Count reduction in unique names
    original_count = df[name_col].nunique()
    normalized_count = normalized_df['normalized_protein'].nunique()
    reduction = original_count - normalized_count
    reduction_pct = 100 * reduction / original_count if original_count > 0 else 0

    print(f"\nNormalized unique protein names: {normalized_count}")
    print(f"Reduction: {reduction} ({reduction_pct:.1f}%)")

    # Find similar protein names that are now treated as the same
    if normalized_count < original_count:
        print("\nExamples of proteins that were normalized to the same name:")
        norm_to_orig = {}

        # Build mapping of normalized names to original names
        for _, row in samples.iterrows():
            norm_name = row['normalized']
            orig_name = row[name_col]

            if norm_name not in norm_to_orig:
                norm_to_orig[norm_name] = []

            if orig_name not in norm_to_orig[norm_name]:
                norm_to_orig[norm_name].append(orig_name)

        # Display examples where multiple original names map to the same normalized name
        examples_shown = 0
        for norm_name, orig_names in norm_to_orig.items():
            if len(orig_names) > 1 and examples_shown < 5:
                print(f"\nNormalized to: {norm_name}")
                for orig in orig_names:
                    print(f"  - {orig}")
                examples_shown += 1

    return normalized_df

#==================BRENDA===========================
import logging
logging.basicConfig(level=logging.INFO)

def parse_brenda_file():
    """Parse BRENDA database file for detailed enzyme information"""
    paths = setup_paths()
    enzyme_brenda_path = paths['enzyme_brenda']

    ec_detailed_info = {}
    current_ec = None
    in_enzyme_entry = False

    try:
        with open(enzyme_brenda_path, 'r') as f:
          for line in f:
              line = line.strip()

              # Skip empty lines
              if not line:
                  continue

              # Check for the end of an entry
              if line == "///":
                  current_ec = None
                  in_enzyme_entry = False
                  continue

              # Process ID line - identify enzyme entries
              if line.startswith('ID\t'):
                  current_ec = line.split('\t')[1]

                  # Skip "spontaneous" and other non-EC entries
                  if not any(c.isdigit() for c in current_ec):
                      current_ec = None
                      in_enzyme_entry = False
                      continue

                  # Initialize proper EC entry
                  ec_detailed_info[current_ec] = {
                      'metals': [],
                      'cofactors': [],
                      'reactions': [],
                      'substrates': [],
                      'inhibitors': []
                  }
                  in_enzyme_entry = True

              # Only process other lines if we're in a valid enzyme entry
              elif in_enzyme_entry and current_ec:
                  if line.startswith('ME\t'):
                      # Extract metal information
                      metal_info = line.split('\t')[1]
                      ec_detailed_info[current_ec]['metals'].append(metal_info)

                  elif line.startswith('CF\t'):
                      # Extract cofactor information
                      cofactor_info = line.split('\t')[1]
                      ec_detailed_info[current_ec]['cofactors'].append(cofactor_info)

                  elif line.startswith('RE\t'):
                      # Extract detailed reaction information
                      reaction_info = line.split('\t')[1]
                      ec_detailed_info[current_ec]['reactions'].append(reaction_info)

                  elif line.startswith('SP\t') or line.startswith('NSP\t'):
                      # Extract substrate information
                      substrate_info = line.split('\t')[1]
                      ec_detailed_info[current_ec]['substrates'].append(substrate_info)

                  elif line.startswith('IN\t'):
                      # Extract inhibitor information
                      inhibitor_info = line.split('\t')[1]
                      ec_detailed_info[current_ec]['inhibitors'].append(inhibitor_info)

        # Verify we have valid EC numbers
        ec_detailed_info = {ec: info for ec, info in ec_detailed_info.items()
                            if ec.count('.') == 3 and all(part.isdigit() for part in ec.split('.'))}

    except Exception as e:
        logging.error("Error parsing BRENDA file: %s", e)
        return {}
    return ec_detailed_info

brenda_data = parse_brenda_file()
#brenda_data
#=============================================
def process_brenda_data(brenda_data):
    """Process BRENDA data to extract clean metal information while keeping other data intact"""
    
    processed_data = {}

    for ec_number, data in brenda_data.items():
        processed_data[ec_number] = {
            'cofactors': data.get('cofactors', []),
            'reactions': data.get('reactions', []),
            'substrates': data.get('substrates', []),
            'inhibitors': data.get('inhibitors', []),
            'raw_metals': data.get('metals', []),
            'clean_metals': [],
            'corrosion_mechanisms': [],
            'pathway_categories': {},
            'organic_processes': {},
            'corrosion_synergies': []
        }

        # Create a single text string for analysis
        all_text = ' '.join([
            ' '.join(data.get('reactions', [])),
            ' '.join(data.get('substrates', [])),
            ' '.join(data.get('cofactors', [])), 
            ' '.join(data.get('inhibitors', []))
        ]).lower()

        # Extract clean metal names first for both paths
        for entry in data.get('metals', []):
            entry_lower = entry.lower()
            for metal in metal_terms:
                metal_lower = metal.lower()
                if metal_lower in entry_lower:
                    if metal not in processed_data[ec_number]['clean_metals']:
                        processed_data[ec_number]['clean_metals'].append(metal)
        # Use the scoring system module for comprehensive scoring
        score_results = cs.calculate_overall_scores(
            all_text, 
            brenda_metals=processed_data[ec_number]['clean_metals'],
            pathways=' '.join(processed_data[ec_number].get('pathways', []))
        )
                
        # Update the processed data with scores from the module
        processed_data[ec_number]['corrosion_mechanisms'] = score_results['corrosion_mechanisms']
        processed_data[ec_number]['organic_processes'] = score_results['organic_processes']
        processed_data[ec_number]['corrosion_synergies'] = score_results['corrosion_synergies']
        
       # Add pathway categories - transform from scoring format to expected format
        for category_info in score_results.get('functional_categories', []):
            category = category_info.get('category')
            if category:
                processed_data[ec_number]['pathway_categories'][category] = True
        
        # Add metal information
        processed_data[ec_number]['metals_involved'] = score_results['metals_involved']
        
        # Define corrosion-relevant metals using cs.metal_terms 
        corrosion_relevant_metals = [metal for metal in cs.metal_terms if metal in cs.corrosion_metals]

        processed_data[ec_number]['corrosion_metals_from_brenda'] = [
            metal for metal in processed_data[ec_number]['clean_metals']
            if metal in corrosion_relevant_metals
        ]
        
        # Calculate pathway scores
        pathways = []  # No pathways directly from BRENDA, but include for API consistency
        enzyme_names = []  # No names from BRENDA, but include for API consistency
        enzyme_class = ""  # No class from BRENDA, but include for API consistency
        
        pathway_score, pathway_category_scores = cs.calculate_pathway_score(
            rec.get('pathways', []), enzyme_names, enzyme_class_text
)
        # Calculate synergy score
        synergy_score = cs.check_metal_organic_synergy(processed_data[ec_number]['clean_metals'], 
            []  # No enzyme names in BRENDA
        )
        
        # Calculate final corrosion relevance score using the module
        corrosion_relevance_score, corrosion_relevance = cs.calculate_corrosion_relevance_score(
            score_results.get('overall_metal_score', 0),
            score_results.get('overall_corrosion_score', 0),
            pathway_score,  score_results.get('overall_organic_process_score', 0),
            score_results.get('overall_keyword_score', 0),
            synergy_score, score_results.get('overall_functional_score', 0)
            )
            
        processed_data[ec_number]['corrosion_relevance_score'] = corrosion_relevance_score
        processed_data[ec_number]['corrosion_relevance'] = corrosion_relevance

        # Calculate comprehensive corrosion relevance score
        corrosion_score = 0
        
        # Score based on metal involvement
        corrosion_score += len(processed_data[ec_number]['corrosion_metals_from_brenda']) * 1.5
        
        # Score based on corrosion mechanisms
        corrosion_score += len(processed_data[ec_number]['corrosion_mechanisms']) * 2.0
        
        # Score based on pathway categories related to corrosion
        for pathway in cs.pathway_categories.keys():
            if pathway in processed_data[ec_number]['pathway_categories']:
                corrosion_score += 1.0
               
        # Score based on organic processes relevant to corrosion
        for process in cs.organic_categories.keys():
            if process in processed_data[ec_number]['organic_processes']:
                corrosion_score += cs.organic_categories[process].get('weight', 0.5)
                        
        # Score based on corrosion synergies
        corrosion_score += len(processed_data[ec_number]['corrosion_synergies']) * 2.0
        
        # Set corrosion relevance category
        if corrosion_score >= 5:
            processed_data[ec_number]['corrosion_relevance'] = 'high'
        elif corrosion_score >= 2:
            processed_data[ec_number]['corrosion_relevance'] = 'medium'
        else:
            processed_data[ec_number]['corrosion_relevance'] = 'low'
        
        # Store the numerical score
        processed_data[ec_number]['corrosion_relevance_score'] = corrosion_score

    return processed_data

#=================================================================
def read_enzyme_names():
    """Read and parse enzyme file to get EC numbers and their names"""
    paths = setup_paths()
    enzyme_path = paths['enzyme']

    ec_to_names = {}  # More descriptive name
    with open(enzyme_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                ec_number = parts[0]
                names = parts[1].split('; ')
                ec_to_names[ec_number] = names

    for ec_number, names_list in ec_to_names.items():

        cleaned_names = [enhanced_clean_protein_name(name) if isinstance(name, str) else str(name) for name in names_list]
        ec_to_names[ec_number] = cleaned_names

    return ec_to_names

#================================
def read_enzyme_class():
    paths = setup_paths()
    ec_file_path = paths['enzyme_class']

    enzyme_class = {}

    with open(ec_file_path, 'r') as f:
        for line in f:
            # Format is like "1. 1. 1.-    With NAD(+) or NADP(+) as acceptor."
            if line.strip() and any(line.startswith(str(i)) for i in range(1, 7)):
                parts = line.strip().split('  ')
                if len(parts) >= 2:
                    ec_id = parts[0].replace(' ', '')
                    desc = parts[1].strip()
                    enzyme_class[ec_id] = desc
    return enzyme_class
#==============================

def read_ec_pathway_mapping():
    """Read EC to pathway mapping file downloaded from KEGG"""
    paths = setup_paths()
    ko_pathway_path = paths['ko_pathway']

    ec_to_pathway = {}

    try:
        with open(ko_pathway_path , 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    ec_id = parts[0].replace('ec:', '')
                    pathway_id = parts[1].replace('path:', '')

                    if ec_id not in ec_to_pathway:
                        ec_to_pathway[ec_id] = []
                    ec_to_pathway[ec_id].append(pathway_id)

        print(f"Loaded pathway mappings for {len(ec_to_pathway)} EC numbers")
    except Exception as e:
        print(f"Error reading EC-pathway mapping: {e}")
        return {}

    return ec_to_pathway

#=======================================
def read_ko_data():
    """Read and parse KEGG KO file"""
    paths = setup_paths()
    ko_file_path = paths['ko']

    ko_info = {}
    with open(ko_file_path, 'r') as f:
        for line in f:
            if line.startswith('K'):
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    ko_info[parts[0]] = {
                        'definition': parts[1],
                        'pathway': parts[2] if len(parts) > 2 else ''
                    }
    return ko_info

#=======================================
def read_ko_hierarchy():
    paths = setup_paths()
    ko_path = paths['ko_hierarchy'] 

    hierarchy = {
        'A': {},  # Top level
        'B': {},  # ko Category
        'C': {},  # Pathway
        'D': {}   # KO/Enzyme
    }

    current = {'A': None, 'B': None, 'C': None}

    with open(ko_path, 'r') as f:
        for line in f:
            if line.startswith('A'):
                parts = line.strip().split()
                id = parts[1]
                name = ' '.join(parts[2:])
                hierarchy['A'][id] = name
                current['A'] = id

            elif line.startswith('B'):
                parts = line.strip().split()
                id = parts[1]
                name = ' '.join(parts[2:])
                hierarchy['B'][id] = {'name': name, 'parent': current['A']}
                current['B'] = id

            elif line.startswith('C'):
                parts = line.strip().split()
                id = parts[1]
                name = ' '.join(parts[2:])
                if '[PATH:' in name:
                    path_parts = name.split('[PATH:')
                    name = path_parts[0].strip()
                    path_id = path_parts[1].split(']')[0]
                else:
                    path_id = None

                hierarchy['C'][id] = {
                    'name': name,
                    'parent': current['B'],
                    'path_id': path_id
                }
                current['C'] = id

            elif line.startswith('D'):
                parts = line.strip().split()
                ko_id = parts[1]
                name = ' '.join(parts[2:])

                # Extract EC numbers if present
                ec_numbers = []
                if '[EC:' in name:
                    ec_part = name.split('[EC:')[1].split(']')[0]
                    ec_numbers = ec_part.split()
                    name = name.split('[EC:')[0].strip()

                hierarchy['D'][ko_id] = {
                    'name': name,
                    'parent': current['C'],
                    'ec_numbers': ec_numbers
                }

    return hierarchy

#=================================
def read_reaction_data():
    paths = setup_paths()
    reaction_file_path = paths['reaction']

    reaction_info = {}

    with open(reaction_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split(None, 1)  # Split on first whitespace
            if len(parts) >= 2:
                rxn_id = parts[0]
                desc_parts = parts[1].split(';')

                # First part is reaction name
                name = desc_parts[0].strip()

                # Rest might contain equation
                equation = desc_parts[1].strip() if len(desc_parts) > 1 else ""

                reaction_info[rxn_id] = {
                    'name': name,
                    'equation': equation
                }

    return reaction_info

#=======================
def read_pathway_data():
    paths = setup_paths()
    pathway_path = paths['pathway']

    pathway_info = {}
    with open(pathway_path, 'r') as f:
          for line in f:
              parts = line.strip().split('\t')
              if len(parts) >= 2:
                  pathway_id = parts[0]
                  pathway_name = parts[1]
                  pathway_info[pathway_id] = pathway_name
    return pathway_info

#=========================
def read_module_data():
    paths = setup_paths()
    module_path = paths['module']
    module_info = {}
    with open(module_path, 'r') as f:
      for line in f:
          parts = line.strip().split('\t')
          if len(parts) >= 2:
              module_id = parts[0]
              module_desc = parts[1]
              module_info[module_id] = module_desc
    return module_info

#==============================
def read_compound_data():
    paths = setup_paths()
    compound_path = paths['compound']
   
    compound_info = {}
    with open(compound_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    compound_id = parts[0]
                    compound_names = parts[1].split('; ')
                    compound_info[compound_id] = {
                        'name': compound_names[0],
                        'synonyms': compound_names[1:] if len(compound_names) > 1 else []
                    }
    return compound_info

start_time_metal_pbd = time.time()
#===============================
def parse_metalpdb_xml():
    """Parse MetalPDB XML file to extract metal-binding information"""
    
    paths = setup_paths()
    metalpdb_path = paths['metalpdb']
    
    metal_binding_data = {}

    try:
        # Use a more tolerant parser
        parser = etree.XMLParser(recover=True)
        tree = etree.parse(metalpdb_path, parser)
        root = tree.getroot()

        # Process each site
        for site in root.findall('.//site'):
            # Extract site information
            site_name = site.findtext('site_name')
            pdb_code = site.findtext('pdb_code')
            site_nuclearity = site.findtext('site_nuclearity')

            # Process each metal in the site
            for metal in site.findall('.//metal'):
                metal_symbol = metal.findtext('periodic_symbol')
                metal_name = metal.findtext('periodic_name')
                coordination_number = metal.findtext('coordination_number')
                geometry = metal.findtext('geometry')

                # Process ligands
                ligands = []
                for ligand in metal.findall('.//ligand'):
                    residue_name = ligand.findtext('residue_name')
                    residue_num = ligand.findtext('residue_pdb_number')
                    chain = ligand.findtext('chain_letter')
                    binding_type = ligand.findtext('endo_exo')

                    # Process donor atoms
                    donors = []
                    for donor in ligand.findall('.//donor'):
                        distance = donor.findtext('distance')
                        atom_name = donor.findtext('atom_pdb_name')
                        atom_symbol = donor.findtext('atom_symbol')
                        interaction_type = donor.findtext('interaction_type')

                        donors.append({
                            'distance': distance,
                            'atom_name': atom_name,
                            'atom_symbol': atom_symbol,
                            'interaction_type': interaction_type
                        })

                    ligands.append({
                        'residue_name': residue_name,
                        'residue_number': residue_num,
                        'chain': chain,
                        'binding_type': binding_type,
                        'donors': donors
                    })

                # Get the protein/molecule information
                site_chains = []
                for chain in site.findall('.//site_chain'):
                    molecule_name = chain.findtext('molecule_name')
                    molecule_type = chain.findtext('molecule_type')
                    chain_letter = chain.findtext('letter')

                    site_chains.append({
                        'molecule_name': molecule_name,
                        'molecule_type': molecule_type,
                        'chain_letter': chain_letter
                    })

                # Create a unique key for this metal site
                metal_site_key = f"{pdb_code}_{site_name}_{metal_symbol}"

                # Store the data
                metal_binding_data[metal_site_key] = {
                    'pdb_code': pdb_code,
                    'site_name': site_name,
                    'site_nuclearity': site_nuclearity,
                    'metal': {
                        'symbol': metal_symbol,
                        'name': metal_name,
                        'coordination_number': coordination_number,
                        'geometry': geometry
                    },
                    'ligands': ligands,
                    'site_chains': site_chains
                }

    except Exception as e:
        logging.error("Error parsing MetalPDB XML: %s", e)
        return {}

    return metal_binding_data

#=================================================================================================

def extract_metal_coordination_patterns(metal_binding_data):
    """Extract metal coordination patterns from MetalPDB data"""

    # Track metal coordination patterns
    metal_coordination = {}
    metal_residue_binding = {}
    # Add a counter to track processing
    processed_sites = 0

    for site_key, site_data in metal_binding_data.items():
        try:
            metal_symbol = site_data.get('metal', {}).get('symbol')
            if not metal_symbol:
                continue
    
            # Track coordination environments
            coord_num = site_data.get('metal', {}).get('coordination_number')
            geometry = site_data.get('metal', {}).get('geometry')

            if coord_num and geometry:
                coord_key = f"{metal_symbol}_{coord_num}_{geometry}"
                metal_coordination[coord_key] = metal_coordination.get(coord_key, 0) + 1
        
              # Track metal-residue binding
            if metal_symbol not in metal_residue_binding:
                metal_residue_binding[metal_symbol] = {}
                
            for ligand in site_data.get('ligands', []):
                residue = ligand.get('residue_name')
                if residue:
                    metal_residue_binding[metal_symbol][residue] = metal_residue_binding[metal_symbol].get(residue, 0) + 1
            
            processed_sites += 1
        except Exception as e:
            print(f"Error processing site {site_key}: {e}")
    
    print(f"Processed {processed_sites} sites out of {len(metal_binding_data)}")
    
    if not metal_coordination:
        print("Warning: No coordination patterns were extracted")
    if not metal_residue_binding:
        print("Warning: No residue binding patterns were extracted")
        
    return {
        'coordination': metal_coordination,
        'residue_binding': metal_residue_binding
    }

#=====================================================================================================================

def create_ec_to_reaction_mapping():  
    ec_to_names = read_enzyme_names()
    reaction_info = read_reaction_data()
    
    # Build keyword-to-ECs reverse index
    keyword_to_ecs = {}
    for ec, names in ec_to_names.items():
        words = set()
        for name in names:
            for word in str(name).lower().split():
                if len(word) > 4:
                    words.add(word)
        for keyword in words:
            if keyword not in keyword_to_ecs:
                keyword_to_ecs[keyword] = set()
            keyword_to_ecs[keyword].add(ec)
    
    # Preprocess reaction names
    rxn_name_lower = {rxn_id: info['name'].lower() for rxn_id, info in reaction_info.items()}
    
    # Match reactions to ECs
    ec_to_reaction = {}
    for rxn_id, rxn_name in rxn_name_lower.items():
        matched_ecs = set()
        for keyword in keyword_to_ecs:
            if keyword in rxn_name: 
                matched_ecs.update(keyword_to_ecs[keyword])
        for ec in matched_ecs:
            if ec not in ec_to_reaction:
                ec_to_reaction[ec] = set()
            ec_to_reaction[ec].add(rxn_id)
    
    return {ec: list(rxns) for ec, rxns in ec_to_reaction.items()}  


#====================================================================================================================================
def consolidate_metal_terms(brenda_metals, text_detected_metals):
    """
    Consolidates metal names from BRENDA and text mining into standardized symbols.
    
    Parameters:   brenda_metals (list of str): Metals obtained from BRENDA data.
    text_detected_metals (list of str): Metals detected from text mining.
    
    Returns:  list: Consolidated list of unique, standardized metal symbols.
    """
    # Import metal_mapping from global terms module
   
    consolidated = set()
    all_metals = (brenda_metals or []) + (text_detected_metals or [])
    
    for metal in all_metals:
        metal_norm = metal.strip().lower()
        # Check if the normalized term matches any key in the standard mapping
        for key, symbol in cs.metal_mapping.items():
            if key in metal_norm:
                consolidated.add(symbol)
                break
        else:
            # If no mapping is found, add the original
            consolidated.add(metal.strip())
    return list(consolidated)

#================================ MAIN FUNCTION ======================================
def create_metabolism_database(sample_size=None):
    """
    Build a list of dictionaries, each representing a single EC record.
    
    This function builds a comprehensive database of enzyme records with detailed information
    about their potential relevance to corrosion processes. It incorporates data from various
    sources, including BRENDA, MetalPDB, KEGG, and others.
    
    The function can use either:
    1. The external scoring_system_py module (preferred for consistency and maintainability)
    2. The original inline scoring logic (as a fallback if the module import fails)
    
    Args:   sample_size (int, optional): If specified, random sample of EC numbers to process
        
    Returns:  list: List of dictionaries, each containing information about an enzyme
    """
    
    try:
        # Read all necessary files
        print("Loading data sources...")
        ec_to_names = read_enzyme_names() or {}
        enzyme_class = read_enzyme_class() or {}
        reaction_equation = read_reaction_data() or {}
        ko_ec = read_ko_data() or {}
        ko_hierarchy = read_ko_hierarchy() or {}
        pathway_data = read_pathway_data() or {}
        module_info = read_module_data() or {}
        compound_info = read_compound_data() or {}
        brenda_en = process_brenda_data(brenda_data) or {}
        metal_binding_data = parse_metalpdb_xml() or {}
        metal_patterns = extract_metal_coordination_patterns(metal_binding_data)
        ec_pathway_mapping = read_ec_pathway_mapping() or {}

        # If sample_size is specified, take a random sample
        if sample_size and isinstance(sample_size, int) and sample_size > 0:
            print(f"Running with sample of {sample_size} EC numbers")
            ec_keys = list(ec_to_names.keys())
            if len(ec_keys) > sample_size:
                ec_sample_keys = random.sample(ec_keys, sample_size)
                ec_to_names = {k: ec_to_names[k] for k in ec_sample_keys}

        print(f"Loaded: {len(ec_to_names)} enzymes, {len(pathway_data)} pathways, {len(brenda_en)} BRENDA entries")

        # 1. Get EC to reaction mapping 
        print("Creating EC to reaction mapping...")
        ec_to_rxn = create_ec_to_reaction_mapping()
        
        # 2. extract protein information
        print("Preparing protein database...")
        protein_database = []
        # From BRENDA
        if brenda_en:
            for ec_number, data in brenda_en.items():
                # Use EC number to find enzyme names from ec_to_names
                if ec_number in ec_to_names:
                    protein_name = ec_to_names[ec_number][0] if ec_to_names[ec_number] else None
                    if protein_name:
                        protein_database.append({
                            'source': 'BRENDA',
                            'ec_number': ec_number,
                            'protein_name': protein_name,
                            'brenda_data': data
                        })
        # 3. From MetalPDB
        if metal_binding_data:
            for site_key, site_data in metal_binding_data.items():
                for chain in site_data.get('site_chains', []):
                    if 'molecule_name' in chain and chain['molecule_type'] == 'protein':
                        protein_database.append({
                            'source': 'MetalPDB',
                            'pdb_code': site_data.get('pdb_code'),
                            'protein_name': chain['molecule_name'],
                            'metal_binding': site_data.get('metal', {})
                        })

    except Exception as e:
        print(f"Error loading data sources: {e}")
    
        return []
    # 4. Track statistics for validation
    stats = {
        'total_enzymes': 0,
        'with_brenda_data': 0,
        'with_reactions': 0,
        'with_pathways': 0,
        'with_ko': 0,
        'with_metal_involvement': 0,
        'with_corrosion_mechanisms': 0
    }
    
    start_time_record = time.time()
    # 5. Prepare a list to store all records
    ec_records = []
    stats['total_enzymes'] = len(ec_to_names)

    print("Pre-computing enzyme class lookups...")
    ec_class_lookup = {}
    for ec_number in ec_to_names.keys():
        try:
            ec_prefix = '.'.join(ec_number.split('.')[:2])
            # Try exact match first
            if ec_prefix in enzyme_class:
                ec_class_lookup[ec_number] = enzyme_class[ec_prefix].strip()
            # Then try pattern match
            else:
                pattern_key = f"{ec_prefix}.-.-"
                if pattern_key in enzyme_class:
                    ec_class_lookup[ec_number] = enzyme_class[pattern_key]
        except Exception:
            pass  # Skip if error occurs
    
    # 6. Pre-compute KO lookups - do this once
    print("Pre-computing KO lookups...")
    ko_lookup = {}
    for ko, data in ko_ec.items():
        if isinstance(data, dict) and 'definition' in data:
            for ec_number in ec_to_names.keys():
                if f"[EC:{ec_number}]" in data['definition']:
                    if ec_number not in ko_lookup:
                        ko_lookup[ec_number] = []
                    ko_lookup[ec_number].append(ko)
    
    # 7. Populate from ec_to_names for basic enzyme names
    print("Building enzyme records...")
    for ec_number, names in ec_to_names.items():
        # Data validation for EC number format
        if not (ec_number.count('.') == 3 and all(part.isdigit() for part in ec_number.split('.'))):
            print(f"Warning: Invalid EC number format: {ec_number}")
            continue

        record ={
            'ec_number': ec_number,
            'enzyme_names': [str(name).strip() for name in (names if isinstance(names, list) else [names]) if name], # Ensure names are strings and stripped
            'protein_name': str(names[0]).strip() if isinstance(names, list) and names else str(names).strip() if names else None, # Strip protein name
            'enzyme_class': ec_class_lookup.get(ec_number),
            'pathways': [],
            'hierarchy': [],
            'ko': ko_lookup.get(ec_number, []),
            'reactions': [],
            'compounds': [],
            'modules': [],
            'metals_from_brenda': [],
            'corrosion_metals_from_brenda': []
        }
        # 8. Add pathways from EC-pathway mapping
        if ec_number in ec_pathway_mapping:
            for pathway_id in ec_pathway_mapping[ec_number]:
                # Standardize to map prefix
                std_id = pathway_id
                if pathway_id.startswith('ec'):
                    std_id = 'map' + pathway_id[2:]

                # Look up the pathway name
                if std_id in pathway_data:
                    pathway_name = pathway_data[std_id]
                    if pathway_name not in record['pathways']:
                        record['pathways'].append(pathway_name)

        # 9. Add pathways from KO data
        if ec_number in ko_ec and isinstance(ko_ec[ec_number], list):
            for path in ko_ec[ec_number]:
                if path not in record['pathways']:
                    record['pathways'].append(path)
        elif ec_number in ko_ec and isinstance(ko_ec[ec_number], dict) and 'pathway' in ko_ec[ec_number]:
            path = ko_ec[ec_number]['pathway']
            if path not in record['pathways']:
                record['pathways'].append(path)
        # Add KO IDs
        ko_ids = []
        for ko, data in ko_ec.items():
            if isinstance(data, dict) and 'definition' in data and f"[EC:{ec_number}]" in data['definition']:
                ko_ids.append(ko)
        record['ko'] = ko_ids

        if ko_ids:
            stats['with_ko'] += 1

        # 10. Build reaction list
        rxns = ec_to_rxn.get(ec_number, [])
        for rxn_id in rxns:
            if rxn_id in reaction_equation:
                eqn = reaction_equation.get(rxn_id, {}).get('equation', 'Unknown')
                record['reactions'].append({'id': rxn_id, 'equation': eqn})

                # Add compounds involved in this reaction
                for compound_id in reaction_equation.get(rxn_id, {}).get('compounds', []):
                    if compound_id in compound_info:
                        compound_data = compound_info[compound_id]
                        if compound_data not in record['compounds']:
                            record['compounds'].append(compound_data)

        if rxns:
            stats['with_reactions'] += 1

        # 11. Add module information
        for module_id, module_desc in module_info.items():
            if f"[EC:{ec_number}]" in module_desc:
                record['modules'].append({'id': module_id, 'description': module_desc})

        # 12. Reconcile metals from BRENDA with text mining
        record['metals_from_brenda'] = []
        record['corrosion_metals_from_brenda'] = []

        # 13. Add BRENDA metal information
        if brenda_en and ec_number in brenda_en:
            record['metals_from_brenda'] = brenda_en[ec_number].get('clean_metals', [])
            record['corrosion_metals_from_brenda'] = brenda_en[ec_number].get('corrosion_metals_from_brenda', [])
            stats['with_brenda_data'] += 1

        # 14. Check if EC number is valid
        ec_number = record['ec_number']
        has_valid_ec = ec_number.count('.') == 3 and all(part.isdigit() for part in ec_number.split('.'))
        # Initialize protein_name to None as a default
        protein_name = None
        # If no valid EC number or it failed any checks
        if not has_valid_ec and 'protein_name' in record and record['protein_name']:
            protein_name = record['protein_name'].lower()    

        ec_records.append(record)
        
    elapsed_time_record = time.time() - start_time_record
    print(f"Processing took {elapsed_time_record:.2f} seconds")

     # 15. Process records in batches to avoid memory issues
    print("Processing protein names and calculating scores...")
    batch_size = 1000  # Adjust based on available memory
    for i in range(0, len(ec_records), batch_size):
        batch = ec_records[i:i+batch_size]
        
        # 16. Process protein names in batch
        for rec in batch:
            # Process protein name matches
            try:
                if rec.get('protein_name') is None and rec.get('enzyme_names'):
                    rec['protein_name'] = rec['enzyme_names'][0]
                # Use indexing/hashing to speed up matching    
                protein_name = rec.get('protein_name', '')
                if not protein_name or len(protein_name) <= 2:
                    continue  
                
                protein_name_lower = protein_name.lower()
                for db_rec in protein_database:
                    db_protein_name = db_rec.get('protein_name', '')
                    if not db_protein_name:
                        continue
                        
                    db_protein_lower = db_protein_name.lower()
                    if (db_protein_lower == protein_name_lower or 
                        db_protein_lower in protein_name_lower or
                        protein_name_lower in db_protein_lower):
                            # Found a matching protein, merge data
                            if 'brenda_data' in db_rec:
                                # Process BRENDA data
                                brenda_data = db_rec['brenda_data']
                                if 'metals' in brenda_data and 'metals' not in rec:
                                    rec['metals'] = brenda_data['metals']
                                    
                            if 'metal_binding' in db_rec:
                                # Process MetalPDB data
                                if 'metal_binding_info' not in rec:
                                    rec['metal_binding_info'] = {}
                                    metal = db_rec['metal_binding'].get('symbol')
                                    if metal:
                                        rec['metal_binding_info'][metal] = db_rec['metal_binding']
            except Exception as e:
                print(f"Error processing protein name {rec.get('protein_name')}: {e}")

        # 17. Calculate scores for this batch
        for rec in batch:
            # Build text for scoring once
            enzyme_names = rec.get('enzyme_names', []) or []
            class_text = rec.get('enzyme_class', '') or ''
            pathways = rec.get('pathways', []) or []
            names_text = ' '.join(enzyme_names) + ' ' + class_text + ' ' + ' '.join(pathways)
            
            # Add reaction text for more context
            reaction_text = ""
            for r in rec.get('reactions', []):
                if isinstance(r, dict) and 'equation' in r:
                    reaction_text += " " + r['equation']

            # 18. Combined text for all scoring
            all_text = f"{names_text} {class_text} {' '.join(pathways)} {reaction_text}".lower()
            using_imported_modules = True
            # Score metals, mechanisms, etc. all at once
            try:
                if using_imported_modules:
                    # Use the scoring system module if available
                    score_results = cs.calculate_overall_scores(
                        all_text, 
                        brenda_metals=rec.get('metals_from_brenda', []),
                        pathways=' '.join(rec.get('pathways', []))
                    )
                    
                    # Merge the score results into the record
                    rec.update(score_results)
                    
                    # 19. Consolidate metals using the module function
                    rec['metals_consolidated'] = cs.consolidate_metal_terms(
                        rec.get('metals_from_brenda', []),
                        rec.get('metals_involved', [])
                    )
                else:
                    # 20. Fall back to original inline scoring logic# Score metals
                    metal_score, metal_matches = cs.score_keyword_matches(all_text, cs.metal_terms)
                    for metal in rec.get('metals_from_brenda', []):
                        if metal not in metal_matches:
                            metal_matches[metal] = 1.0
                    
                    rec['metals_involved'] = list(metal_matches.keys())
                    rec['metal_scores'] = metal_matches
                    rec['overall_metal_score'] = float(metal_score)
                    
                    # 21.Score corrosion mechanisms  
                    corrosion_score, corrosion_matches = cs.score_keyword_matches(all_text, cs.corrosion_mechanisms)
                    rec['corrosion_mechanisms'] = list(corrosion_matches.keys())
                    rec['corrosion_mechanism_scores'] = corrosion_matches
                    rec['overall_corrosion_score'] = float(corrosion_score)
                    
                    # 22. Score synergies
                    synergy_score, synergy_matches = cs.score_keyword_matches(all_text, cs.corrosion_synergies)
                    rec['corrosion_synergies'] = list(synergy_matches.keys())
                    rec['corrosion_synergy_scores'] = synergy_matches
                    rec['overall_synergy_score'] = float(synergy_score)
                    
                    # 23. Score functional categories
                    functional_terms = {cat: details['terms'] for cat, details in cs.functional_categories.items()}
                    func_score, func_matches = cs.score_keyword_matches(all_text, functional_terms)
                    weighted_func_matches = {}
                    for cat, match_score in func_matches.items():
                        original_weight = functional_categories[cat]['score']
                        weighted_func_matches[cat] = match_score * original_weight

                    rec['functional_categories'] = [{"category": cat, "score": score} 
                                                 for cat, score in weighted_func_matches.items()]
                    rec['overall_functional_score'] = float(sum(weighted_func_matches.values()))
                    
                    # 24. Score organic processes
                    organic_score, organic_matches = cs.score_keyword_matches(all_text, cs.organic_categories)
                    rec['organic_processes'] = list(organic_matches.keys())
                    rec['organic_process_scores'] = organic_matches
                    rec['overall_organic_process_score'] = float(organic_score)
                    
                    # 25. Score keyword groups
                    keyword_score, keyword_matches = cs.score_keyword_matches(all_text, cs.corrosion_keyword_groups)
                    rec['corrosion_keyword_groups'] = list(keyword_matches.keys())
                    rec['corrosion_keyword_scores'] = keyword_matches
                    rec['overall_keyword_score'] = float(keyword_score)
                    
                    # 26. Consolidate metals
                    rec['metals_consolidated'] = cs.consolidate_metal_terms(
                        rec.get('metals_from_brenda', []),
                        rec.get('metals_involved', [])
                    )
                
                # 27. Update statistics (same for both paths)
                if rec['metals_consolidated']:
                    stats['with_metal_involvement'] += 1
            
                if rec['corrosion_mechanisms']:
                    stats['with_corrosion_mechanisms'] += 1
                    
            except Exception as e:
                print(f"Error scoring data for {rec.get('ec_number')}: {e}")
                
    # 28. Process metal binding in a separate pass
    print("Processing metal binding information...")
    try:
        # Add metal binding information to records
        for rec in ec_records:
            rec['metal_binding_info'] = {}

            # Check metals from BRENDA or detected in text
            all_metals = set(rec.get('metals_consolidated', []))

            for metal in all_metals:
                # Try to map to standard symbol
                for metal_name, symbol in cs.metal_mapping.items():
                    if metal_name in metal.lower() or symbol.lower() in metal.lower():
                        # Check if we have binding data for this metal
                        if symbol in metal_patterns.get('residue_binding', {}):
                            # Get top binding residues
                            residue_counts = metal_patterns['residue_binding'][symbol]
                            top_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)[:5]

                            rec['metal_binding_info'][symbol] = {'common_residues': [res for res, count in top_residues],  'binding_count': sum(residue_counts.values())}
                                                                  
                # Check for common coordination geometries
                geometries = [key.split('_')[2:] for key, count in metal_patterns['coordination'].items() if key.startswith(f"{symbol}_")]
                if geometries:
                    rec['metal_binding_info'][symbol]['common_geometries'] = geometries[:3]

    except Exception as e:
        print(f"Error processing metal binding data: {e}")

    # 29- Process KO hierarchy information
    print("Processing KO hierarchy...")
    if 'D' in ko_hierarchy:
        for ko, info in ko_hierarchy.get('D', {}).items():
            for ec in info.get('ec_numbers', []):
                # find matching records
                for rec in ec_records:
                    if rec['ec_number'] == ec:
                        try:
                            parent_c = info.get('parent')
                            if parent_c and 'C' in ko_hierarchy and parent_c in ko_hierarchy['C']:
                                path_info = ko_hierarchy['C'][parent_c]
                                parent_b = path_info.get('parent')
                                if parent_b and 'B' in ko_hierarchy and parent_b in ko_hierarchy['B']:
                                    hi_category = ko_hierarchy['B'][parent_b].get('name', '')
                                    pathway = path_info.get('name', '')

                                    # Use sets to efficiently track unique values
                                    if pathway and pathway not in rec['pathways']:
                                        rec['pathways'].append(pathway)

                                    hierarchy = f"{hi_category} > {pathway}"
                                    if hierarchy and hierarchy not in rec['hierarchy']:
                                        rec['hierarchy'].append(hierarchy)
                        except Exception as e:
                            print(f"Error processing KO hierarchy for {rec.get('ec_number')}: {e}")
    
    # 30. Calculate pathway statistics
    print("Calculating pathway statistics...")
    for rec in ec_records:
        if rec['pathways']:
            stats['with_pathways'] += 1
    
    # 31. Calculate final corrosion relevance scores using the scoring module
    print("Calculating corrosion relevance scores...")
    for rec in ec_records:
        try:
            enzyme_names = rec.get('enzyme_names', []) or []
            enzyme_class_text = rec.get('enzyme_class', '') or ''
            pathways = rec.get('pathways', []) or []
            
            if using_imported_modules:
                # 32. Calculate pathway scores using the module 
                pathway_score, pathway_category_scores = cs.calculate_pathway_score(
                    rec.get('pathways', []), enzyme_names, enzyme_class_text
                )
                # Extract mechanisms from pathways
                if rec.get('pathways'):
                    pathway_text = '\n'.join(rec.get('pathways', []))
                    pathway_mechanisms = cs.assign_mechanism_from_pathway(pathway_text)
                    existing_mechanisms = rec.get('corrosion_mechanisms', [])
                    rec['corrosion_mechanisms'] = list(set(existing_mechanisms + pathway_mechanisms))

                # 33. Check for metal-organic synergy
                synergy_score = cs.check_metal_organic_synergy(
                    rec.get('metals_consolidated', []), 
                    enzyme_names, pathways
                )
                rec['metal_organic_synergy_score'] = synergy_score
                
                # 34. Calculate final corrosion relevance score using the module function
                corrosion_relevance_score, corrosion_relevance = cs.calculate_corrosion_relevance_score(
                    rec.get('overall_metal_score', 0),
                    rec.get('overall_corrosion_score', 0),
                    pathway_score,
                    rec.get('overall_organic_process_score', 0),
                    rec.get('overall_keyword_score', 0),
                    synergy_score,
                    rec.get('overall_functional_score', 0)
                )
            
                # 45.Add pathway scores to record (same for both paths)
                rec['pathway_category_scores'] = pathway_category_scores
                rec['overall_pathway_category_score'] = sum(pathway_category_scores.values())
                rec['pathway_score'] = float(pathway_score)
                
                # 46. Store final scores in record (same for both paths)
                rec['corrosion_relevance_score'] = corrosion_relevance_score
                rec['corrosion_relevance'] = corrosion_relevance
            
        except Exception as e:
            print(f"Error calculating corrosion score for {rec.get('ec_number')}: {e}")
            rec['corrosion_relevance_score'] = 0
            rec['corrosion_relevance'] = 'unknown'        

    # 47. Filter records without content
    print("Filtering records...")
    filtered_ec_records = []
    for rec in ec_records:
        protein_name = (rec.get('protein_name') or "").lower()
        enzyme_names = rec.get('enzyme_names', [])
        enzyme_names = [] if enzyme_names is None else enzyme_names
        ec_number = rec.get('ec_number', "")

        # Condition 1: At least one valid identifier must be present
        has_valid_protein = "uncharacterized" not in protein_name and len(protein_name) > 2
        has_valid_enzyme = any(len(name) > 2 for name in enzyme_names)
        has_valid_ec = ec_number.count('.') == 3 and all(part.isdigit() for part in ec_number.split('.'))
      
        # Condition 2: Check for valuable data that should be preserved
        has_mechanisms = len(rec.get('corrosion_mechanisms', [])) > 0
        has_pathways = len(rec.get('pathways', [])) > 0
        has_metals_consolidated = len(rec.get('metals_consolidated', [])) > 0

        # Include record if it meets either condition
        if (has_valid_protein or has_valid_enzyme or has_valid_ec) or (has_mechanisms or has_pathways or has_metals_consolidated):
            filtered_ec_records.append(rec)

    # 48. Replace original list with filtered version
    ec_records = filtered_ec_records

    # Print summary statistics
    print("\nMetabolism Database Summary:")
    print(f"Total enzyme records: {stats['total_enzymes']}")
    print(f"Records with BRENDA data: {stats['with_brenda_data']} ({stats['with_brenda_data']/stats['total_enzymes']*100:.1f}%)")
    print(f"Records with reactions: {stats['with_reactions']} ({stats['with_reactions']/stats['total_enzymes']*100:.1f}%)")
    print(f"Records with pathways: {stats['with_pathways']} ({stats['with_pathways']/stats['total_enzymes']*100:.1f}%)")
    print(f"Records with KO terms: {stats['with_ko']} ({stats['with_ko']/stats['total_enzymes']*100:.1f}%)")
    print(f"Records with corrosion mechanisms: {stats['with_corrosion_mechanisms']} ({stats['with_corrosion_mechanisms']/stats['total_enzymes']*100:.1f}%)")

    # 49. Validate the data - check for missing essential fields
    validation_issues = []
    for i, rec in enumerate(ec_records):
        if not rec.get('ec_number'):
            validation_issues.append(f"Record {i} missing EC number")
        if not rec.get('enzyme_names'):
            validation_issues.append(f"EC {rec.get('ec_number')} missing enzyme names")

    if validation_issues:
        print("\nValidation Issues:")
        for issue in validation_issues[:10]:  # Show first 10 issues
            print(f"- {issue}")
        if len(validation_issues) > 10:
            print(f"...and {len(validation_issues) - 10} more issues")
    else:
        print("\nValidation: All records have essential fields.")

    return ec_records
#=========================================  RUNNING THE PIPELINE ================================================ 

from pathlib import Path
import json

if __name__ == "__main__":
    ec_records = create_metabolism_database(sample= 1500) 
    #json_path = Path("/kaggle/working/")
    output_large = Path("/home/beatriz/MIC/output_large") 
    json_path = output_large / "ec_records.json"
    
    with open(json_path, 'w', encoding="utf-8") as f:
        json.dump(ec_records, f, ensure_ascii=False, indent=2)
