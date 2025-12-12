#!/usr/bin/env python
# coding: utf-8

# # 7 Database Functional Notebook
# This notebook was originally part of notebook 6, because the thematic was part of it,  but the notebook was too heavy and hence was divided.
# This part comprises the creation of a compendium database with later feeds the picrust2 result data. The last pipeline is the classification and filtering of the protein-bacteria to narrow down the most important pairs on corrosion systems.

# _________________________
# ## importing scoring system
# ___________________________

# In[1]:


import os
import sys
from pathlib import Path
sys.path.append(os.path.abspath('..'))  # Ensures the project root is in Python's search path

if Path("/kaggle").exists():
    print ("Running Kaggle environment")
    # Create directory structure
    get_ipython().system('mkdir -p corrosion_scoring_v3')
    get_ipython().system('wget -O corrosion_scoring_v3/__init__.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/__init__.py')
    get_ipython().system('wget -O corrosion_scoring_v3/global_terms.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/global_terms.py')
    get_ipython().system('wget -O corrosion_scoring_v3/term_processor.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/term_processor.py')
    get_ipython().system('wget -O corrosion_scoring_v3/config.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/config.py')
    get_ipython().system('wget -O corrosion_scoring_v3/score_calculator.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/score_calculator.py')
    get_ipython().system('wget -O corrosion_scoring_v3/synergy_detector.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/synergy_detector.py')
    get_ipython().system('wget -O corrosion_scoring_v3/name_utils.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/name_utils.py')
    get_ipython().system('wget -O corrosion_scoring_v3/utils_ec.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/utils_ec.py')
    get_ipython().system('wget -O corrosion_scoring_v3/validators.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/validators.py')
    get_ipython().system('wget -O corrosion_scoring_v3/exceptions.py https://raw.githubusercontent.com/MagicAlex238/2_Micro/main/corrosion_scoring_root/corrosion_scoring_v3/exceptions.py')
    # Add current directory to path
    sys.path.append(os.getcwd())
    print("Running in Kaggle environment")  
else:
    print("Running in local (VSCode) environment")# Silencing the imports after stable package
    os.system("pip uninstall -y corrosion_scoring_v3 || true")
    #os.system("pip cache purge")
    #os.system("pip install --force-reinstall git+https://github.com/MagicAlex238/2_Micro.git@refactor-scoring-system#subdirectory=corrosion_scoring_root/corrosion_scoring_v3")
    #os.system("pip install git+https://github.com/MagicAlex238/2_Micro.git@refactor-scoring-system#subdirectory=corrosion_scoring_root/corrosion_scoring_v3")
    # ensuring the path is set correctly:
    sys.path.insert(0, "/home/beatriz/MIC/2_Micro/corrosion_scoring_root") 

import corrosion_scoring_v3 as cs
from corrosion_scoring_v3.term_processor import TermProcessor
from corrosion_scoring_v3.name_utils import enhanced_clean_protein_name, clean_protein_name
from corrosion_scoring_v3.utils_ec import normalize_ec_id, strip_all_ec_tokens , normalize_listlike, standardize_metals_list, standardize_metal_symbol
from corrosion_scoring_v3.validators import ValidationError
from corrosion_scoring_v3.exceptions import ScoringError, TextMiningError, SynergyDetectionError
from corrosion_scoring_v3.global_terms import (metal_terms_dict, functional_categories_dict, corrosion_synergies_dict,
    mechanisms_dict, pathway_dict, operational_environmental_factors_dict, metal_mapping)

# create processors (use cs.<name> for consistency)
fc_processor = TermProcessor(cs.functional_categories_dict)
metal_processor = TermProcessor(cs.metal_terms_dict)
synergy_processor = TermProcessor(cs.corrosion_synergies_dict)
mechanisms_processor = TermProcessor(cs.mechanisms_dict)
pathway_processor = TermProcessor(cs.pathway_dict)
ope_processor = TermProcessor(cs.operational_environmental_factors_dict)

processors =   {'fc_processor': fc_processor, 'metal_processor': metal_processor,
    'synergy_processor': synergy_processor}

# Initialize v3 components with existing processors
# ---- Initialize v3 system components ----
config = cs.ScoringConfig()
text_miner = cs.TextMiner(config)
text_miner.processors = processors 

score_calculator = cs.ScoreCalculator(config)
synergy_detector = cs.SynergyDetector(config)
metal_mapping = cs.metal_mapping 


# __________________
# ## Libraries Imports
# _______________

# In[2]:


# Standard library imports
import sys
import os
from pathlib import Path
import ast
import subprocess
import logging
import time
from datetime import datetime
import shutil
from io import StringIO
import psutil
import re
from IPython import get_ipython
from IPython.display import display

# Data processing and analysis
import pandas as pd
from pandas.api.types import is_categorical_dtype
from pandas.api.types import CategoricalDtype
import numpy as np
import openpyxl
import seaborn as sns
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib.colors import to_rgba, LinearSegmentedColormap
import matplotlib.patches as mpatches
from collections import defaultdict
from collections import Counter 
from itertools import chain
from lxml import etree

# Machine learning and statistical analysis
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
from sklearn.decomposition import PCA, NMF
from sklearn.cluster import AgglomerativeClustering, KMeans
from scipy.cluster import hierarchy
from sklearn.manifold import TSNE
import umap
import scipy
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
import scipy.cluster.hierarchy as sch
from statsmodels.stats.multitest import multipletests
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr, kruskal, mannwhitneyu
from kneed import KneeLocator
from scipy.signal import savgol_filter
from community import community_louvain
from joblib import Parallel, delayed

# Bioinformatics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biom import Table, load_table
from biom.util import biom_open
import xml.etree.ElementTree as ET

# Utility libraries
import gzip
import random
from natsort import natsorted
from typing import Tuple, Set, Iterable, Union
from typing import Any, Optional, Dict, List
import pickle
import gc
import joblib
import h5py
import os
import csv
import json
import pyarrow.parquet as pq
from collections import defaultdict
import warnings

import logging
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

os.environ['DISPLAY'] = ':0'


# ____________________
# ## Paths
# _______________

# In[3]:


# Initial environment detection and package installation
is_colab = "google.colab" in sys.modules
is_kaggle = Path("/kaggle").exists()
is_vscode = not (is_colab or is_kaggle)

if is_colab:
    print("Running in Google Colab environment")
    get_ipython().run_line_magic('pip', 'install psutil, biopython, biom-format, umap-learn, fuzzywuzzy, lxml pandas, pyarrow, openpyxl, scipy, python-Levenshtein, -U kaleido, statsmodels, kneed, natsort, adjustText')
    from google.colab import drive
    drive.mount('/content/drive')
    base_dir  = Path("/content/drive/My Drive")
    os.chdir(base_dir) # Change directory to base_dir
    # Directory to output large files # eccontris, compilated db
    large_dir = base_dir / "MIC"
    data_ref = large_dir / "2_Micro/data_ref" #input dir 
    data_qiime = large_dir / "2_Micro/data_qiime" #input dir 
    data_picrust = large_dir/ "2_Micro/data_picrust" #main output dir
    data_picrust.mkdir(parents=True, exist_ok=True) 

elif is_kaggle:
    print("Running in Kaggle environment")
    get_ipython().system('pip install psutil, biopython, biom-format, umap-learn, fuzzywuzzy, lxml pandas, pyarrow, openpyxl, scipy, python-Levenshtein, -U kaleido, statsmodels, kneed, natsort, adjustText')
    base_dir  = Path("/kaggle/input") 
    data_Ref  = base_dir / "place dataref" # revisit in kaggle "2_Micro/data_ref" #input dir
    data_qiime  = base_dir / "crust" #  needs revisiting
    # Directory to output large files # eccontris, compilated db
    large_dir =  Path("/kaggle/working/")  
    data_picrust = large_dir / "data_picrust" # directory assignment no creation

else:
    print("Running in VSCode/local environment")
    base_dir = Path("/home/beatriz")
    large_dir  = base_dir /"MIC"  
    data_qiime = large_dir / "2_Micro" / "data_qiime"  # this already come from former notebook and it is not to make
    data_Ref  = large_dir /  "2_Micro" / "data_Ref" # this already come from former notebook and it is not to make
    data_picrust  = large_dir / "2_Micro" / "data_picrust" 
    data_picrust.mkdir(parents=True, exist_ok=True)


output_dir = data_picrust
#datasets large galaxies and databases
input_galaxy = large_dir  / "data_galaxies"
db_dir = large_dir / "Databases"
newick_path = input_galaxy / "Galaxy7-PICRUSt2-Full-pipeline-on-data-2-and-data-1-Tree-reference-study-16S-sequences.newick" 
# Directory to output large files # eccontris, compilated db
output_large = large_dir / "output_large"
output_large.mkdir(parents=True, exist_ok=True)

abundance_excel= data_Ref / "merged_to_sequence.xlsx" #C:\home\beatriz\MIC\2_Micro\data_Ref\merged_to_sequence.xlsx
fasta_file_final = data_qiime / "results_match_gg/final_sequences_gg.fasta"
aligned_fasta = data_qiime / "results_match_gg/aligned-dna-sequences_gg.fasta"
# Create output directory if it doesn't exist
output_base = output_dir / "output_base" # directory4
output_base.mkdir(parents=True, exist_ok=True)
# Print the paths for verification
print(f"Using base_dir: {base_dir}")
print(f"Using abundance_excel: {abundance_excel}")
print(f"Using fasta_file_final: {fasta_file_final}")
print(f"Using output_base: {output_base}")
print(f"Using large_dir: {large_dir}")
print(f"Using db_dir: {db_dir}")
print(f"Using input_galaxy: {input_galaxy}")
print(f"Using output_large: {output_large}")


# In[4]:


# Integrated taxa from origin genus as headers with levels 6 for the genera, 7 for the GID, muss be cleaned
Integrated_T = pd.read_excel(abundance_excel, sheet_name='core_check_usual_taxa', header=[0,1,2,3,4,5,6,7], engine ='openpyxl')
# Drop first row (index 0) and first column in one chain
Integrated_T = Integrated_T.drop(index=0).drop(Integrated_T.columns[0], axis=1)
Integrated_T= Integrated_T.astype({'Sites': str})
Integrated_T['Sites'] = Integrated_T['Sites'].fillna('Source')
# Remove 'Unnamed' level names
Integrated_T.columns = Integrated_T.columns.map(lambda x: tuple('' if 'Unnamed' in str(level) else level for level in x))
# Changing dtypes to category whiles respecting structure
Integrated_T["Category"] = Integrated_T["Category"].astype("Int64")
Integrated_T= Integrated_T.set_index("Sites")
pre_Integrated = Integrated_T.T


# In[5]:


# Define category dict outside
category_dict = Integrated_T.T.iloc[0, 0:-1].to_dict()


# In[6]:


#temporal function
def notify_complete():

    # Visual notification
    print(f"\n🔔 TASK FINISHED at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} 🔔")

    # Try audio methods (no installation required)
    try:
        # Method 1: PowerShell (Windows/WSL)
        os.system('powershell.exe -c "[console]::beep(4800,8800)"')
    except:
        try:
            # Method 2: Terminal bell
            print('\a\a\a')
            print('\a\a\a')
            print('\a\a\a')
            sys.stdout.flush()
        except:
            # Method 3: Just visual
            print("🚨🚨🚨 DONE 🚨🚨🚨")


# ____________
# # 1. Building a Dictionary from Databases
# ____________
# The data from picrust2 was insuficient to properly identify the protein-genus pairs and it was necesary to extend the search by including autoritative databases to feed into the data. Some steps were done for best practice:
# Data Normalization and Mapping: It was ensured that all protein/EC data (including EC numbers, KO numbers, and reaction IDs) were parsed and cleaned. This identifiers were mapped to their corresponding metabolic pathways using databases such as KEGG, MetaCyc, and BioCyc.
# 
# Identifying Metal-Related Proteins: cross-reference proteins with metal-related databases (BRENDA, MetalPDB, TransportDB) where crossreferenced to flag those with direct metal-binding or metal-transport roles. Then consolidate similar metal terms (e.g., “iron”, “Fe”, “ferric”) into a unified field to improve consistency in later analyses.
# 
# Final Assembly: Data was compiled into a final dictionary/table that includes all relevant columns (Protein, EC/KO, Metabolism, Pathway, Metal Interaction, MIC Function). This allows to search programmatically for the functional roles of proteins that are influential in corrosion studies.  

# ____________
# ### Databases Paths
# ___________

# In[7]:


def setup_paths() -> Dict[str, Path]:
    """Set up paths for database access"""

    # Database paths
    db_paths = {
        'enzyme': db_dir / 'enzyme',
        'enzyme_class': db_dir / 'enzclass.txt',
        'enzyme_brenda' : db_dir/ 'brenda_2024.txt',
        'ko': db_dir / 'kegg_ko.txt',
        'ko_hierarchy': db_dir / 'ko_hierarchy.txt',
        'ko_brite': db_dir / 'ko00001.keg',
        'pathway': db_dir / 'pathway',
        'module': db_dir / 'module',
        'reaction': db_dir / 'reaction',
        'compound': db_dir / 'compound',
        'metalpdb': db_dir / 'flat_db_file.xml', #MetalPDB
        'ko_pathway': db_dir / 'ec_pathway.list' #Mapping between EC numbers and KEGG pathway IDs
    }

    return db_paths

#  Calling the paths
if __name__ == "__main__":
    paths = setup_paths()
    # Print paths to verify
    for db_name, path in paths.items():
        print(f"{db_name}: {path}")
        print(f"Exists: {path.exists()}")


# '''https://www.brenda-enzymes.org/download.php
# 
# Chang A., Jeske L., Ulbrich S., Hofmann J., Koblitz J., Schomburg I., Neumann-Schaal M., Jahn D., Schomburg D.
# BRENDA, the ELIXIR core data resource in 2021: new developments and updates. (2021), Nucleic Acids Res., 49:D498-D508.
# DOI: 10.1093/nar/gkaa1025 PubMed: 33211880
# Brenda Enzyme Database. (n.d.). In BRENDA. Retrieved from https://www.brenda-enzymes.org (APA citation format)
# 
# wget https://www.enzyme-database.org/downloads/enzyme-database.sql.gz
# wget -O ec_pathway.list "https://rest.kegg.jp/link/pathway/ec",
# Pathways rsync -avz rsync://rest.kegg.jp/kegg/pathway/ 
# Reactions !wget -c "ftp://ftp.genome.jp/pub/kegg/reaction/reaction.tar.gz" # This I can not use because the id are no compatible
# Chemical compounds database.wget https://biocyc.org/download.shtml. wget https://www.brenda-enzymes.org/download.php
# MetalPDB in 2018: a database ofbL., Andreini C.
# Nucleic Acids Res. 2018 Jan;46(D1):D459-D464. [PMID: 29077942]
# MetalPDB: a database of metal sites in biological macromolecular structures. (n.d.). Retrieved from http://metalpdb.cerm.unifi.it (APA citation format)
# '''

# ## 1.1. Importing and Parsing Authoritative Databases
# __Brenda enrichent sources__
# | Final Column             | BRENDA Tag/Field                  | Used for Scoring? | How to Populate/Notes                                                             | Normalization/Function                  |
# |--------------------------|-----------------------------------|-------------------|-----------------------------------------------------------------------------------|-----------------------------------------|
# | ec                       | ID (normalized EC)                | Linkage only      | Parse EC from ID; do not overwrite canonical EC from Uniprot                      | normalize_ec_id                         |
# | protein_name             | DE                                | No (ref only)     | Parse DE; keep as BRENDA reference; do not overwrite Uniprot protein_name    | enhanced_clean_protein_name (optional)  |
# | enzyme_class             | Derived from EC (first digit)     | No (ref only)     | Map first EC digit to class name (e.g., 1 oxidoreductases)                        | ec_class_mapping                        |
# | functional_categories    | DE, enzyme_class, SP/NSP, IN      | Yes               | Run fc_processor on these texts only; exclude GI/RE/PM from FC scoring            | TermProcessor(functional_categories)    |
# | metals                   | ME (+ detected categories)        | Yes               | Consolidate raw ME terms with detected metal categories                           | consolidate_metal_terms                 |
# | corrosion_synergies      | From FC results (preferred)       | Yes               | Use FC co-occurrence; fall back to synergy keywords when needed                   | scoring_system synergy logic            |
# | corrosion_mechanisms     | GI and RE lines containing “mechanism” | No            | Extract mechanism snippets from GI/RE; also de-duplicate                          | _extract_mechanism_snippet + cleaning   |
# | organism                 | OS                                | No                | Parse OS; use primarily for genus-level matching during enrichment                | as-is (optionally standardize genus)    |
# | organism_classification  | OC                                | No                | Parse OC for context                                                              | as-is                                   |
# | gene_name                | GN                                | No                | Parse GN; treat as ambiguous unless paired with OS (genus)                        | as-is                                   |
# | reference_number         | RN                                | No                | Trace/debug only                                                                  | as-is                                   |
# | reference_position       | RP                                | No                | Trace/debug only                                                                  | as-is                                   |
# | physiological_functions  | GI entries mentioning “physiological function” | No      | Keep as context; can help human review                                            | as-is                                   |
# | protein_modifications    | PM                                | No                | Keep for context (e.g., PTMs)                                                     | as-is                                   |
# | reactions                | RE                                | No (not import)|                 | —                                       |
# 
# Notes
# - Keep Uniprot (EC, protein_name, genus) canonical; do not overwrite them with BRENDA values.
# - For enrichment, only harvest BRENDA fields when OS genus matches the Uniprot genus. Everything else stays as auxiliary context.
# - FC inputs: DE, enzyme_class, SP/NSP, IN. Exclude GI/RE/PM from FC scoring to respect the design.

# __________________________________________________
# ### Parse Brenda
# _________________________________
# 
# https://www.brenda-enzymes.org/download.php
# 
# Chang A., Jeske L., Ulbrich S., Hofmann J., Koblitz J., Schomburg I., Neumann-Schaal M., Jahn D., Schomburg D.
# BRENDA, the ELIXIR core data resource in 2021: new developments and updates. (2021), Nucleic Acids Res., 49:D498-D508.
# DOI: 10.1093/nar/gkaa1025 PubMed: 33211880

# In[8]:


logger = logging.getLogger(__name__)

# =========================
# Tag parsing & utilities
# =========================

# Accept either:
#   TAG value
# or
#   LONG_SECTION SHORT_TAG value
TAG_LINE_RE = re.compile(r'^([A-Z_]{2,30})(?:\s+([A-Z]{2,10}))?\s+(.*)$')
INDEX_PREFIX_RE = re.compile(r'^#(\d+(?:,\d+)*)#\s*')      # "#1,4,7# text ..."
CITATION_HASH_RE = re.compile(r'#\s*\d+(?:\s*,\s*\d+)*#?')          # "#1#" or "#1,2,3#" "#1,2,3#" and "#1,2,3".
CITATION_ANGLE_RE = re.compile(r'<[^>]*>')                 # "<44>"
WS_RE = re.compile(r'\s+')
MECH_RE = re.compile(r'\bmechanism\b', re.IGNORECASE)
_COFACTOR_LEXICON = ('NADPH', 'NADP+', 'NADP', 'NADH', 'NAD+', 'NAD','FAD', 'FMN', 'HEME', 'BIOTIN', 'SAM', 'TPP', 'PLP', 'THF', 'COA')

# Long section names to short tags
LONG_TO_SHORT = {
    'RECOMMENDED_NAME': 'RN',
    'SYSTEMATIC_NAME': 'SN',
    'SYNONYMS': 'SY',
    'COFACTOR': 'CF',
    'ACTIVATING_COMPOUND': 'AC',
    'INHIBITORS': 'IN',
    'METALS_IONS': 'ME',
    'SUBSTRATE_PRODUCT': 'SP',
    'NATURAL_SUBSTRATE_PRODUCT': 'NSP',
    'GENERAL_INFORMATION': 'GI',
    'GENERAL_STABILITY': 'GS',
    'STORAGE_STABILITY': 'SS',
    'OXIDATION_STABILITY': 'OS',
    'PH_OPTIMUM': 'PO',
    'PH_STABILITY': 'PS',
    'PH_RANGE': 'PHR',
    'TEMPERATURE_OPTIMUM': 'TO',
    'TEMPERATURE_STABILITY': 'TS',
    'TEMPERATURE_RANGE': 'TR',
    'ORGANIC_SOLVENT_STABILITY': 'OSS',
    'IONIC_STRENGTH': 'IS',
    'PROTEIN_MODIFICATIONS': 'PM',
    'POSTTRANSLATIONAL_MODIFICATION': 'PM',
    'PROTEIN': 'PR',
    'CLONED': 'CL',
    'EXPRESSION': 'EXP',  # ignored
    # Intentionally NOT mapping REFERENCE/REACTION/MW/PI/AP/EN/CR/LO/ST/KM/SA/PU here (ignored)
}

# Short tags we accept from the file
VALID_SHORT_TAGS = {
    'ID', 'PR', 'RN', 'SN',
    'SY',
    'CF', 'AC',
    'ME', 'SP', 'NSP', 'IN',
    'GI', 'GS', 'SS', 'OS', 'PO', 'PS', 'PHO', 'PHS', 'PHR',
    'TO', 'TS', 'TR', 'OSS', 'IS',
    'PM',
    'CL',
    'EXP',
    # tolerated but ignored:
    'RE', 'RF', 'MW', 'PI', 'AP', 'EN', 'CR', 'LO', 'ST', 'KM', 'SA', 'PU', 'RT', 'BR'
}

# Safe aliases (do NOT alias REN→RN to avoid name pollution)
SHORT_TAG_ALIAS = {
    'PHO': 'PO',
    'PHS': 'PS',
}

# Kept environmental tags (captured for later; not used in current scoring)
ENV_TAGS = {'OS', 'GS', 'SS', 'PO', 'PS', 'PHR', 'TO', 'TS', 'TR', 'OSS', 'IS'}

EC_CLASS_MAP = {
    '1': 'oxidoreductases',
    '2': 'transferases',
    '3': 'hydrolases',
    '4': 'lyases',
    '5': 'isomerases',
    '6': 'ligases',
    '7': 'translocases',
}

def enzyme_class_from_ec(ec: str) -> str:
    if not ec:
        return ""
    return EC_CLASS_MAP.get(ec.split('.', 1)[0], "")

def brenda_sanitize_text(text: str) -> str:
    if not isinstance(text, str):
        if isinstance(text, (list, tuple)):
            text = ' '.join(map(str, text))
        else:
            text = '' if text is None else str(text)

    t = CITATION_HASH_RE.sub('', text)       # "#1,2#"
    t = CITATION_ANGLE_RE.sub('', t)         # "<44>"
    t = re.sub(r'\{[^}]*\}', '', t)          # "{r}" blocks
   # Check for chemistry tags and handle parentheses accordingly
    chemistry_tags = ['CF', 'ME', 'SP', 'NSP', 'RE']
    has_chemistry_tags = any(tag in t for tag in chemistry_tags)

    if has_chemistry_tags:
        # For chemistry tags: replace parentheses with spaces
        t = re.sub(r'[()]', ' ', t)
    else:
        # For non-chemistry tags: remove parentheses and their contents entirely
        t = re.sub(r'\([^)]*\)', '', t)
    t = WS_RE.sub(' ', t).strip()
    t = re.sub(r'^\s*[,;]?\s*\d[\d,\s]*\s*', '', t)  # leading enumerations

    # strip trailing section ghost labels
    t = re.sub(
        r'\s+(?:METALS_IONS|PH_RANGE|INHIBITORS|SUBSTRATE_PRODUCT|NATURAL_SUBSTRATE_PRODUCT|'
        r'PH_OPTIMUM|PH_STABILITY|TEMPERATURE_OPTIMUM|TEMPERATURE_STABILITY|SPECIFIC_ACTIVITY|'
        r'SUBUNITS|MOLECULAR_WEIGHT)\s*$',
        '', t
    )
    t = re.sub(r'\s*[:;,\-–]\s*$', '', t)
    assert isinstance(t, str), f"sanitize produced non-str: {type(t)}"
    return t

def is_valid_text(text):
    if not text or len(text) < 10:
        return False
    if '...' in text or '?' in text:
        return False
    if text.lower() in ('unknown', 'n/a', 'none'):
        return False
    return True

def is_valid_equation(eq: str) -> bool:
    # Basic filter for chemical equations
    if not eq or len(eq) < 10:
        return False
    if eq.startswith('+') or '...' in eq or eq.count('=') > 1:
        return False
    if not any(token in eq for token in ['+', '=', '→', '⇌']):
        return False
    return True

def extract_indices_and_value(value: str) -> Tuple[Optional[List[str]], str]:
    """Extract leading #idx[,idx]#; return (indices, rest)."""
    if not value:
        return None, ""
    m = INDEX_PREFIX_RE.match(value)
    if not m:
        return None, value.strip()
    idxs = m.group(1).split(',')
    rest = value[m.end():].strip()
    return idxs, rest

def extract_genus(species: str) -> str:
    """Return genus from species name; handle 'Candidatus X ...'."""
    if not species:
        return ""
    tokens = species.strip().split()
    if not tokens:
        return ""
    if tokens[0].lower() == 'candidatus' and len(tokens) > 1:
        return tokens[1]
    return tokens[0]

def contains_mechanism(text: str) -> bool:
    return bool(text and MECH_RE.search(text))

def mechanism_snippet(text: str) -> str:
    if not text:
        return ""
    m = MECH_RE.search(text)
    return text[m.start():].strip() if m else text.strip()

def _is_gene_token(tok: str) -> bool:
    if not tok:
        return False
    t = tok.strip().strip('.,;:()[]{}').replace('"','').replace("'",'')
    # Reject obvious non-genes (add/remove terms as needed)
    BAD = {'gene', 'locus', 'tag', 'orf', 'protein', 'domain', 'subunit', 'fragment', 'hypothetical', 'functional', 'cluster'}
    if t.lower() in BAD:
        return False
    # Typical gene tokens: short alnum with optional mixed case, underscores or dashes
    if not re.match(r'^[A-Za-z][A-Za-z0-9._-]{1,10}$', t):
        return False
    return True

# Conservative gene capture patterns from CL (cloned) lines
CL_GENE_PATTERNS = [
    re.compile(r'\bgene\s+([A-Za-z0-9_\-\.]{2,20})', re.IGNORECASE),
    re.compile(r'\b([A-Za-z0-9_\-\.]{2,20})\s+gene\b', re.IGNORECASE),
    re.compile(r'\blocus\s+tag\s+([A-Za-z0-9_\-\.]{2,20})', re.IGNORECASE),
    re.compile(r'\b(ORF[A-Za-z0-9_\-\.]{1,19})\b', re.IGNORECASE)]

# =========================
# Parser
# =========================

def parse_brenda_file() -> Dict[str, dict]:
    """
    Trimmed BRENDA parser (concise, no bloat):
    - protein_name: RN (fallback SN). No REN aliasing.
    - organisms: PR indices -> {species, genus}
    - by_genus routing via indices for SP/NSP, IN, ME, CF, AC, GI(mechanism)
    - corrosion_mechanisms: from GI 'mechanism...' snippets ONLY (no physiological_functions)
    - operational_environmental_factors: OS/GS/SS/PO/PS/PHR/TO/TS/TR/OSS/IS (captured for later; not used here)
    - gene_name: from SY "(gene name ...)" AND enriched from CL lines (without retaining cloning_info)
    - Ignored entirely: reactions (RE), references (RF), applications (AP), engineering (EN),
      structural info (CR/LO/ST/SU), biochemical kinetics/properties (KM/SA/MW/PI),
      purification info (PU), expression (EXP)
    """
    # parse_1 Resolve paths 
    paths = setup_paths()  # must return {'enzyme_brenda': path}
    enzyme_brenda_path = paths['enzyme_brenda']

    if not os.path.exists(enzyme_brenda_path):
        logger.error("BRENDA file not found at path: %s", enzyme_brenda_path)
        return {}

    ec_info: Dict[str, dict] = {}
    current_ec: Optional[str] = None
    in_entry = False

    # parse_2 Continuations
    last_tag: Optional[str] = None
    last_ec_field_list_name: Optional[str] = None
    last_by_genus_targets: Optional[List[str]] = None

    # parse_3 GI buffering
    gi_buffer_text: Optional[str] = None
    gi_buffer_genera: Optional[List[str]] = None

    kept, skipped = 0, 0
    missing_headers = set()
    # parse_4 Helper functions ensure_ec, append_ec_list, append_by_genus, indices_to_genera, flush_gi_buffer
    ###############
    DEFAULT_GENUS_FIELDS = ('substrates','inhibitors','metals','cofactors','activators',
                            'corrosion_mechanisms','reaction_equation','operational_environmental_factors','gene_names')
    def _make_genus_bucket():
        '''Avoids repeating same dict construction'''
        return {k: [] for k in DEFAULT_GENUS_FIELDS}

    def ensure_ec(ec: str):
        if ec not in ec_info:
            ec_info[ec] = {
                'ec_number': ec,
                'protein_name': "",
                'enzyme_class': enzyme_class_from_ec(ec),

                # Core lists kept
                'substrates': [],
                'inhibitors': [],
                'metals': [],
                'cofactors': [],
                'compounds': [],
                'activators': [],
                'corrosion_mechanisms': [],
                'reaction_equation': [],  # from RE lines only
                'operational_environmental_factors': [],
                'brenda_metals': [],
                '_metal_text': [],
                '_mech_text': [], 
                # Gene aggregation
                'gene_name': [],
                # Organisms and per-genus routing
                'organisms': {},  # idx -> {'species','genus'}
                'by_genus': defaultdict(_make_genus_bucket),
            }
    _METAL_PATTERNS = [(re.compile(rf'(?<![A-Za-z0-9]){re.escape(k)}(?![A-Za-z0-9])'), v)
                        for k, v in cs.metal_mapping.items()]
    # parse_5 append_ec_list, append_by_genus, indices_to_genera, flush_gi_buffer
    ###############
    def append_ec_list(ec: str, field: str, value: str):
        if not value:
            return
        v = brenda_sanitize_text(value)
        if not v:
            return
        # Make this resilient: create the list if it doesn't exist
        if field not in ec_info[ec]:
            ec_info[ec][field] = []
        lst = ec_info[ec][field]
        if v not in lst:
            lst.append(v)
    # parse_6 append by genus
    def append_by_genus(ec: str, genera: List[str], field: str, value: str):
        if not value or not genera:
            return
        v = brenda_sanitize_text(value)
        if not v:
            return
        for g in genera:
            gen_bucket = ec_info[ec]['by_genus'][g]
            if field not in gen_bucket:
                gen_bucket[field] = []
            if v not in gen_bucket[field]:
                gen_bucket[field].append(v)
    # parse_7
    def indices_to_genera(ec: str, idxs: Optional[List[str]]) -> List[str]:
        if not idxs:
            return []
        gens, seen = [], set()
        for idx in idxs:
            org = ec_info[ec]['organisms'].get(idx)
            g = org['genus'] if org else None
            if g and g not in seen:
                seen.add(g)
                gens.append(g)
        return gens
    # parse_8 for general information
    def flush_gi_buffer():
        nonlocal gi_buffer_text, gi_buffer_genera
        if not (current_ec and gi_buffer_text):
            gi_buffer_text = None
            gi_buffer_genera = None
            return
        text = brenda_sanitize_text(gi_buffer_text)
        # parse_9
        if text:
            ec_info[current_ec]['_mech_text'].append(text)
        gi_buffer_text = None
        gi_buffer_genera = None
    # parse_12 Main file loop
    ###############
    try:
        with open(enzyme_brenda_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for raw_line in fh:
                # Skip star banners; do not treat as delimiters
                if raw_line.startswith('*' * 20):
                    continue

                line = raw_line.rstrip('\n')
                if not line.strip():
                    continue

                # Entry delimiter
                if line.strip() == "///":
                    flush_gi_buffer()
                    current_ec = None
                    in_entry = False
                    last_tag = None
                    last_ec_field_list_name = None
                    last_by_genus_targets = None
                    continue

                m = TAG_LINE_RE.match(line)
                if not m:
                    # Continuation: extend last list item and GI buffer only
                    if not (in_entry and current_ec and last_tag):
                        continue
                    cont = brenda_sanitize_text(line.strip())
                    if not cont:
                        continue
                    # parse_13 Continuation handling    
                    # Extend EC-wide list item
                    if last_ec_field_list_name:
                        lst = ec_info[current_ec][last_ec_field_list_name]
                        if lst:
                            # only allow continuation if the continuation line is short (e.g., < 100 chars) and the last line isn't already long.
                            if len(cont) < 100 and len(lst[-1]) < 200:
                                lst[-1] = brenda_sanitize_text(f"{lst[-1]} {cont}")
                            else:   
                                lst.append(cont)
                        else:
                            lst.append(cont)

                    # Extend per-genus list item
                    if last_by_genus_targets and last_ec_field_list_name:
                        for g in last_by_genus_targets:
                            gen_bucket = ec_info[current_ec]['by_genus'][g]
                            if last_ec_field_list_name not in gen_bucket:
                                gen_bucket[last_ec_field_list_name] = []
                            lst = gen_bucket[last_ec_field_list_name]
                            if lst:
                                if len(cont) < 100 and len(lst[-1]) < 200:
                                    lst[-1] = brenda_sanitize_text(f"{lst[-1]} {cont}")
                                else:
                                    lst.append(cont)
                            else:
                                lst.append(cont)

                    # Extend GI buffer
                    if last_tag == "GI":
                        gi_buffer_text = f"{gi_buffer_text or ''} {cont}".strip()
                    # Also collect continuation text for metal text-mining (exclude CF)
                    if last_tag in ("RN","SN","SP","NSP","IN","AC","RE","GI") and current_ec:
                        ec_info[current_ec]['_metal_text'].append(cont)
                        ec_info[current_ec]['_mech_text'].append(cont)
                    continue

                # parse_14 New tagged line handling → flush GI if we were buffering
                if last_tag == "GI":
                    flush_gi_buffer()

                tag1, maybe_tag2, value = m.group(1), m.group(2), m.group(3)
                tag1 = tag1.strip()
                val = value.strip()

                # Resolve tag
                if maybe_tag2 and maybe_tag2 in VALID_SHORT_TAGS:
                    tag = maybe_tag2
                elif tag1 in VALID_SHORT_TAGS:
                    tag = tag1
                else:
                    tag = LONG_TO_SHORT.get(tag1)

                if tag in SHORT_TAG_ALIAS:
                    tag = SHORT_TAG_ALIAS[tag]

                if tag is None:
                    # Unknown section header; track and continue (do not reset entry)
                    missing_headers.add(tag1)
                    last_tag = None
                    last_ec_field_list_name = None
                    last_by_genus_targets = None
                    continue

                # parse_15 Start of a new EC entry
                if tag == "ID":
                    norm = normalize_ec_id(val)
                    if not norm:
                        current_ec = None
                        in_entry = False
                        last_tag = None
                        skipped += 1
                        continue
                    current_ec = norm
                    in_entry = True
                    ensure_ec(current_ec)
                    ec_info[current_ec]['enzyme_class'] = enzyme_class_from_ec(current_ec)
                    kept += 1
                    last_tag = tag
                    last_ec_field_list_name = None
                    last_by_genus_targets = None
                    continue

                if not (in_entry and current_ec):
                    continue

                # Reset continuation anchors for this new tagged line
                last_tag = tag
                last_ec_field_list_name = None
                last_by_genus_targets = None

                # Parse_17 indices and value
                idxs, val_wo_idx = extract_indices_and_value(val)
                clean_val = brenda_sanitize_text(val_wo_idx)
                genera = indices_to_genera(current_ec, idxs) if idxs else []
                # collect text for metal text-mining (exclude CF)
                if clean_val and tag in ("RN","SN","SP","NSP","IN","AC","RE","GI"):
                    ec_info[current_ec]['_metal_text'].append(clean_val)
                    ec_info[current_ec]['_mech_text'].append(clean_val)

                # parse_18 Tag handlers (trimmed to requested set)
                if tag == "PR":
                    if not idxs:
                        continue
                    species = brenda_sanitize_text(val_wo_idx)
                    genus = extract_genus(species)
                    for idx in idxs:
                        ec_info[current_ec]['organisms'][idx] = {'species': species, 'genus': genus}

                elif tag == "RN":
                    if clean_val:
                        ec_info[current_ec]['protein_name'] = clean_protein_name(enhanced_clean_protein_name(clean_val))

                elif tag == "SN":
                    if clean_val and not ec_info[current_ec]['protein_name']:
                        ec_info[current_ec]['protein_name'] = clean_protein_name(enhanced_clean_protein_name(clean_val))

                elif tag =="RE":
                        # Extract detailed reaction information
                    if '=' in clean_val and is_valid_text(clean_val):
                        append_ec_list(current_ec, 'reaction_equation', clean_val)

                elif tag in ("SP", "NSP"):
                    if is_valid_text(clean_val):
                        append_ec_list(current_ec, 'substrates', clean_val)
                        append_ec_list(current_ec, 'compounds', clean_val)
                        if genera:
                            append_by_genus(current_ec, genera, 'substrates', clean_val)
                    last_ec_field_list_name = 'substrates'
                    last_by_genus_targets = genera

                elif tag == "IN":
                    if is_valid_text(clean_val):
                        append_ec_list(current_ec, 'inhibitors', clean_val)
                        append_ec_list(current_ec, 'compounds', clean_val)
                        if genera:
                            append_by_genus(current_ec, genera, 'inhibitors', clean_val)
                    last_ec_field_list_name = 'inhibitors'
                    last_by_genus_targets = genera

                elif tag == "ME":
                    if is_valid_text(clean_val):
                        lc = clean_val.lower()
                        metals_found = set()
                        # match mapping keys with safe boundaries to avoid 'fe' in 'female', 'na' in 'nadh'
                        for pat, sym in _METAL_PATTERNS:
                            # word-ish boundary: not alnum on either side of the key
                            if pat.search(lc):
                                metals_found.add(sym)
                        for sym in metals_found:
                            append_ec_list(current_ec, 'metals', sym)
                            if genera:
                                append_by_genus(current_ec, genera, 'metals', sym)
                    last_ec_field_list_name = 'metals'
                    last_by_genus_targets = genera

                elif tag == "CF":
                    if is_valid_text(clean_val):
                        # remove inline citations, angle refs, (...) blocks, concentrations
                        t = CITATION_HASH_RE.sub('', clean_val)
                        t = CITATION_ANGLE_RE.sub('', t)
                        t = re.sub(r'\([^)]*\)', '', t)
                        t = re.sub(r'\b\d[\d.\s]*\s*(?:mM|µM|uM|nM|pM|M)\b', '', t, flags=re.IGNORECASE)
                        t = re.sub(r'\s{2,}', ' ', t).strip()
                        up = t.upper()

                        tok = next((k for k in _COFACTOR_LEXICON if k in up), None)
                        if not tok:
                            tok = t.split(';',1)[0].split(',',1)[0].split(' ',1)[0]

                        append_ec_list(current_ec, 'cofactors', tok)
                        append_ec_list(current_ec, 'compounds', clean_val)
                        if genera:
                            append_by_genus(current_ec, genera, 'cofactors', clean_val)
                        tokens = set(re.split(r'[\s,;()]+', up))
                        hits = [c for c in _COFACTOR_LEXICON if c in tokens]
                        for cof in hits:
                            append_ec_list(current_ec, 'cofactors', cof)
                            # keep compounds comment out if too noisy
                            #append_ec_list(current_ec, 'compounds', clean_val)
                            if genera:
                                append_by_genus(current_ec, genera, 'cofactors', cof)    

                    last_ec_field_list_name = 'cofactors'
                    last_by_genus_targets = genera

                elif tag == "AC":
                    if is_valid_text(clean_val):
                        append_ec_list(current_ec, 'activators', clean_val)
                        append_ec_list(current_ec, 'compounds', clean_val)
                        if genera:
                            append_by_genus(current_ec, genera, 'activators', clean_val)
                    last_ec_field_list_name = 'activators'
                    last_by_genus_targets = genera

                elif tag == "SY":
                    # Capture only explicit "(gene name ...)" annotations.
                    if is_valid_text(val_wo_idx):
                        # Use a lighter clean to preserve parentheses content needed for the pattern
                        sy_raw = CITATION_HASH_RE.sub('', val_wo_idx)
                        sy_raw = CITATION_ANGLE_RE.sub('', sy_raw).strip()
                        m_gene = re.search(r'\(gene name[:\s]*([A-Za-z0-9._-]{2,20})\)', sy_raw, flags=re.IGNORECASE)
                        gene = m_gene.group(1) if m_gene else None

                        if gene:
                            if genera:
                                for g in genera:
                                    gb = ec_info[current_ec]['by_genus'][g]
                                    if 'gene_names' not in gb:
                                        gb['gene_names'] = []
                                    if gene not in gb['gene_names']:
                                        gb['gene_names'].append(gene)
                            else:
                                # fall back to top-level if no indices
                                if gene not in ec_info[current_ec]['gene_name']:
                                    ec_info[current_ec]['gene_name'].append(gene)

                elif tag == "CL":
                    # Enrich gene_name from CL text without keeping cloning_info
                    if is_valid_text(clean_val) and genera:
                        found = set()
                        for pat in CL_GENE_PATTERNS:
                            for m in pat.finditer(clean_val):
                                token = m.group(1).strip()
                                if _is_gene_token(token):
                                    found.add(token)
                        if found:
                            for g in genera:
                                gb = ec_info[current_ec]['by_genus'][g]
                                if 'gene_names' not in gb:
                                    gb['gene_names'] = []
                                for token in found:
                                    if token not in gb['gene_names']:
                                        gb['gene_names'].append(token)

                elif tag == "GI":
                    # Buffer for mechanism-only extraction
                    gi_buffer_text = clean_val
                    gi_buffer_genera = genera

                elif tag in ENV_TAGS:
                    if is_valid_text(clean_val):
                        labeled = f"{tag}: {clean_val}"
                        append_ec_list(current_ec, 'operational_environmental_factors', labeled)
                        if genera:
                            append_by_genus(current_ec, genera, 'operational_environmental_factors', labeled)
                    last_ec_field_list_name = 'operational_environmental_factors'
                    last_by_genus_targets = genera

                # Everything else is ignored to keep output concise
                elif tag in ("RE", "RF", "MW", "PI", "AP", "EN", "CR", "LO", "ST", "KM", "SA", "PU", "EXP", "RT", "BR"):
                    last_ec_field_list_name = None
                    last_by_genus_targets = None

                else:
                    # Unknown/unused tag
                    last_ec_field_list_name = None
                    last_by_genus_targets = None
        # parse_19 Final flush at EOF
        flush_gi_buffer()

        # # parse_20 deduplicate reactions dict fromkeys and sort by lenght prioritize equation with = sign                  
        for ec, info in ec_info.items():                    
            reactions = info.get('reaction_equation', [])
            # Sort reactions: those containing '=' first, keep top 10 by length reaction[:10]
            sorted_reactions = sorted(reactions, key=lambda s: (('=' not in s), -len(s)))
            # Deduplicate while preserving order
            deduped_reactions = list(dict.fromkeys(sorted_reactions))[:10]
            info['reaction_equation'] = deduped_reactions

        # parse_21 Finalize: dedupe lists and aggregate gene names to top-level
        for ec, info in ec_info.items():
            # Build top-level gene_name from per-genus gene_names
            agg_gene = set(info.get('gene_name', []))

            # Clean and dedupe per-genus
            by_genus_clean = {}
            for g, buckets in info['by_genus'].items():
                cleaned = {k: list(dict.fromkeys(v)) for k, v in buckets.items()}
                by_genus_clean[g] = cleaned
                for gn in cleaned.get('gene_names', []):
                    agg_gene.add(gn)

            info['by_genus'] = by_genus_clean
            info['gene_name'] = sorted(agg_gene)

            # parse_21b: Mechanism classification from buffered text - EXTRACT CHILD TERMS
            mech_corpus = ' '.join(info.get('_mech_text', []))

            mechanisms_processor = TermProcessor(mechanisms_dict)
            mech_hits = mechanisms_processor.find_all_matches(mech_corpus)

            # FIXED: Extract child terms instead of subcategory keys
            child_terms_found = []
            subcategories_involved = []
            for subcategory, child_terms in mech_hits.items():
                if child_terms:  # Only if we actually found terms
                    subcategories_involved.append(subcategory)
                    child_terms_found.extend(child_terms)

            # Store child terms (for text mining evidence)
            info['corrosion_mechanisms'] = sorted(set(child_terms_found))

            # Optionally store subcategories too (for analysis/debugging)
            info['corrosion_mechanism_categories'] = sorted(set(subcategories_involved))


            info.pop('_mech_text', None)   # drop buffer

            # parse 22 Finalize: metal text-mining#brenda_metals = ME ∪ mined-text (exclude CF)
            corpus = ' '.join(info.get('_metal_text', [])).lower()
            mined_syms = set()
            if corpus:
                for pat, sym in _METAL_PATTERNS: #metal processor incorporated
                    if pat.search(corpus):
                        mined_syms.add(sym)  

            # Symbols already captured from ME tag (parser put canonical symbols there)
            me_syms = set(info.get('metals', []))
            # Whitelist to canonical symbols only, as defined by the mapping values
            allowed_symbols = set(cs.metal_mapping.values())
            merged = [s for s in dict.fromkeys(list(me_syms | mined_syms)) if s in allowed_symbols]

            # parse 23 merging of all the metals from brenda
            info['brenda_metals'] = merged # corresponds to all the metals that are found in brenda
            info.pop('_metal_text', None)
            info.pop('metals', None)

            # Deduplicate EC-wide lists
            for fld, val in list(info.items()):
                if fld == 'by_genus':
                    continue
                if isinstance(val, list):
                    info[fld] = list(dict.fromkeys(val))

        # parse_24 Stats   
        ec_with_org = sum(1 for v in ec_info.values() if v['organisms'])
        ec_with_name = sum(1 for v in ec_info.values() if v['protein_name'])
        ec_with_metals = sum(1 for v in ec_info.values() if v['brenda_metals'])
        ec_with_mech = sum(1 for v in ec_info.values() if v['corrosion_mechanisms'])
        ec_with_compounds = sum(1 for v in ec_info.values() if v['compounds'])
        ec_with_reaction = sum(1 for v in ec_info.values() if v['reaction_equation'])
        ec_with_env = sum(1 for v in ec_info.values() if v['operational_environmental_factors'])
        ec_with_cofactors = sum(1 for v in ec_info.values() if v['cofactors'])
        ec_with_activators = sum(1 for v in ec_info.values() if v['activators'])
        ec_with_gene = sum(1 for v in ec_info.values() if v['gene_name'])

        total = max(kept, 1)
        logger.info("BRENDA parsing (trimmed): kept=%d, skipped=%d | organisms=%d (%.1f%%), names=%d (%.1f%%), brenda_metals=%d (%.1f%%), mechanisms=%d (%.1f%%), env=%d (%.1f%%), cofactors=%d (%.1f%%), activators=%d (%.1f%%), gene=%d (%.1f%%), reactions=%d (%.1f%%)", kept, skipped, ec_with_org, 100*ec_with_org/total, ec_with_name, 100*ec_with_name/total, ec_with_metals, 100*ec_with_metals/total, ec_with_mech, 100*ec_with_mech/total, ec_with_env, 100*ec_with_env/total, ec_with_cofactors, 100*ec_with_cofactors/total, ec_with_activators, 100*ec_with_activators/total, ec_with_gene, 100*ec_with_gene/total, ec_with_reaction, 100*ec_with_reaction/total) # ec_with_pathway, 100*ec_with_pathway/total,

        if missing_headers:
            logger.debug("Unmapped section headers encountered (ignored): %s", sorted(missing_headers))

        return ec_info

    except Exception as e:
        logger.exception("Error parsing BRENDA file: %s", e)

    return {}

#brenda_data = parse_brenda_file()


# _____________________
# ### Process Brenda Data
# _______________

# In[9]:


def process_brenda_data(brenda_data, processors, cs, config):
    """
    Redesign: build brenda_en with full v3 scoring and NO field loss.
    - Preserves ALL original parser fields (no pops/deletes).
    - Adds/overwrites the important scoring fields so they are always present.
    - Uses parser's `brenda_metals` as-is (name preserved).
    """
    if not isinstance(brenda_data, dict):
        logger.error("Invalid brenda_data; expected dict")
        return {}

    REQUIRED_DEFAULTS = {
        'functional_categories': [],
        'functional_child_terms': [],
        'functional_matches_detailed': {},
        'fc_cooccurrence_synergy_hit': False,
        'synergy_child_terms_found': [],
        'synergy_categories_involved': [],
        'synergy_description': '',
        'synergy_type': 'none',
        'functional_score': 0.0,
        'overall_functional_score': 0.0,
        'metal_score': 0.0,
        'overall_metal_score': 0.0,
        'synergy_score': 0.0,
        'overall_synergy_score': 0.0,
        'corrosion_relevance_score': 0.0,
        'metal_count': 0,
        'metal_categories_detected': []
    }

    processed = {}

    # local config passthrough—use given config as-is
    local_cfg = config if isinstance(config, cs.ScoringConfig) else cs.ScoringConfig()

    for ec_number, rec in brenda_data.items():
        if not isinstance(rec, dict):
            continue

        # --- 1) Preserve ALL original parser fields (shallow copy)
        out = dict(rec)

        # --- 2) Build analysis text (same evidence recipe as before)
        def _is_valid_equation(eq):
            if not eq or len(eq) < 10: return False
            if eq.startswith('+') or '...' in eq or eq.count('=') > 1: return False
            return any(tok in eq for tok in ['+', '=', '→', '⇌'])

        subs  = ' '.join(out.get('substrates', []) or [])
        cof   = ' '.join(out.get('cofactors', []) or [])
        inhib = ' '.join(out.get('inhibitors', []) or [])
        comp  = ' '.join(out.get('compounds', []) or [])
        env   = ' '.join(out.get('operational_environmental_factors', []) or [])
        rxn   = ' '.join([r for r in (out.get('reaction_equation', []) or []) if _is_valid_equation(r)])
        prot  = str(out.get('protein_name', '') or '')
        mechs = ' '.join(out.get('corrosion_mechanisms', []) or [])  # child terms from parser

        parts = [subs, cof, inhib, comp, env, rxn, prot, mechs]
        analysis_text = ' '.join(p for p in parts if p).strip()

        # Metals from parser (keep name/content)
        brenda_metals = list(dict.fromkeys(out.get('brenda_metals', []) or []))

        # --- 3) Score with v3 (processors dict)
        try:
            results = cs.calculate_overall_scores(
                text=analysis_text,
                processors=processors,     # {'fc_processor':..., 'metal_processor':..., 'synergy_processor':...}
                config=local_cfg,
                brenda_metals=brenda_metals
            )
        except Exception as e:
            logger.error(f"V3 scoring failed for EC {ec_number}: {e}")
            results = {}

        # --- 4) Normalize synergy child terms to plain list
        syn_terms = results.get('synergy_child_terms_found', [])
        if isinstance(syn_terms, set):
            syn_terms = sorted(syn_terms)
        elif isinstance(syn_terms, list):
            syn_terms = list(dict.fromkeys(syn_terms))
        else:
            syn_terms = []

        # --- 5) Inject/overwrite ONLY the important fields (no deletions of originals)
        out.update({
            # keep parser’s identifiers
            'ec_number': out.get('ec_number', ec_number),
            'protein_name': out.get('protein_name', ''),
            'enzyme_class': out.get('enzyme_class', ''),

            # keep parser metals name exactly
            'brenda_metals': brenda_metals,

            # functional mining outputs
            'functional_categories': results.get('functional_categories_detected', REQUIRED_DEFAULTS['functional_categories']),
            'functional_child_terms': results.get('functional_child_terms', REQUIRED_DEFAULTS['functional_child_terms']),
            'functional_matches_detailed': results.get('functional_matches_detailed', REQUIRED_DEFAULTS['functional_matches_detailed']),

            # synergy outputs
            'fc_cooccurrence_synergy_hit': results.get('fc_cooccurrence_synergy_hit', REQUIRED_DEFAULTS['fc_cooccurrence_synergy_hit']),
            'synergy_child_terms_found': syn_terms,
            'synergy_categories_involved': results.get('synergy_categories_involved', REQUIRED_DEFAULTS['synergy_categories_involved']),
            'synergy_description': results.get('synergy_description', REQUIRED_DEFAULTS['synergy_description']),
            'synergy_type': results.get('synergy_type', REQUIRED_DEFAULTS['synergy_type']),

            # scores
            'functional_score': results.get('functional_score', REQUIRED_DEFAULTS['functional_score']),
            'overall_functional_score': results.get('overall_functional_score', REQUIRED_DEFAULTS['overall_functional_score']),
            'metal_score': results.get('metal_score', REQUIRED_DEFAULTS['metal_score']),
            'overall_metal_score': results.get('overall_metal_score', REQUIRED_DEFAULTS['overall_metal_score']),
            'synergy_score': results.get('synergy_score', REQUIRED_DEFAULTS['synergy_score']),
            'overall_synergy_score': results.get('overall_synergy_score', REQUIRED_DEFAULTS['overall_synergy_score']),
            'corrosion_relevance_score': results.get('corrosion_relevance_score', REQUIRED_DEFAULTS['corrosion_relevance_score']),

            # metal detection metadata
            'metal_count': results.get('metal_count', REQUIRED_DEFAULTS['metal_count']),
            'metal_categories_detected': results.get('metal_categories_detected', REQUIRED_DEFAULTS['metal_categories_detected']),
        })

        # NO pops. NO deletions. Everything from parser stays; important fields are guaranteed.
        processed[ec_number] = out

    logger.info("Processed %d enzyme records with v3 scoring (preserve-all).", len(processed))
    return processed


# In[10]:


'''# process data sample
{'1.1.1.1': {'ec_number': '1.1.1.1',
  'protein_name': 'alcohol-dehydrogenase',
  'enzyme_class': 'oxidoreductases',
  'substrates': ['acetaldehyde + NADH + H+ = ethanol + NAD+ ( cells with an extra copy of ADH1 display chronological life-span extension. Antioxidant enzymes are induced in 2xADH1 cells. Strains carrying an extra ADH1 copy show extended replicative life span and',
   'increased Sir2p activity )',
   'methylglyoxal + NADH + H+ = acetol + NAD+',
'''


# __output brenda_en__
# Temporal brenda sample creation  and their associated test for testing purposes

# In[11]:


'''# Temporal brenda sample creation for testing purposes # 20 min at 11 GB and 3 cores 
import random
def create_brenda_sample(brenda_data, sample_size=1000):
    """Create a small sample of BRENDA data for testing"""
    if not brenda_data:
        return {}

    # Taking first N entries
    sample_keys = list(brenda_data.keys())[:sample_size]
    return {k: brenda_data[k] for k in sample_keys}
# Create a sample of BRENDA data 6 min
brenda_data = create_brenda_sample(brenda_data, sample_size=1500
#brenda_en = process_brenda_data(brenda_data, processors, cs, config)
# #brenda_en_df = pd.DataFrame.from_dict(brenda_en, orient='index')
#brenda_en_df['brenda_metals'].explode().unique()'''


# ___________________
# ### Enzyme names
# ____________________
# The database containing enzyme names and EC numbers
# 
# wget https://www.enzyme-database.org/downloads/enzyme-database.sql.gz

# In[12]:


def read_enzyme_names(unique_ecs_to_filter=None):
    """
    Read and parse enzyme file to get EC numbers and their names.
    - Normalizes EC keys (strip 'EC ', accept hyphens)
    - Applies enhanced_clean_protein_name to each name
    - Optionally filters to a set of (normalized) ECs
    """
    paths = setup_paths()
    enzyme_path = paths['enzyme']

    # Normalize the filter set up front
    ec_filter_norm = None
    if unique_ecs_to_filter is not None:
        ec_filter_norm = {normalize_ec_id(ec) for ec in unique_ecs_to_filter}
        ec_filter_norm.discard(None)

    ec_to_names = {}
    with open(enzyme_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t', 1)
            if len(parts) < 2:
                continue
            raw_ec, raw_names = parts[0], parts[1]
            ec_norm = normalize_ec_id(raw_ec)
            if not ec_norm:
                continue
            if ec_filter_norm is not None and ec_norm not in ec_filter_norm:
                continue
            if raw_names.lower().startswith("transferred to"):
                continue
            # Robust split on semicolons with optional spaces
            names = [n.strip() for n in raw_names.split(';') if n.strip()]
            cleaned = [enhanced_clean_protein_name(n) for n in names]
            cleaned = [clean_protein_name(n)for n in cleaned]
            ec_to_names[ec_norm] = cleaned

    return ec_to_names
#ec_to_names = read_enzyme_names(ECcontri_Uniprot["ec"].unique())# THis is just for testing purposes this is to be silenced

'''ec_to_names
{'1.1.1.1': ['alcohol-dehydrogenase',
  'aldehyde-reductase',
  'adh',
  'alcohol-dehydrogenase',
  'aliphatic alcohol-dehydrogenase',
  'ethanol-dehydrogenase',
  'nad-dependent alcohol-dehydrogenase',
  'nad-specific aromatic alcohol-dehydrogenase',
  'nadh-alcohol-dehydrogenase',
  'nadh-aldehyde-dehydrogenase',
  'primary alcohol-dehydrogenase',
  'yeast alcohol-dehydrogenase'],'''
#ec_to_names = read_enzyme_names()


# ___________
# ### Enzyme Class
# ________________
# The enzyme classification system (text-based hierarchy).

# In[13]:


def read_enzyme_class():
    """Parse enzclass.txt exactly as provided (handles spaces inside the code)."""
    paths = setup_paths()
    ec_file_path = paths['enzyme_class']
    enzyme_class = {}
    with open(ec_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or not s[0].isdigit():
                continue
            # split on 2+ spaces → [code-part, description]
            parts = re.split(r'\s{2,}', s, maxsplit=1)
            if len(parts) != 2:
                continue
            code_raw, desc = parts[0], parts[1].strip()
            ec_id = code_raw.replace(' ', '')  # "1. 1. 1.-" → "1.1.1.-", "1. -. -.-" → "1.-.-.-"
            if ec_id:
                enzyme_class[ec_id] = desc
    return enzyme_class
#enzyme_class =read_enzyme_class()
#enzyme_class


# _________________________________________
# ### Reaction Data
# ____________________________________________
#  Reaction-level information.
# 
# !wget -c "ftp://ftp.genome.jp/pub/kegg/reaction/reaction.tar.gz" # This reaction data wont be used in the creation of the db

# In[14]:


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
'''
{'R00001': {'name': 'polyphosphate polyphosphohydrolase',
  'equation': 'Polyphosphate + n H2O <=> (n+1) Oligophosphate'},
 'R00002': {'name': 'reduced ferredoxin:dinitrogen oxidoreductase (ATP-hydrolysing)',
  'equation': '16 ATP + 16 H2O + 8 Reduced ferredoxin <=> 8 e- + 16 Orthophosphate + 16 ADP + 8 Oxidized ferredoxin'},
 'R00004': {'name': 'diphosphate phosphohydrolase'}'''
#reaction_info = read_reaction_data()


# ________________________________
# ### Create reaction mapping
# _____________________

# In[15]:


def create_ec_to_reaction_mapping():  
    ec_to_names = read_enzyme_names()
    reaction_info = read_reaction_data()
    # join the reaction_brenda on ec numbers
    reaction_brenda = pd.DataFrame(list(reaction_info.values()), index=reaction_info.keys())

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
#ec_to_reaction = create_ec_to_reaction_mapping()


# ______________________________
# ### reaction curated
# ____________________________

# In[16]:


def build_reactions_db_curated(record, brenda_en, ec_to_reaction, reaction_info,
                               max_items: int = 10, max_str_len: int = 220):
    """
    Mutates `record` to add:
      record['reaction_db'] : list[str] (prefer DB "name: equation"; fallback BRENDA equations)

    - Dedupe (preserve order)
    - Priority: entries containing '=' first, then length desc
    - Cap to `max_items`
    - Truncate each string to `max_str_len` chars
    """
    if not isinstance(record, dict):
        raise TypeError("record must be a dict with 'ec_number'")
    ec = normalize_ec_id(record.get('ec_number', '') or '')
    if not ec:
        record['reaction_db'] = []
        return record

    out = []

    # 1) DB reactions via ec_to_reaction → reaction_info
    for rid in (ec_to_reaction.get(ec, []) or []):
        info = reaction_info.get(rid, {}) or {}
        nm = (info.get('name') or '').strip()
        eq = (info.get('equation') or '').strip()
        if nm and eq:
            s = f"{nm}: {eq}"
        elif nm:
            s = nm
        elif eq:
            s = eq
        else:
            continue
        s = ' '.join(s.split())  # compress whitespace
        if len(s) > max_str_len:
            s = s[:max_str_len - 1] + '…'
        out.append(s)

    # 2) Fallback — BRENDA reaction_equation (only if DB list is empty)
    if not out and brenda_en and ec in brenda_en:
        for s in (brenda_en.get(ec,  {}).get('reaction_equation', []) or []):
            if not isinstance(s, str):
                continue
            s = ' '.join(s.split()).strip()
            if not s:
                continue
            if len(s) > max_str_len:
                s = s[:max_str_len - 1] + '…'
            out.append(s)

    # 3) Dedupe, prioritize '=', then cap
    seen = set()
    dedup = []
    for s in out:
        if s not in seen:
            seen.add(s)
            dedup.append(s)

    dedup.sort(key=lambda x: (('=' not in x), -len(x)))
    record['reaction_db'] = dedup[:max_items]
    return record
#reaction_db # single EC test
#rec = {"ec_number": "1.1.1.1"}
#rec = build_reactions_db_curated(rec, brenda_en_co, ec_to_reaction, reaction_info)
# rec['reaction_db'] -> compact list[str]


# ____________________
# ### ec_to_pathway mapping
# _______________________
# 
# wget -O ec_pathway.list "https://rest.kegg.jp/link/pathway/ec"

# In[17]:


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
                    # normalise ec
                    normalized_ec = normalize_ec_id(ec_id)

                    pathway_id = parts[1].replace('path:', '')

                    if not normalized_ec:
                        continue
                    if normalized_ec not in ec_to_pathway:
                        ec_to_pathway[normalized_ec] = []
                    ec_to_pathway[normalized_ec].append(pathway_id)

        print(f"Loaded pathway mappings for {len(ec_to_pathway)} EC numbers")
    except Exception as e:
        print(f"Error reading EC-pathway mapping: {e}")
        return {}

    return ec_to_pathway
#ec_pathway_mapping = read_ec_pathway_mapping()
'''
ec_pathway_mapping
{'1.1.1.1': ['map00010',
  'ec00010',
  'map00071',
  'ec00071',
  'map00260',
  'ec00260',
  'map00350',
  'ec00350'''


# _________________
# ### Ko Database
# ____________________
#  A new variable mapping  KO numbers to EC numbers from KEGG KO file.
# rsync -avz rsync://rest.kegg.jp/kegg/ko/ . 

# In[18]:


def read_ko_data() -> dict:
    """Read and parse KEGG KO file into a dict of {Kxxxxx: {definition, pathway}}"""
    paths = setup_paths()
    ko_file_path = paths['ko'] # from 'ko': db_dir / 'kegg_ko.txt',

    ec_pattern = re.compile(r'\[EC:([0-9\.\-\s]+)\]')

    ko_info = {}
    with open(ko_file_path, 'r') as f:
        for line in f:
            if line.startswith('K'):
                parts = line.strip().split('\t')
                ko_id = parts[0]
                raw_definition = parts[1] if len(parts) > 1 else ''
                # keep raw for EC extraction
                definition_raw = raw_definition
                # also keep a cleaned version without EC tokens 
                definition = strip_all_ec_tokens(raw_definition)

                # extract ALL EC numbers from raw definition
                ec_blob = ec_pattern.findall(definition_raw)  # list like ["1.1.1.1", "3.2.-.-"]
                ec_numbers = []
                for blob in ec_blob:
                    # blobs sometimes have space-separated ECs within the brackets
                    for token in blob.split():
                        ec_norm = normalize_ec_id(token)
                        if ec_norm:
                            ec_numbers.append(ec_norm)
                path_ids = re.findall(r'\bmap\d{5}\b', definition_raw)  # also find any inline mapXXXXX
                ko_info[ko_id] = {
                    'definition_raw': definition_raw,
                    'definition': definition,
                    'pathway': sorted(set(path_ids)),
                    'ec_numbers': sorted(set(ec_numbers)),
                }

    return ko_info
#ko_ec = read_ko_data()
'''
ko_ec sample 
{'K00001': {'definition_raw': 'E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]',
  'definition': 'E, adh; alcohol dehydrogenase [EC:]',
  'pathway': '',
  'ec_numbers': ['1.1.1.1']},
 'K00002': {'definition_raw': 'AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]',
  'definition': 'AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:]',
  'pathway': '','''


# _______________
# ### ko Hierarchi
# ________________
# The KEGG Orthology system is structured hierarchically into broad categories (A), subcategories (B), pathways (C), and then individual KO entries (D) which have functional annotations and often include enzyme commission (EC) numbers(ncbi)
# Hierarchy of KO numbers (helps in pathway mapping).The BRITE hierarchy (ko00001.keg) is part of KEGG's comprehensive classification system representing various biological processes and protein families. B: Subcategories or BRITE categories - include "Protein families: metabolism" or "Protein families: signaling and cellular processes" which group functionally related protein families.
# The strategy to use here is to chose only from the categories A and B, Table 1. but also only allow some categories, later the inclusion factors are given in table 2
# Table 1. Kegg Ortology levels A and B
# 
# |Level A	Example |Level B (Subcategory)|
# |--|--|
# |Metabolism |	Nucleotide metabolism, Amino acid metabolism, Carbohydrate metabolism, Energy metabolism|
# |Genetic Information Processing|	Transcription, Translation, Replication and repair, Chromosome|
# |Environmental Information Processing|	Membrane transport, Signaling molecules and interaction|
# |Cellular Processes|	Cell motility, Cell growth and death, Cellular community eukaryotes|
# |Organismal Systems|	Immune system, Nervous system, Endocrine system|
# |Human Diseases	|Infectious disease (bacterial, parasitic), Cancer, Neurodegenerative disease|
# |BRITE Hierarchies|	Protein families: metabolism, Protein families: signaling and cellular processes|
# 
# Table 2. Kegg Ortology allowed Categories taken as Hierarchi for the present study
# 
# |Category Type	|Reason for inclusion|
# |--|--|
# |Protein families: metabolism|	Core metabolic enzymes linked to microbial life|
# |Xenobiotics biodegradation and metabolism|	Microbial pathways for breaking down industrial chemicals, pollutants|
# |Transporters|	Microbial substrate uptake and efflux, key to biofilm survival
# |Signal transduction|	Microbial response to environmental stresses
# |Cellular community - prokaryotes|	Biofilm formation and quorum sensing
# |Energy metabolism|	Core energy pathways, important in microbial corrosion
# |Metabolism of cofactors and vitamins|	Support enzymatic functions, influence microbial activity
# |Protein families: signaling and cellular processes	|Immune evasion, stress response—may influence persistence
# |Aging and stress response pathways|	Microbial adaptation to hostile conditions
# |Infectious disease (bacterial & parasitic)|	Possible indicators of bacterial species present
# 

# In[19]:


def read_ko_hierarchy():
    paths = setup_paths()
    ko_path = paths['ko_hierarchy']   # from db_dir / 'ko_hierarchy.txt',
    allowed_categories = [
    'Protein families: metabolism', 'Metabolism of cofactors and vitamins', 'Xenobiotics biodegradation and metabolism',
    'Cellular community - prokaryotes','Signal transduction','Protein families: signaling and cellular processes',
    'Transport and catabolism','Energy metabolism','Nucleotide metabolism','Amino acid metabolism',
    'Glycan biosynthesis and metabolism','Carbohydrate metabolism', 'Lipid metabolism','Metabolism of terpenoids and polyketides',
    'Biosynthesis of other secondary metabolites',
    ]

    all_hierarchy = {
        'A': {},  # Top level Metabolism, Genetic, Environmental, etc.
        'B': {},  # ko Category detailed subcategories (Nucleotide metabolism, Signal transduction, etc.).
        'C': {},  # Pathway specific biochemical or signaling pathways within subcategories.
        'D': {}   # KO/Enzyme  Individual orthologous genes, proteins, enzymes, and their functional annotations.
    } 

    current = {'A': None, 'B': None, 'C': None, 'D': None}
    rows=[]

    with open(ko_path, 'r') as f:
        for line in f:
            if line.startswith('B'):
                parts = line.strip().split()
                b_id = parts[1]
                hierarchy = ' '.join(parts[2:])
                if hierarchy in allowed_categories:
                    current['B'] = b_id
                    current['B_name'] = hierarchy
                else:
                    current['B'] = None
                    current['B_name'] = None

            elif line.startswith('C') and current['B'] is not None:
                parts = line.strip().split()
                c_id = parts[1]
                pathway = ' '.join(parts[2:])
                if '[PATH:' in pathway:
                    pathway = pathway.split('[PATH:')[0].strip()
                current['C'] = c_id
                current['C_name'] = pathway

            elif line.startswith('D') and current['C'] is not None and current['B'] is not None:
                parts = line.strip().split()
                ko_id = parts[1]
                protein_name = ' '.join(parts[2:])
                ec_number = []
                if '[EC:' in protein_name:
                    ec_part = protein_name.split('[EC:')[1].split(']')[0]
                    ec_number = ec_part.split()
                    protein_name = protein_name.split('[EC:')[0].strip()
                ec = ';'.join(ec_number) if ec_number else None

                rows.append({
                        'ko': ko_id,
                        'hierarchy': current['B_name'],
                        'protein_name': protein_name,
                        'ec': ec,
                        'pathway': current['C_name']
                    })

    df_hierarchy = pd.DataFrame(rows)
    return df_hierarchy
hierarchy_all =  read_ko_hierarchy()

'''{'A': {'Metabolism': '',
  'Genetic': 'Information Processing',
  'Environmental': 'Information Processing',
  'Cellular': 'Processes',
  'Organismal': 'Systems',
  'Human': 'Diseases',
  'Brite': 'Hierarchies',
  'Not': 'Included in Pathway or Brite'},
 'B': {'09101': {'name': 'Carbohydrate metabolism', 'parent': 'Metabolism'},
  '09102': {'name': 'Energy metabolism', 'parent': 'Metabolism'},
  '09103': {'name': 'Lipid metabolism', 'parent': 'Metabolism'},
  '09104': {'name': 'Nucleotide metabolism', 'parent': 'Metabolism'},
  '09105': {'name': 'Amino acid metabolism', 'parent': 'Metabolism'}'''


# I decided against getting all the KEGG ortogonal classification because it was too much and it was necesary just the broad picture, so the Brite database was instead retrieved. The brite database had also more EC and protein_names retrievable from the eccontri.

# In[20]:


paths = setup_paths()
brite_path = paths['ko_brite'] # from 'ko_brite': db_dir / 'ko00001.keg',

allowed_categories = [
    'Protein families: metabolism', 'Metabolism of cofactors and vitamins', 'Xenobiotics biodegradation and metabolism',
    'Cellular community - prokaryotes','Signal transduction','Protein families: signaling and cellular processes',
    'Transport and catabolism','Energy metabolism','Nucleotide metabolism','Amino acid metabolism','Infectious disease: bacterial',
    'Glycan biosynthesis and metabolism','Carbohydrate metabolism', 'Lipid metabolism','Metabolism of terpenoids and polyketides',
    'Biosynthesis of other secondary metabolites', 'Metabolism of other amino acids', 'Drug resistance: antimicrobial','Environmental adaptation'
    ]
def parse_ko_brite_filtered(brite_path: str, allowed_b_names: Iterable[str]):
    """
    Parse brite database and gets level B as hierarchy, C as pathways and D as protein_name.
    Builds a flat DataFrame with columns: ko, hierarchy, protein_name, ec, pathway.
    """

    # "D      K08041  ADCY1; adenylate cyclase 1 [EC:4.6.1.1]"
    pattern = re.compile(r'^D\s+(\w+)\s+([^\[]+)(?:\[EC:(.*?)\])?')

    brite_rows = []
    with open(brite_path, 'r', encoding='utf-8') as f:
        current_category = None
        current_pathway = None
        for line in f:
            if line.startswith('B'):
                parts = line.strip().split()
                b_id = parts[1] if len(parts) > 2 else None
                hierarchy = ' '.join(parts[2:]) if len(parts) > 2 else ""
                current_category = hierarchy if hierarchy in allowed_b_names else None

            elif line.startswith('C'):
                payload = line[1:].strip()  # drop leading "C"
                payload = re.sub(r'\[PATH:[^\]]+\]', '', payload).strip()
                current_pathway = re.sub(r'^\d+\s+', '', payload)

            elif line.startswith('D'):
                m = pattern.match(line)
                if m and current_category:
                    ko = m.group(1)
                    names = m.group(2).strip()  # "ADCY1; adenylate cyclase 1"
                    ec = m.group(3).strip() if m.group(3) else ""
                    for name in re.split(r'[;,/]', names):
                        brite_rows.append({
                            'ko': ko,
                            'protein_name': name.strip(),
                            'ec_number': ec,
                            'hierarchy': current_category,
                            'pathway': current_pathway
                        })
    hierarchy_brite = pd.DataFrame(brite_rows)
    #hierarchy_brite["protein_name"] = hierarchy_brite["protein_name"].apply(clean_protein_name)
    hierarchy_brite["protein_name"] = hierarchy_brite["protein_name"].apply(enhanced_clean_protein_name)
    return hierarchy_brite

#hierarchy_brite = parse_ko_brite_filtered(brite_path, allowed_categories)

'''    ko protein_name       ec                hierarchy  \
0  K00844           hk  2.7.1.1  Carbohydrate metabolism   
1  K00844   hexokinase  2.7.1.1  Carbohydrate metabolism   
2  K12407          gck  2.7.1.2  Carbohydrate metabolism   
3  K12407  glucokinase  2.7.1.2  Carbohydrate metabolism '''


# ___________________________
# ### Cleaning pathways function
# ______________________________
#  The following function applies the correct pandas series logic Series.apply behavior and regular expression design (Pandas documentation, 2024).
# For regex practice
# Friedl, J. E. F. (2006). Mastering Regular Expressions (3rd ed.). O'Reilly Media.
# 

# In[21]:


def clean_pathway_strings(pathway_series: pd.Series) -> pd.Series:
    """
    Clean pathway strings by removing common headers and database references that don't represent specific pathways.
    Parameters:pathway_series : pandas.Series containing pathway strings with pathways separated by semicolons
    Returns:  pandas.Series with cleaned pathway strings
    """
    # Terms to remove - add or remove based on specific needs
    terms_to_remove = [
        'Enzymes with EC numbers', # anotations mistake
        'proteins [BR:ko00194]', # anotations mistake
        #'Metabolic pathways', # to broad term
        'Other Metabolic Processes',
        #'Biosynthesis of secondary metabolites', # to broad term
        'Microbial metabolism in diverse environments','Microbial metabolism in diverse environments', 'Metabolic pathways', # to broad term
        'Exosome \[BR:ko04147\]', # mostly eucariotic
        'photosynthesis', 'Photosynthesis', 'Photosynthesis  proteins [BR:ko00194]', # mostly plants
        'Photosynthesis proteins', # mostly plants
    ]

    # Compile regex pattern
    pattern = re.compile(
        '|'.join([re.escape(term) for term in terms_to_remove]),
        flags=re.IGNORECASE
    )

    def clean_single_entry(entry: str) -> str:
        if pd.isna(entry) or not isinstance(entry, str) or entry.strip() == '[]':
            return ''

        # Remove matched unwanted terms
        cleaned = pattern.sub('', entry)

        # Clean up formatting issues
        cleaned = pattern.sub('', entry)
        cleaned = re.sub(r';\s*;', ';', cleaned)
        cleaned = re.sub(r'[;\s]+$', '', cleaned)
        cleaned = re.sub(r'^[;\s]+', '', cleaned)

        return cleaned

    return pathway_series.apply(clean_single_entry)


# ______________________
# ### Pathway Database
# ________________
# Chemical compounds database.wget https://biocyc.org/download.shtml. wget https://www.brenda-enzymes.org/download.php
# 

# In[22]:


def read_pathway_data():
    paths = setup_paths()
    pathway_path = paths['pathway']

    pathway_data = {}
    with open(pathway_path, 'r') as f:
          for line in f:
              parts = line.strip().split('\t')
              if len(parts) >= 2:
                  pathway_id = parts[0]
                  pathway_name = parts[1]
                  pathway_data[pathway_id] = pathway_name
    return pathway_data
#pathway_data =read_pathway_data()
'''pathway_data
{'map01100': 'Metabolic pathways',
 'map01110': 'Biosynthesis of secondary metabolites',
 'map01120': 'Microbial metabolism in diverse environments',
 'map01200': 'Carbon metabolism',
 'map01210': '2-Oxocarboxylic acid metabolism',
 'map01212': 'Fatty acid metabolism','''


# ______________________
# ### Mapping ko to ec
# ______________________

# In[23]:


def build_ec_to_ko_map(ko_ec:dict) -> dict:
    """
    Create a map {ec_number: [KO IDs]} from the parsed KO dict, using ec_numbers collected from definition_raw.
    """
    ec_to_ko = defaultdict(list)
    for ko_id, data in ko_ec.items():
        for ec in data.get('ec_numbers', []):
            ec_to_ko[ec].append(ko_id)
    # make lists unique and consistent
    return {k: sorted(set(v)) for k, v in ec_to_ko.items()}

#ko_ec = read_ko_data() # previous function called that feed this one
#ec_to_ko = build_ec_to_ko_map(ko_ec)
#ec_to_ko_sample = {k: ec_to_ko[k] for k in list(ec_to_ko.keys())[:10]}
''' {'1.1.1.1': ['K00001',
  'K11440',
  'K13951',
  'K13952',
  'K13953',
  'K13954',
  'K13980',
  'K18857'],
 '1.1.1.2': ['K00002', 'K13979'],'''


# _________________
# ### Pathway Curated   
# 
# This master assembler takes the outputs from various parsers (BRENDA, pathway_info/data, pathway_mapping) and combines them.

# In[24]:


def _std_map_id(pid: str) -> str:
    """Standardize pathway ID format - FAIL if invalid."""
    if not isinstance(pid, str):
        raise TypeError(f"Expected string, got {type(pid)}")

    m = re.match(r'^(map|ec)(\d{5})$', pid.strip(), flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Invalid pathway ID format: {pid}")

    return f"map{m.group(2)}"

def build_pathways_db_curated(record, ec_pathway_mapping, pathway_data, ko_ec, ec_to_ko):
    """
    Builds the 'pathways_db' list using a prioritized hierarchy of data sources.
    Hierarchy: 1. BRENDA -> 2. KEGG DB -> 3. KEGG KO

    Args:
        record (dict): The enzyme record being built.
        brenda_en (dict): processed data from BRENDA.
        ec_pathway_mapping (dict): map from {ec: [map_id_list]}.
        pathway_data (dict):  map from {map_id: description}.
        ko_ec (dict):  map from {ko_id: data}, used for KO lookups.
    """
    normalized_ec = normalize_ec_id(record.get('ec_number', ''))
    if not normalized_ec:
        record['pathways_db'] = []
        return record

    found_pathways = set()

    # --- Priority 2: structured KEGG EC→map mapping
    if normalized_ec in ec_pathway_mapping:
        for path_id in ec_pathway_mapping[normalized_ec]:
            std_id = path_id.replace('ec', 'map')
            name = pathway_data.get(std_id)
            if name:
                found_pathways.add(name)

    # --- Priority 3: KO as fallback
    if not found_pathways:
        ko_ids = ec_to_ko.get(normalized_ec, [])
        for ko_id in ko_ids:
            ko_rec = ko_ec.get(ko_id, {}) or {}

            # 3.a direct 'pathway' field
            ko_paths = ko_rec.get('pathway', [])
            if isinstance(ko_paths, str):
                iter_paths = [ko_paths] if ko_paths.strip() else []
            elif isinstance(ko_paths, list):
                iter_paths = ko_paths
            else:
                iter_paths = []

            for path in iter_paths:
                p = str(path).strip()
                if not p:
                    continue
                pid = p.replace('ec', 'map')  # standardize id form
                found_pathways.add(pathway_data.get(pid, pid))

            # 3.b parse definition_raw for [PATH:mapxxxxx]
            def_raw = ko_rec.get('definition_raw', '') or ''
            for pid in re.findall(r'\[PATH:(map\d{5})\]', def_raw):
                found_pathways.add(pathway_data.get(pid, pid))

    record['pathways_db'] = sorted(list(found_pathways))
    record['pathways_db'] = clean_pathway_strings(pd.Series(record['pathways_db'])).tolist()

    return record


# In[25]:


'''# test
ec = str(ECcontri_Uniprot["ec"].iloc[0])
rec = {"ec_number": ec}
rec = build_pathways_db_curated(rec, ec_pathway_mapping, pathway_data, ko_ec, ec_to_ko)
rec
{'ec_number': '1.1.1.1',
 'pathways_db': ['Biosynthesis of secondary metabolites',
  'Chloroalkane and chloroalkene degradation',
  'Drug metabolism - cytochrome P450'''


# In[26]:


pd.set_option("display.max_colwidth", True)
pd.reset_option('display.max_colwidth')


# _________________
# ### Module Database
# _________________

# In[27]:


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
#module_info = read_module_data()
'''M00001': 'Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate',
 'M00002': 'Glycolysis, core module involving three-carbon compounds',
 'M00003': 'Gluconeogenesis, oxaloacetate => fructose-6P',
 'M00004': 'Pentose phosphate pathway (Pentose p'''


# _____________________
# ### Compound Database
# _____________________

# In[28]:


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

#compound_info= read_compound_data()
#compound_info
'''{'C00001': {'name': 'H2O', 'synonyms': ['Water']},
 'C00002': {'name': 'ATP', 'synonyms': ["Adenosine 5'-triphosphate"]},
 'C00003': {'name': 'NAD+',
  'synonyms': ['NAD',
   'Nicotinamide adenine dinucleotide',
   'DPN',
   'Diphosphopyridine nucleotide',
   'Nadide',
   'beta-NAD+']},'''


# ___________
# ### Metal pdb
# ____________
# MetalPDB in 2018: a database of metal sites in biological macromolecular structures.
# Putignano V., Rosato A., Banci L., Andreini C.
# Nucleic Acids Res. 2018 Jan;46(D1):D459-D464. [PMID: 29077942]
# 

# In[29]:


def parse_metalpdb_xml():
    paths = setup_paths()
    metalpdb_path = paths['metalpdb']

    def norm_symbol(raw: str) -> str | None:
        # Only output normal periodic symbols, never species
        if not raw: return None
        r = raw.strip()
        # If mapped, just use the base symbol (no charge)
        mapped = CLEAN_METAL_MAPPING.get(r.lower()) if 'CLEAN_METAL_MAPPING' in globals() else cs.metal_mapping.get(r.lower()) #mapped = cs.metal_mapping.get(r.lower())
        if mapped:
            m = re.match(r'^([A-Z][a-z]?)\d?[+-]?$', mapped)
            if m:
                return m.group(1)
            return mapped
        # fallback: proper-case
        return r[0].upper()+r[1:].lower() if len(r)>=2 else r.upper()

    out = {}
    context = etree.iterparse(metalpdb_path, events=('end',), tag='site', recover=True)
    for _ev, site in context:
        try:
            site_name = (site.findtext('site_name') or '').strip()
            pdb_code  = (site.findtext('pdb_code')  or '').strip()

            ec_numbers = set()
            for sc in site.findall('.//site_chain') or site.findall('.//{*}site_chain'):
                ec = (sc.findtext('ec_number') or sc.findtext('{*}ec_number') or '').strip()
                if ec:
                    ec_numbers.add(ec)

            for m_idx, metal in enumerate(site.findall('.//metal') or site.findall('.//{*}metal'), start=1):
                sym = norm_symbol(metal.findtext('periodic_symbol') or metal.findtext('{*}periodic_symbol') or '')
                if not sym: 
                    continue
                cn_txt = (metal.findtext('coordination_number') or metal.findtext('{*}coordination_number') or '').strip()
                try: cn = int(cn_txt)
                except Exception: cn = None

                ligs = metal.findall('ligands/ligand') or metal.findall('.//ligand') or metal.findall('.//{*}ligand')
                residues, seen = [], set()
                for lig in ligs:
                    def gl(tag): return (lig.findtext(tag) or lig.findtext(f'{{*}}{tag}') or '').strip()
                    key = (gl('residue_name').upper(), gl('residue_pdb_number'), gl('chain_letter'), gl('endo_exo') or None)
                    if key in seen: 
                        continue
                    seen.add(key)
                    residues.append({
                        'residue_name': key[0],
                        'residue_number': key[1],
                        'chain': key[2],
                        'binding_type': key[3],
                    })

                site_key = f"{pdb_code}_{site_name}_{sym}_{m_idx:02d}"
                out[site_key] = {
                    'pdb_code': pdb_code,
                    'site_name': site_name,
                    'metal': {'symbol': sym, 'coordination_number': cn},
                    'residues': residues,
                    'ec_numbers': sorted(ec_numbers),
                }
        finally:
            site.clear()
            while site.getprevious() is not None:
                del site.getparent()[0]

    return out
metal_binding_data = parse_metalpdb_xml()
'''metal_binding_data
{'101d_101d_1_Mg_01': {'pdb_code': '101d',
  'site_name': '101d_1',
  'metal': {'symbol': 'Mg', 'coordination_number': 3},
  'residues': [{'residue_name': 'HOH',
    'residue_number': '61',
    'chain': 'A',
    'binding_type': 'exogenous'},
   {'residue_name': 'HOH',
    'residue_number': '60'''


# _____________________
# ### metals consolidated
# _____________________

# In[30]:


# Regex for detecting ionic species in EITHER format: Mg+2 OR Mg2+ (malformed BRENDA)
SPECIES_RE = re.compile(r'^[A-Z][A-Za-z0-9]*[+-]\d*$')  # Matches: Fe+2, PO4-3, MoO4-2, HCO3-, etc.

# Mapping for anions found in PDB residues (not in metal.symbol field)
ANION_RESIDUE_MAP = {
    'PO4': 'PO4-3','SO4': 'SO4-2','CO3': 'CO3-2','NO2': 'NO2-','NO3': 'NO3-','S2O3': 'S2O3-2','SO3': 'SO3-2','PO3': 'PO3-3'
}
# Elements to completely exclude from consolidated_metals
EXCLUDE_ELEMENTS = {'as', 'se', 'sr', 'pt', 'au', 'ag', 'w', 'li', 'u'}

def _normalize_ion_charge(tok: str) -> str:
    """   Convert malformed BRENDA notation to correct format.
    Mg2+ → Mg+2 Ba2+ → Ba+2 Fe3+ → Fe+3 Cd2+ → Cd+2 Al3+ → Al+3    Cr3+ → Cr+3
    Args: tok: Ion string (may be malformed)    
    Returns: Corrected ion string with Sign+Number format    """
    tok = tok.strip()
    m = re.match(r'^([A-Z][A-Za-z0-9]*?)(\d+)([+\-])$', tok)
    if m:
        return f"{m.group(1)}{m.group(3)}{m.group(2)}"
    return tok
def _proper_case_symbol(sym: str) -> str:
    """   Proper-case a symbol WITHOUT breaking compounds.
    - Simple elements (1-2 letters, no digits): fe → Fe, mg → Mg, al → Al
    - Compounds (contains digits or >2 letters): PO4 → PO4, SO4 → SO4, MoO4 → MoO4
    Args:    sym: Element or compound symbol
    Returns:  Properly cased symbol    """
    s = (sym or '').strip()
    if not s:
        return s

    # Only proper-case simple 1-2 letter elements without digits
    # Compounds stay as-is to preserve correct notation
    if len(s) <= 2 and not any(c.isdigit() for c in s):
        return s[0].upper() + s[1:].lower() if len(s) == 2 else s.upper()

    return s  # Compound or complex formula - keep as-is

def _to_symbol(val: str) -> str | None:
    """
    Extract base element/compound from any ion notation.
    Examples:Fe+2 → Fe  Mg+2 → Mg   PO4-3 → PO4  V5+ → V    Cl- → Cl
    Args: val: Ion or element string
    Returns: Base element/compound symbol, or None if invalid
    """
    if not val:
        return None
    # ADDED check for 'nan' string content:
    clean_val = val.strip()
    if clean_val.lower() in ('nan', 'none', 'na', ''):
        return None 
    # Match: base element/compound + optional charge
    # Handles: Fe, Fe+2, PO4-3, V5+, Cl-, etc.
    m = re.match(r'^([A-Z][A-Za-z0-9]*)([+-]\d*)?$', clean_val)
    # Check if the captured base element itself is 'nan'
    if m and m.group(1).lower() in ('nan', 'none', 'na'):
        return None

    return m.group(1) if m else None

def _norm_token_any(tok: str, mapping: dict) -> tuple[str | None, str | None]:
    """Normalize any metal/ion token using mapping and return (base_symbol, species).
    Process:
    1. Normalize malformed ions (Mg2+ → Mg+2)
    2. Apply mapping (words like 'magnesium' → 'Mg')
    3. Detect if result is a species (has charge)
    4. Extract base symbol

    Examples:'magnesium' → ('Mg', None) 'Mg+2' → ('Mg', 'Mg+2')
        'Fe3+' → ('Fe', 'Fe+3')  [normalized]
        'PO4-3' → ('PO4', 'PO4-3')        'zinc' → ('Zn', None)

    Args: tok: Raw token from BRENDA or PDB
          mapping: metal_mapping dictionary

    Returns: (base_symbol, species_or_None)"""
    raw = (tok or '').strip()
    if not raw:
        return (None, None)

    # Step 1: Normalize malformed ions first
    normalized = _normalize_ion_charge(raw)

    # Step 2: Apply mapping using lowercase key (handles words like 'magnesium')
    mapped = mapping.get(normalized.lower(), normalized).strip()
    mapped = _normalize_ion_charge(mapped)  # Ensure mapped value is also normalized

    # Step 3: Check if mapped result is a species (has charge notation)
    if SPECIES_RE.match(mapped):
        return (_to_symbol(mapped), mapped)

    # Step 4: Check if normalized raw token is a species
    if SPECIES_RE.match(normalized):
        return (_to_symbol(normalized), normalized)

    # Step 5: Not a species - treat as base element/word
    base = _to_symbol(mapped) or _to_symbol(normalized)
    return (_proper_case_symbol(base) if base else None, None)

# Global cache for EC → PDB symbols index
_EC2PDB: dict[str, list[str]] | None = None

def _build_ec2pdb_symbols_index(metal_binding_data: dict) -> dict[str, list[str]]:
    """
    Build EC → list of ionic species (metals + anions) from PDB data.
    Extracts: - Base metals from metal.symbol field (Mg, Zn, Fe, etc.)
              - Polyatomic anions from residues field (PO4, SO4, CO3, etc.)
    Returns lists (not sets) to preserve insertion order for deterministic results.
    Args:    metal_binding_data: Parsed MetalPDB XML structure
    Returns: Dictionary mapping EC numbers to lists of proper-cased symbols
    """
    idx: dict[str, list[str]] = {}

    for site in (metal_binding_data or {}).values():
        # Get EC numbers for this site
        ec_numbers = site.get('ec_numbers', []) or []
        if not ec_numbers:
            continue

        symbols_for_site = []

        # Extract metal symbol
        m = site.get('metal') or {}
        sym = m.get('symbol')
        if sym:
            # Normalize PDB metal symbol to proper element base (never species)
            base, _ = _norm_token_any(sym, {})  # No mapping needed for PDB
            if base and base.lower() not in EXCLUDE_ELEMENTS:
                symbols_for_site.append(base)

        # Extract anions from residues
        for residue in site.get('residues', []) or []:
            res_name = (residue.get('residue_name') or '').strip().upper()
            if res_name in ANION_RESIDUE_MAP:
                canonical_anion = ANION_RESIDUE_MAP[res_name]
                # Extract base from anion (PO4-3 → PO4)
                anion_base = _to_symbol(canonical_anion)
                if anion_base:
                    symbols_for_site.append(_proper_case_symbol(anion_base))

        # Add symbols to all EC numbers for this site
        for ec in ec_numbers:
            if not ec:
                continue
            if ec not in idx:
                idx[ec] = []
            # Add each symbol if not already present (preserve order, no duplicates)
            for symbol in symbols_for_site:
                if symbol not in idx[ec]:
                    idx[ec].append(symbol)

    return idx

def aggressive_clean(text: str) -> str:
    """Removes non-printable characters and non-breaking spaces."""
    if not isinstance(text, str):
        return ""
    # Remove non-printable control characters (ASCII 0-31) and non-breaking spaces (CHAR(160))
    cleaned = re.sub(r'[\x00-\x1F\x80-\x9F\xa0]', '', text)
    return cleaned.strip()

def build_metals_db_curated(record: dict, brenda_en: dict, allowed_symbols=None) -> dict:
    """
    Build consolidated list of ionic species (metals, cations, anions) relevant to corrosion.

    NOTE ON TERMINOLOGY: We use "metals" as shorthand for all ionic species/electrolytes affecting water chemistry
    and corrosion. This includes:
    - Metal cations (Fe+2, Cu+, Zn+2, Al+3, etc.)- Non-metal cations (H+, Na+, K+, Ca+2, Mg+2, Ba+2)
    - Anions (Cl-, F-, S-2)- Polyatomic anions (PO4-3, SO4-2, CO3-2, NO3-, NO2-, S2O3-2, SO3-2, MoO4-2)

    RATIONALE: Solid alloys in system components interact with their environment by converting to ionic species or being assimilated as organic complexes, affecting water chemistry
    and corrosion in heating/cooling systems. Cations cannot exist without anions in solution, so we treat all ionic species together as "metals" for simplicity.

    PROCESS:
    1. Extract raw metals list from BRENDA for this EC
    2. Normalize malformed BRENDA ions (Mg2+ → Mg+2, Fe3+ → Fe+3)
    3. Extract metals and anions from PDB for this EC
    4. Apply canonical charge states from DEFAULT_CANONICAL when base element detected
    5. Deduplicate by base element, preserve insertion order (BRENDA first, then PDB)
    6. Filter by allowed_symbols if provided (checks base element only)

    Args:
        record: EC record dict (must have 'ec_number')
        brenda_en: Processed BRENDA data {ec: {'brenda_metals': [...]}}
        allowed_symbols: Optional set/list to filter results (checks base element)

    Returns:
        record with added/updated fields:
        - 'brenda_metals': audit copy of raw BRENDA list
        - 'consolidated_metals': final deduplicated list of ionic species (species if available, else base)
    """
    global _EC2PDB, metal_binding_data

    # Lazy-build EC → PDB symbols index from current metal_binding_data
    if _EC2PDB is None:
        _EC2PDB = _build_ec2pdb_symbols_index(metal_binding_data or {})

    # Canonical charge states for base elements (when BRENDA only gives base element)
    DEFAULT_CANONICAL = {
        'cu': 'Cu+', 'co': 'Co+2', 'cd': 'Cd+2', 'fe': 'Fe+2', 'mn': 'Mn+2', 
        'ni': 'Ni+2', 'zn': 'Zn+2', 'pb': 'Pb+2', 'hg': 'Hg+2', 'al': 'Al+3',
        'cr': 'Cr+3', 'mg': 'Mg+2', 'ca': 'Ca+2', 'ba': 'Ba+2', 'na': 'Na+',
        'k': 'K+', 'h': 'H+',
        # Anions
        'f': 'F-', 'cl': 'Cl-', 's': 'SO4-2', 
        # Polyatomic anions
        'mo': 'MoO4-2', 'po4': 'PO4-3', 'so4': 'SO4-2', 'po3': 'PO3-3', 
        'co3': 'CO3-2', 'no2': 'NO2-', 'no3': 'NO3-', 's2o3': 'S2O3-2', 
        'so3': 'SO3-2', 'v': 'V+5', 'hco3': 'HCO3-'
    }

    metal_mapping = cs.metal_mapping  # Use corrected mapping from config

    # Build allowed base elements set (for filtering)
    allowed_base = None
    if allowed_symbols:
        allowed_base = {_proper_case_symbol(_to_symbol(s)) for s in allowed_symbols if _to_symbol(s)}
        allowed_base = {b for b in allowed_base if b and b.lower() not in EXCLUDE_ELEMENTS}
    # Get EC number and BRENDA data
    ec = normalize_ec_id(record.get('ec_number', '') or '')
    be = brenda_en.get(ec, {}) if brenda_en else {}
    br_raw = []
    for s in (be.get('brenda_metals') or []):
        # 1. Aggressively convert to string, handling non-strings and invisible chars
        try:
            # Coerce to string first
            s_str = aggressive_clean(str(s))
        except Exception:
            continue # Skip if conversion or cleaning fails unexpectedly
        is_nan_like = False
        if not s_str or s_str.lower() in ('nan', 'none', 'na'):
            is_nan_like = True

        try:
            if pd.isna(s): # Check if the original object was NaN
                is_nan_like = True
        except Exception:
            pass # Ignore errors if pd.isna can't handle 's'

        if is_nan_like:
            continue

        # 3. Append the final clean token
        br_raw.append(s_str)

    # !!! ADD THIS PRINT STATEMENT !!!
    print(f"DEBUG: br_raw contents before processing: {br_raw}")


    # Store audit copy of raw BRENDA list
    record['brenda_metals'] = list(dict.fromkeys(br_raw))

    # Get PDB symbols for this EC (already base symbols, proper-cased)
    pdb_syms = _EC2PDB.get(ec, []) or [] 
    # Track ALL unique species (allows Cu+ and Cu+2 to coexist)
    seen_species: set[str] = set()
    order: list[str] = []  # Preserves insertion order: BRENDA first, then PDB

    def observe(token: str):
        """
        Process a single token and update best species for its base element.

        Logic:
        1. Parse token to get (base_symbol, species)
        2. If species detected, store it
        3. If no species but base in DEFAULT_CANONICAL, apply canonical species
        4. Store first-seen base in order list
        """
        if not isinstance(token, str):
            return
        # Skip empty or NaN strings . new
        token = token.strip()
        if not token or token.lower() in ('nan', 'none', '', 'na'):
            return
        try:
            if pd.isna(token):
                return
        except:
            pass

        base, species = _norm_token_any(token, metal_mapping)
        if not base:
            return
        # Proper-case base consistently
        base = _proper_case_symbol(base)
        k = base.lower()
        # Exclude unwanted elements (ALL forms - base and species)
        if k in EXCLUDE_ELEMENTS:
            return
        # 🆕 ALSO check first element of compound (e.g., 'As' in 'AsO3')
        first_element_match = re.match(r'^([A-Z][a-z]?)', base)
        if first_element_match:
            first_elem = first_element_match.group(1).lower()
            if first_elem in EXCLUDE_ELEMENTS:
                print(f"  Excluding '{token}' (compound of excluded element '{first_elem}')")
                return
        # Filter by allowed base elements if provided
        if allowed_base and base not in allowed_base:
            return
        # Determine final species
        if species:
            # Already has charge - keep as-is (Fe+3 stays Fe+3, Cu+2 stays Cu+2)
            final = species
        else:
            # Base only - convert using DEFAULT_CANONICAL (Fe → Fe+2, S → SO4-2)
            final = DEFAULT_CANONICAL.get(k)
            if not final:
                return  # Skip if not in DEFAULT_CANONICAL

        # RIGHT BEFORE the final if statement:
        if final:  # Add this diagnostic
            print(f"  Token '{token}' → base='{base}', species='{species}' → final='{final}'")

        # Add if not duplicate
        final_lower = final.lower()
        if final_lower not in seen_species:
            seen_species.add(final_lower)
            order.append(final)

    # Process BRENDA first (preserves priority), then PDB
    for t in br_raw:
        observe(t)
    for t in pdb_syms:
        if t and t.lower() not in ('nan', 'none', '', 'na'): 
            observe(t)

    record['consolidated_metals'] = order

    return record


# In[31]:


'''rec = {'ec_number': '1.1.1.1'}
out = build_metals_db_curated(dict(rec), brenda_en if 'brenda_en' in globals() else brenda_data, allowed_symbols=None)
print("brenda_metals:", out.get('brenda_metals'))
print("brenda_metals_standardized:", out.get('brenda_metals_standardized'))
print("consolidated_metals:", out.get('consolidated_metals'))
print("pdb ligands for EC:", (_EC2PDB_LIGANDS.get('1.1.1.1') if '_EC2PDB_LIGANDS' in globals() else None))'''


# ⚠️ Token 'As' → base='As', species=None, NOT in DEFAULT_CANONICAL
# brenda_metals: ['Ca+2', 'Co+2', 'As', 'Fe+2', 'Mg+2', 'Cl-', 'Ni+2', 'PO4-3', 'SO4-2', 'Cr+3', 'S-2', 'Cu+2', 'F-', 'K+', 'Al+3', 'S', 'Zn+2', 'Pb+2', 'Mn+2', 'H+', 'Ba+2', 'V5+', 'Hg+2', 'Na+']
# brenda_metals_standardized: None
# consolidated_metals: ['Ca+2', 'Co+2', 'Fe+2', 'Mg+2', 'Cl-', 'Ni+2', 'PO4-3', 'SO4-2', 'Cr+3', 'S-2', 'Cu+2', 'F-', 'K+', 'Al+3', 'Zn+2', 'Pb+2', 'Mn+2', 'H+', 'Ba+2', 'V+5', 'Hg+2', 'Na+', 'Cd+2']
# pdb ligands for EC: None

# ____________________________
# ### Integrating pathways to the ECcontri_Uniprot Data 
# __________________________

# In[32]:


# Reading loading 10 min
output_path = output_large / "chunks_pathway"
chunk_files = sorted([os.path.join(output_path, f) for f in os.listdir(output_path) if f.startswith('chunk_') and f.endswith('.pkl')])
ECcontri_pathway = pd.concat((pd.read_pickle(f) for f in chunk_files), ignore_index=True)
# Extract the number from 'site_X'
ECcontri_pathway['sample_num'] = ECcontri_pathway['Sites'].str.extract('site_(\d+)').astype(int)
# Sort by the numerical part
ECcontri_pathway = ECcontri_pathway.sort_values(by='sample_num')
# Drop the temporary column
ECcontri_pathway = ECcontri_pathway.drop(columns=['sample_num'])
ECcontri_pathway = ECcontri_pathway.drop_duplicates(subset=["ec", "Sites", "Genus"])
ECcontri_pathway = ECcontri_pathway.copy()
ECcontri_pathway["Sites"] = ECcontri_pathway["Sites"].astype("category")
ECcontri_pathway["Genus"] = ECcontri_pathway["Genus"].astype("category")
ECcontri_pathway["pathway"] = ECcontri_pathway["pathway"].astype("category")
ECcontri_pathway["ipath"] = ECcontri_pathway["ipath"].astype("category")


# In[33]:


# Read the df ECcontri_Uniprot
ECcontri_Uniprot_path = output_large / 'pre_ECcontri_Uniprot.parquet'
ECcontri_Uniprot = pd.read_parquet(ECcontri_Uniprot_path)
print(f"DataFrame loaded from {ECcontri_Uniprot_path} with shape {ECcontri_Uniprot.shape}")
print(f"Memory usage after loading: {ECcontri_Uniprot.memory_usage(deep=True).sum() / 1024**2:.2f} MB") 
#32.97


# In[34]:


ECcontri_Uniprot = pd.merge(ECcontri_Uniprot, ECcontri_pathway, on=['Sites', 'Genus', 'ec'], how='left')


# In[35]:


# Making sure the nans are Nans
ECcontri_Uniprot['pathway']= ECcontri_Uniprot['pathway'].astype('str')
ECcontri_Uniprot['ipath']= ECcontri_Uniprot['ipath'].astype('str')
ECcontri_Uniprot['protein_name'] = ECcontri_Uniprot['protein_name'].astype('str')
ECcontri_Uniprot['pathway'] = ECcontri_Uniprot['pathway'].replace('nan', np.nan)
ECcontri_Uniprot['protein_name'] = ECcontri_Uniprot['protein_name'].replace('nan', np.nan)
ECcontri_Uniprot['ipath'] = ECcontri_Uniprot['ipath'].replace('nan', np.nan)
ECcontri_Uniprot["pathway"].isna().sum()


# In[36]:


del ECcontri_pathway
gc.collect()  # Run garbage collection to free up memory


# ______________
# ### Consolidated pathway column primary_pathway, creation
# _________
# In order to use the pathway to retrieve data from the db, the programatic ipath is no suitable since the db do no contain data wiht pwy id so the best option is to find them with the descriptive pathways obtained from pathway column. However when there are nans the programatic column would be used to fill the gaps when ever posible.

# In[37]:


# It takes the value from 'pathway'. If that value is null, it fills it with the value from 'ipath' from the exact same row.
ECcontri_Uniprot['pathways'] = ECcontri_Uniprot['pathway'].fillna(ECcontri_Uniprot['ipath'])
# Remove the ipath and pathway to let only one comprehensive column pathways
ECcontri_Uniprot = ECcontri_Uniprot.drop(columns=["pathway", "ipath"])
# Renaming the pathways column to pathway_primary in order to recall its origin
ECcontri_Uniprot = ECcontri_Uniprot.rename(columns={"pathways": "pathway_primary"})
# probe that the number of nans has decresed
ECcontri_Uniprot["pathway_primary"].isna().sum()


# In[38]:


ECcontri_Uniprot['pathway_primary']= ECcontri_Uniprot['pathway_primary'].astype('category')
ECcontri_Uniprot['protein_name'] = ECcontri_Uniprot['protein_name'].astype('category')
ECcontri_Uniprot['protein_name'] = ECcontri_Uniprot['protein_name'].apply(clean_protein_name)


# In[39]:


# saving the df
ECcontri_Uniprot_path = output_large / 'ECcontri_Uniprot_pathway.parquet'
ECcontri_Uniprot.to_parquet(ECcontri_Uniprot_path)


# In[40]:


# I decided that taking KO information was too memory intensive and it would take a year to process on my computer or pay an external platform to do it.
'''# Reintegrating KO information
ko_per_gensite_path = output_base / "ko_per_gensite_tabular.csv"

# reading the csv as pandas df
ko_per_gensite_df = pd.read_csv(ko_per_gensite_path)
# convert the ko column of semicolon separate string to a column of list
ko_per_gensite_df['KOs'] = ko_per_gensite_df['KOs'].apply(lambda x: x.split(';') if pd.notna(x) else [])

# back into the dictionary format {('Sites', 'Genus', 'taxon_abun'): [KOs list], ...}
ko_per_gensite = {}
for index, row in ko_per_gensite_df.iterrows():
    # Create the tuple key for the dictionary
    key_tuple = (row['Sites'], row['Genus'], row['taxon_abun'])
    # Assign the list of KOs to this key
    ko_per_gensite[key_tuple] = row['KOs']
# Add KOs via lookup
def get_ko_for_row(row):
    key = (row['Sites'], row['Genus'], row['taxon_abun'])
    return ko_per_gensite.get(key, [])

ECcontri_Uniprot['KOs'] = ECcontri_Uniprot.apply(get_ko_for_row, axis=1)
ECcontri_Uniprot = ECcontri_Uniprot.rename(columns={"KOs": "iKO"})'''


# ## 1.2. Main Function of Creating a Consolidated DataBase

# In[41]:


global brenda_data, metal_binding_data, _EC2PDB, _EC2SITES
brenda_data = parse_brenda_file() or {}
metal_binding_data = parse_metalpdb_xml() or {}
_EC2PDB = None
_EC2SITES = None
def create_metabolism_database(ECcontri_Uniprot, fc_processor, metal_processor, synergy_processor):
    """
    Build a list of dictionaries, each representing a single EC record.    
    This function builds a comprehensive database of enzyme records with detailed information
    about their potential relevance to corrosion processes. It incorporates data from various
    sources, including BRENDA, MetalPDB, KEGG, and own experimental data proccesed by Picrust2.

    Args:     ECcontri_Uniprot: Main dataframe with EC, protein, abundances, genus (1.5M rows)       
    Returns:  list: List of dictionaries, each containing information about an enzyme
    """
    global brenda_data, metal_binding_data, _EC2SITES
    try:
        # 1 Read all necessary files
        print(f"[{time.ctime()}] Function started. Processing {len(ECcontri_Uniprot)} rows.")
        print("Globals ready:",

        f"brenda_data={len(brenda_data)} entries,",
        f"metal_binding_data={len(metal_binding_data)} sites")
        print("Loading data sources...")
        ec_to_names = read_enzyme_names() or {}
        print(f"[{time.ctime()}] Loading ec_to_names...names={len(ec_to_names)}")
        enzyme_class = read_enzyme_class() or {}
        print(f"[{time.ctime()}] Loading enzyme_class...classes={len(enzyme_class)}")
        ko_ec = read_ko_data() or {}
        print(f"[{time.ctime()}] Loading ko_data...KO={len(ko_ec)}")
        ec_to_ko = build_ec_to_ko_map(ko_ec)
        print(f"[{time.ctime()}] Building ec_to_ko map...ec_to_ko ={len(ec_to_ko)}")
        hierarchy_brite = parse_ko_brite_filtered(brite_path, allowed_categories) 
        if hierarchy_brite is None or (hasattr(hierarchy_brite, 'empty') and hierarchy_brite.empty):
            hierarchy_brite =pd.DataFrame()
        print(f"[{time.ctime()}] Loading hierarchy_brite ...hierarchy_brite = {len(hierarchy_brite)}")
        pathway_data = read_pathway_data() or {}
        print(f"[{time.ctime()}] Loading pathway_data...pathways={len(pathway_data)}")
        reaction_info = read_reaction_data() or {}
        print(f"[{time.ctime()}] Loading reaction_info...rxns={len(reaction_info)}")
        ec_to_reaction = create_ec_to_reaction_mapping() or {}  
        print(f"[{time.ctime()}] Building ec_to_reaction map...ec→rxn={len(ec_to_reaction)}")
        module_info = read_module_data() or {}
        print(f"[{time.ctime()}] Loading module_info...")
        compound_info = read_compound_data() or {}
        print(f"[{time.ctime()}] Loading compound_info...")       
        brenda_en = process_brenda_data(brenda_data, processors, cs, config)  or {}
        print(f"[{time.ctime()}] Loading brenda_en...brenda={len(brenda_en)}")
        metal_binding_data = parse_metalpdb_xml() or {}
        _EC2PDB = None  # reset cache
        _EC2SITES = None 
        print(f"[{time.ctime()}] MetalPDB sites: {len(metal_binding_data)}")
        ec_pathway_mapping = read_ec_pathway_mapping() or {}
        print(f"[{time.ctime()}] Extracting ec_pathway_mapping...ec→path={len(ec_pathway_mapping)}")
        _ALLOWED = None   

        # 2. The EC that will serve to retrieve the data are selected
        unique_ecs = set()
        for ec in ECcontri_Uniprot['ec'].unique():
            normalized_ec = normalize_ec_id(ec)
            if normalized_ec:
                unique_ecs.add(normalized_ec)

        unique_ecs_list = list(unique_ecs)  # Make a list for slicing
        ec_to_names = read_enzyme_names(unique_ecs_to_filter=unique_ecs) or {}

        # 3. extract protein information
        print("Preparing protein database...")
        protein_database = []
        # From BRENDA
        if brenda_en:
            for ec_number, data in brenda_en.items():
                # NORMALIZE EC NUMBER
                normalized_ec = normalize_ec_id(ec_number)
                if not normalized_ec:
                    continue
                # 4.Use EC number to find enzyme names from ec_to_names
                if normalized_ec in ec_to_names:
                    protein_name = ec_to_names[normalized_ec][0] if ec_to_names[normalized_ec] else None

                    # 5. Apply clean_protein_name here for the protein_name being added to protein_database
                    protein_name_cleaned = clean_protein_name(protein_name) if protein_name else None
                    if protein_name_cleaned:
                        protein_database.append({
                            'source': 'BRENDA',
                            'ec_number': normalized_ec,
                            'protein_name': protein_name_cleaned,
                            'brenda_en': data
                        })  

        # 6. From MetalPDB
        if metal_binding_data:
            for site_key, site_data in metal_binding_data.items():
                for chain in site_data.get('site_chains', []):
                    if 'molecule_name' in chain and chain['molecule_type'] == 'protein':
                        protein_name = chain['molecule_name']
                        # 7.Apply clean_protein_name here for MetalPDB names
                        protein_name_cleaned = clean_protein_name(protein_name) if protein_name else None
                        if protein_name_cleaned:
                            protein_database.append({
                            'source': 'MetalPDB',
                            'pdb_code': site_data.get('pdb_code'),
                            'protein_name': protein_name_cleaned, 
                            'metal_binding': site_data.get('metal', {})
                        })

    except Exception as e:
        print(f"Error loading data sources: {e}")

        return []

    # 8. Track statistics for validation
    stats = {
        'total_enzymes': 0,
        'with_brenda_en': 0,
        'with_pathways_db': 0,
        'with_ko': 0,
        'with_metal_involvement': 0,
        'with_corrosion_mechanisms': 0,
        'with_reaction_db': 0
    }

    start_time_record = time.time()
    # 9. Prepare a list to store all records
    ec_records = []
    stats['total_enzymes'] = len(unique_ecs)

    print("Pre-computing enzyme class lookups...")
    ec_class_lookup = {}
    for normalized_ec in unique_ecs:
        try:
            ec_prefix = '.'.join(normalized_ec.split('.')[:2])
            # 10.Try exact match first
            if ec_prefix in enzyme_class:
                ec_class_lookup[normalized_ec] = enzyme_class[ec_prefix].strip()
            # 11.Then try pattern match
            else:
                pattern_key = f"{ec_prefix}.-.-"
                if pattern_key in enzyme_class:
                    ec_class_lookup[normalized_ec] = enzyme_class[pattern_key]

        except Exception:
            pass  # Skip if error occurs

    # 12. Create lookup dictionary for faster matching
    print("Building protein lookup..")
    protein_lookup = {}
    for db_rec in protein_database:
        protein_name = db_rec.get('protein_name')
        protein_name_cleaned = clean_protein_name(protein_name) if isinstance(protein_name, str) else None
        protein_name_key = enhanced_clean_protein_name(protein_name_cleaned) if isinstance(protein_name_cleaned, str) else None
        if protein_name_key:
            protein_lookup[protein_name_key] = {**db_rec, 'protein_name': protein_name_cleaned}

    # 15.Populate from ec_to_names for basic enzyme names
    print("Building enzyme records...")

    chunk_size = 500 # EC records per batch - balances memory vs API calls
    for chunk_start in range(0, len(unique_ecs_list), chunk_size):
        chunk_ecs = unique_ecs_list[chunk_start:chunk_start + chunk_size]
        print(f"Processing chunk {chunk_start//chunk_size + 1}")
        for ec_number in chunk_ecs:
            # normalise ec number
            normalized_ec = normalize_ec_id(ec_number)
            if not normalized_ec:
                print(f"Skipping invalid EC number: {ec_number}")
                continue

            names = ec_to_names.get(normalized_ec, [])

            if names:  # checks if the list is non-empty
                first_name = names[0]
            else:
                first_name = ''
            # 16. clean protein name
            #protein_name_cleaned = clean_protein_name(first_name) if isinstance(first_name, str) else None
            names_list = names if isinstance(names, list) else [names]
            names_clean = [clean_protein_name(str(n).strip()) for n in names_list if n]
            # Choose protein_name from first valid name (or empty string)
            protein_name_cleaned = names_clean[0] if names_clean else ''

            # 17. create a record for each EC number
            record ={
                'ec_number': normalized_ec,
                'enzyme_names':  names_clean,
                'protein_name': protein_name_cleaned if protein_name_cleaned else '',   
                'enzyme_class': ec_class_lookup.get(normalized_ec, ''),
                'pathways_db': [],
                'reaction_db':[],
                'metal_binding_info': {},  # from MetalPDB
                'hierarchy': [],
                'ko': [],
                'reaction': [],                       
                'consolidated_metals': [],                 
                'corrosion_mechanisms': [],
                'operational_environmental_factors': [],
                'residue_name': [],
                'coordination_number': []
            }
            # 21. Add pathways_db the curated pathway resulting from the consolidation of the db info
            record = build_pathways_db_curated(record, ec_pathway_mapping, pathway_data, ko_ec, ec_to_ko)
            if record['pathways_db']:
                stats['with_pathways_db'] += 1

            # 22 Add KO IDs
            record['ko'] = ec_to_ko.get(normalized_ec, [])
            if record['ko']:
                stats['with_ko'] += 1

            # 23. Add module information
            #for module_id, module_desc in module_info.items():
            #   if f"[EC:{normalized_ec}]" in module_desc:
            #      record['modules'].append({'id': module_id, 'description': module_desc})

            # 24. Add reaction information from brenda and from reaction_info dataframe that we joined back on reaction mapping
            record = build_reactions_db_curated(record=record, brenda_en=brenda_en,
                ec_to_reaction=ec_to_reaction, reaction_info=reaction_info)
            if record['reaction_db']:
                stats['with_reaction_db'] += 1

            # 25. Add BRENDA metal information
            record = build_metals_db_curated(record, brenda_en, allowed_symbols=None)

            # 26 Directly extend the record with the mechanisms already processed in `process_brenda_data`.
            precalculated_mechanisms = brenda_en.get(normalized_ec,{}).get('corrosion_mechanisms', [])
            if precalculated_mechanisms:
                record['corrosion_mechanisms'].extend(precalculated_mechanisms)
            stats['with_brenda_en'] += 1  
            # Copy text fields for richer scoring if present
            record['operational_environmental_factors'] = brenda_en.get(normalized_ec,{}).get('operational_environmental_factors', [])

            ec_records.append(record)

    elapsed_time_record = time.time() - start_time_record
    print(f"Processing took {elapsed_time_record:.2f} seconds")

     # 27. Process records in batches to avoid memory issues
    print("Processing protein names and calculating scores...")
    batch_size = 1000  # Protein records per batch - optimized for 16GB RAM
    for i in range(0, len(ec_records), batch_size):
        batch = ec_records[i:i+batch_size]

        # 28. Process protein names in batch
        for rec_item in batch:
            # Process protein name matches
            try:
                protein_name_norm = enhanced_clean_protein_name(clean_protein_name(rec_item.get('protein_name', ''))) 

                if len(protein_name_norm) > 2:
                    # 29 Direct match first
                    if protein_name_norm in protein_lookup:
                        db_rec_match = protein_lookup[protein_name_norm] # Assign matched record to a variable
                        # 31 Process MetalPDB data, copy metalpdb info if present        
                        if 'metal_binding' in db_rec_match: 
                            rec_item.setdefault('metal_binding_info', {})
                            metal = db_rec_match['metal_binding'].get('symbol')
                            if metal:
                                rec_item['metal_binding_info'][metal] = db_rec_match['metal_binding']

                    else:  # Fuzzy match on the same normalized keys
                        for db_name_key, db_rec_val in protein_lookup.items():
                            if db_name_key in protein_name_norm or protein_name_norm in db_name_key:
                                # copy BRENDA metals if present
                                if 'brenda_en' in db_rec_val:
                                    seed = db_rec_val['brenda_en'].get('brenda_metals', [])

                                # copy MetalPDB info if present
                                if 'metal_binding' in db_rec_val:
                                    rec_item.setdefault('metal_binding_info', {})
                                    m = db_rec_val['metal_binding'].get('symbol')
                                    if m:
                                        rec_item['metal_binding_info'][m] = db_rec_val['metal_binding']
                                        break

            except Exception as e:
                print(f"Error processing protein name {rec_item.get('protein_name')}: {e}") # Use rec_item

        # 32. Calculate scores for this batch
        for rec in batch:
            # Build text for scoring once
            try:
                enzyme_names = rec.get('enzyme_names', []) or []
                class_text = rec.get('enzyme_class', '') or ''
                pathways_txt = ' '.join(rec.get('pathways_db', []) or [])
                rxn_text     = ' '.join(rec.get('reaction_db', []) or [])

                # 33. Combined text for all scoring
                subs = []; inhib = []; cofac = []; env = []
                ec_local = rec.get('ec_number')
                if brenda_en and ec_local in brenda_en:
                    be = brenda_en.get(ec_local, {})
                    subs  = be.get('substrates', []) or []
                    inhib = be.get('inhibitors', []) or []
                    cofac = be.get('cofactors', []) or []
                    comp  = be.get('compounds')  or []
                    env   = be.get('operational_environmental_factors', []) or []

                text_parts = [
                    ' '.join(enzyme_names),
                    class_text,
                    pathways_txt,
                    rxn_text,
                    rec.get('protein_name', '') or '',
                    ' '.join(subs), ' '.join(inhib), ' '.join(cofac), ' '.join(env)
                ]
                text = ' '.join([t for t in text_parts if t]).lower()

                # 34. Use the scoring system module
                score_results = cs.calculate_overall_scores(text,  processors={'fc_processor': fc_processor,'metal_processor': metal_processor,
                                                                'synergy_processor': synergy_processor  },
                                 config=config, brenda_metals=rec.get('consolidated_metals', []))

                # 35. Merge the score results into the record.
                # Preserve mechanisms from BRENDA (score_results may contain an empty list)
                existing_mechanisms = list(rec.get('corrosion_mechanisms', []))
                score_results.pop('corrosion_mechanisms', None)  
                rec.update(score_results)
                if existing_mechanisms:
                    rec['corrosion_mechanisms'] = list({*existing_mechanisms, *rec.get('corrosion_mechanisms', [])})

            except Exception as e:
                print(f"Error scoring data for {rec.get('ec_number')}: {e}")

    # 40. Process metal binding in a separate pass one-time index: EC -> list of PDB sites
    def _build_ec2sites_index(_mbd: dict):
        idx = defaultdict(list)
        for site in (_mbd or {}).values():
            ecs = site.get('ec_numbers', []) or site.get('ecs', [])
            if not ecs:
                continue
            for ec in ecs:
                if ec:
                    idx[ec].append(site)
        return dict(idx)

    print("Processing metal binding information from MetalPDB ...")
    try:
        if _EC2SITES is None:
            _EC2SITES = _build_ec2sites_index(metal_binding_data)

        for rec in ec_records:
            ec = normalize_ec_id(rec.get('ec_number') or '')
            if not ec:
                continue

            # accumulate per-metal residue counts and coordination numbers
            res_counts_by_metal = defaultdict(Counter)  # metal -> Counter(AA)
            cn_list_by_metal = defaultdict(list)        # metal -> [cn,...]

            for site in _EC2SITES.get(ec, []):
                m = site.get('metal', {}) or {}
                sym = m.get('symbol')
                if not sym:
                    continue
                cn = m.get('coordination_number')
                if cn is not None:
                    cn_list_by_metal[sym].append(cn)

                # collect residue 3-letter codes only
                for r in site.get('residues', []) or []:
                    aa = (r.get('residue_name') or '').strip().upper()
                    if aa:
                        res_counts_by_metal[sym][aa] += 1

            # build per-metal summary: top residues + modal coordination number
            metal_binding_info = {}
            overall_residues = Counter()
            overall_cns = []

            for sym in sorted(set(list(res_counts_by_metal.keys()) + list(cn_list_by_metal.keys()))):
                # top 3 residues by frequency for this metal
                top_res = [aa for aa, _cnt in res_counts_by_metal[sym].most_common(3)]
                # modal coordination number for this metal (break ties by larger count, then higher cn)
                cns = cn_list_by_metal[sym]
                if cns:
                    cn_mode = Counter(cns).most_common(1)[0][0]
                else:
                    cn_mode = None

                metal_binding_info[sym] = {
                    'residue_names': top_res,
                    'coordination_number': cn_mode 
                }

                overall_residues.update(res_counts_by_metal[sym])
                if cn_mode is not None:
                    overall_cns.append(cn_mode)

            # attach to record
            rec['metal_binding_info'] = metal_binding_info
            # overall per-EC summaries (we asked for these fields explicitly)
            rec['residue_name'] = [aa for aa, _ in overall_residues.most_common(5)]  # top 5 overall AAs
            rec['coordination_number'] = sorted(Counter(overall_cns).keys())         # unique modes across metals

    except Exception as e:
        print(f"Error processing metal binding data (PDB aggregation): {e}")

    # 45 Process KO hierarchy information
    print("Processing KO hierarchy...")
    # hierarchy_brite: DataFrame with columns ['ec_number', 'hierarchy']

    # ec_records is a list of dicts
    for rec in ec_records:
        ec = rec.get('ec_number')
        # Find all hierarchy entries for that ec in hierarchy_brite
        hierarchies = []
        if not hierarchy_brite.empty:
            mask = hierarchy_brite['ec_number'].astype(str).str.contains(# rf'(^|[;\s]){re.escape(ec)}($|[;\s])'
                 rf'(?:^|[;\s]){re.escape(ec)}(?:$|[;\s])', regex=True, na=False
            )

            hierarchies = hierarchy_brite.loc[mask, 'hierarchy'].unique()
        # Ensure 'hierarchy' field is present and a list
        rec.setdefault('hierarchy', [])
        # Add each unique hierarchy if not already present
        for hier in hierarchies:
            if hier and hier not in rec['hierarchy']:
                rec['hierarchy'].append(hier)

    # 49. Calculate final corrosion relevance scores using the scoring module
    print("Calculating corrosion relevance scores...")

    for rec in ec_records:
        try:
            # First, ensure rec is a dictionary
            if not isinstance(rec, dict):
                print(f"Warning: Found non-dict record: {type(rec)} - {rec}")
                continue

            # 50. Assign corrosion mechanisms based on all available pathway data
            existing_mechanisms = set(rec.get('corrosion_mechanisms', []))
            #existing_mechanisms.update(pathway_mechanisms)
            rec['corrosion_mechanisms'] = list(existing_mechanisms)

            # 51. Calculate final relevance score using the correct overall scores
            if 'corrosion_relevance_score' not in rec:
                rec['corrosion_relevance_score'] = (
                    float(rec.get('overall_functional_score', 0.0)) +
                    float(rec.get('overall_metal_score', 0.0)) +
                    float(rec.get('overall_synergy_score', 0.0))
                )           

        except Exception as e:
                print(f"Error calculating corrosion score for {rec.get('ec_number', 'unknown')}: {e}")
                rec['corrosion_relevance_score'] = rec.get('corrosion_relevance_score', 0.0)
                continue

    ko_counts = [len(rec.get('ko', [])) for rec in ec_records]
    print("Records with ≥1 KO:", sum(1 for c in ko_counts if c > 0))
        #===============================================================================================

    # 53. Filter records without content (rest of the code remains the same)
    print("Filtering records...")
    filtered_ec_records = []
    for rec in ec_records:
        # Skip non-dict records
        if not isinstance(rec, dict):
            print(f"Skipping non-dict record during filtering: {type(rec)}")
            continue
        val = rec.get('protein_name')
        protein_name_cleaned = enhanced_clean_protein_name(val) if isinstance(val, str) else None
        protein_name_cleaned = protein_name_cleaned or ""

        enzyme_names = rec.get('enzyme_names', [])
        enzyme_names = [] if enzyme_names is None else enzyme_names
        normalized_ec = rec.get('ec_number', "")

        # 54.Condition 1: At least one valid identifier must be present
        has_valid_protein = bool(protein_name_cleaned) and len(protein_name_cleaned) > 2 and ("uncharacterized" not in protein_name_cleaned)
        has_valid_enzyme = any(len(name) > 2 for name in enzyme_names)
        has_valid_ec = normalized_ec.count('.') == 3 and all(part.isdigit() for part in normalized_ec.split('.'))

        # 55.Condition 2: Check for valuable data that should be preserved
        has_mechanisms = len(rec.get('corrosion_mechanisms', [])) > 0
        has_pathways_db = len(rec.get('pathways_db', [])) > 0
        has_metals_consolidated = len(rec.get('consolidated_metals', [])) > 0

        # 56.Include record if it meets either condition
        if (has_valid_protein or has_valid_enzyme or has_valid_ec) or (has_mechanisms or has_pathways_db or has_metals_consolidated):
            filtered_ec_records.append(rec)

    # 57. Replace original list with filtered version
    ec_records = filtered_ec_records
    stats['with_corrosion_mechanisms'] = sum(1 for rec in ec_records if rec.get('corrosion_mechanisms'))
    stats['with_pathways_db'] = sum(1 for rec in ec_records if rec.get('pathways_db')) #the instead silenced

    # 58. Print summary statistics
    print("\nMetabolism Database Summary:")
    print(f"Total enzyme records: {stats['total_enzymes']}")

    # A variable is created for the total to avoid repetition and add a safety check.
    total_records = len(ec_records)

    if total_records > 0:
        n_with_brenda = sum(1 for rec in ec_records if rec.get('ec_number') in brenda_en)
        print(f"Records with BRENDA data: {n_with_brenda} ({(n_with_brenda/total_records*100):.1f}%)")
        print(f"Records with reaction_db: {stats.get('with_reaction_db', 0)}")
        print(f"Total DB Pathways: {sum(len(rec.get('pathways_db', [])) for rec in ec_records)}")
        print(f"Total KO terms: {sum(len(rec.get('ko', [])) for rec in ec_records)}")
        print(f"Total corrosion mechanisms: {sum(len(rec.get('corrosion_mechanisms', [])) for rec in ec_records)}")
        print(f"Total consolidated metals: {sum(len(rec.get('consolidated_metals', [])) for rec in ec_records)}")
        print(f"Total corrosion synergies: {sum(len(rec.get('corrosion_synergies', [])) for rec in ec_records)}")
        print(f"Total functional categories: {sum(len(rec.get('functional_categories', [])) for rec in ec_records)}")
        print(f"Total operational environmental factors: {sum(len(rec.get('operational_environmental_factors', [])) for rec in ec_records)}")
    else:
        print("No enzyme records were processed to generate statistics.")

    print(f"[{time.ctime()}] finishing ec_records...")
    notify_complete()
    return ec_records


# In[ ]:





# In[42]:


# Regex for detecting ionic species in EITHER format: Mg+2 OR Mg2+ (malformed BRENDA)
SPECIES_RE = re.compile(r'^[A-Z][A-Za-z0-9]*[+-]\d*$')  # Matches: Fe+2, PO4-3, MoO4-2, HCO3-, etc.

# Mapping for anions found in PDB residues (not in metal.symbol field)
ANION_RESIDUE_MAP = {
    'PO4': 'PO4-3','SO4': 'SO4-2','CO3': 'CO3-2','NO2': 'NO2-','NO3': 'NO3-','S2O3': 'S2O3-2','SO3': 'SO3-2','PO3': 'PO3-3'
}
# Elements to completely exclude from consolidated_metals
EXCLUDE_ELEMENTS = {'as', 'se', 'sr', 'pt', 'au', 'ag', 'w', 'li', 'u'}

def _normalize_ion_charge(tok: str) -> str:
    """   Convert malformed BRENDA notation to correct format.
    Mg2+ → Mg+2 Ba2+ → Ba+2 Fe3+ → Fe+3 Cd2+ → Cd+2 Al3+ → Al+3    Cr3+ → Cr+3
    Args: tok: Ion string (may be malformed)    
    Returns: Corrected ion string with Sign+Number format    """
    tok = tok.strip()
    m = re.match(r'^([A-Z][A-Za-z0-9]*?)(\d+)([+\-])$', tok)
    if m:
        return f"{m.group(1)}{m.group(3)}{m.group(2)}"
    return tok
def _proper_case_symbol(sym: str) -> str:
    """   Proper-case a symbol WITHOUT breaking compounds.
    - Simple elements (1-2 letters, no digits): fe → Fe, mg → Mg, al → Al
    - Compounds (contains digits or >2 letters): PO4 → PO4, SO4 → SO4, MoO4 → MoO4
    Args:    sym: Element or compound symbol
    Returns:  Properly cased symbol    """
    s = (sym or '').strip()
    if not s:
        return s

    # Only proper-case simple 1-2 letter elements without digits
    # Compounds stay as-is to preserve correct notation
    if len(s) <= 2 and not any(c.isdigit() for c in s):
        return s[0].upper() + s[1:].lower() if len(s) == 2 else s.upper()

    return s  # Compound or complex formula - keep as-is

def _to_symbol(val: str) -> str | None:
    """
    Extract base element/compound from any ion notation.
    Examples:Fe+2 → Fe  Mg+2 → Mg   PO4-3 → PO4  V5+ → V    Cl- → Cl
    Args: val: Ion or element string
    Returns: Base element/compound symbol, or None if invalid
    """
    if not val:
        return None
    # ADDED check for 'nan' string content:
    clean_val = val.strip()
    if clean_val.lower() in ('nan', 'none', 'na', ''):
        return None 
    # Match: base element/compound + optional charge
    # Handles: Fe, Fe+2, PO4-3, V5+, Cl-, etc.
    m = re.match(r'^([A-Z][A-Za-z0-9]*)([+-]\d*)?$', clean_val)
    # Check if the captured base element itself is 'nan'
    if m and m.group(1).lower() in ('nan', 'none', 'na'):
        return None

    return m.group(1) if m else None

def _norm_token_any(tok: str, mapping: dict) -> tuple[str | None, str | None]:
    """Normalize any metal/ion token using mapping and return (base_symbol, species).
    Process:
    1. Normalize malformed ions (Mg2+ → Mg+2)
    2. Apply mapping (words like 'magnesium' → 'Mg')
    3. Detect if result is a species (has charge)
    4. Extract base symbol

    Examples:'magnesium' → ('Mg', None) 'Mg+2' → ('Mg', 'Mg+2')
        'Fe3+' → ('Fe', 'Fe+3')  [normalized]
        'PO4-3' → ('PO4', 'PO4-3')        'zinc' → ('Zn', None)

    Args: tok: Raw token from BRENDA or PDB
          mapping: metal_mapping dictionary

    Returns: (base_symbol, species_or_None)"""
    raw = (tok or '').strip()
    if not raw:
        return (None, None)

    # Step 1: Normalize malformed ions first
    normalized = _normalize_ion_charge(raw)

    # Step 2: Apply mapping using lowercase key (handles words like 'magnesium')
    mapped = mapping.get(normalized.lower(), normalized).strip()
    mapped = _normalize_ion_charge(mapped)  # Ensure mapped value is also normalized

    # Step 3: Check if mapped result is a species (has charge notation)
    if SPECIES_RE.match(mapped):
        return (_to_symbol(mapped), mapped)

    # Step 4: Check if normalized raw token is a species
    if SPECIES_RE.match(normalized):
        return (_to_symbol(normalized), normalized)

    # Step 5: Not a species - treat as base element/word
    base = _to_symbol(mapped) or _to_symbol(normalized)
    return (_proper_case_symbol(base) if base else None, None)

# Global cache for EC → PDB symbols index
_EC2PDB: dict[str, list[str]] | None = None

def _build_ec2pdb_symbols_index(metal_binding_data: dict) -> dict[str, list[str]]:
    """
    Build EC → list of ionic species (metals + anions) from PDB data.
    Extracts: - Base metals from metal.symbol field (Mg, Zn, Fe, etc.)
              - Polyatomic anions from residues field (PO4, SO4, CO3, etc.)
    Returns lists (not sets) to preserve insertion order for deterministic results.
    Args:    metal_binding_data: Parsed MetalPDB XML structure
    Returns: Dictionary mapping EC numbers to lists of proper-cased symbols
    """
    idx: dict[str, list[str]] = {}

    for site in (metal_binding_data or {}).values():
        # Get EC numbers for this site
        ec_numbers = site.get('ec_numbers', []) or []
        if not ec_numbers:
            continue

        symbols_for_site = []

        # Extract metal symbol
        m = site.get('metal') or {}
        sym = m.get('symbol')
        if sym:
            # Normalize PDB metal symbol to proper element base (never species)
            base, _ = _norm_token_any(sym, {})  # No mapping needed for PDB
            if base and base.lower() not in EXCLUDE_ELEMENTS:
                symbols_for_site.append(base)

        # Extract anions from residues
        for residue in site.get('residues', []) or []:
            res_name = (residue.get('residue_name') or '').strip().upper()
            if res_name in ANION_RESIDUE_MAP:
                canonical_anion = ANION_RESIDUE_MAP[res_name]
                # Extract base from anion (PO4-3 → PO4)
                anion_base = _to_symbol(canonical_anion)
                if anion_base:
                    symbols_for_site.append(_proper_case_symbol(anion_base))

        # Add symbols to all EC numbers for this site
        for ec in ec_numbers:
            if not ec:
                continue
            if ec not in idx:
                idx[ec] = []
            # Add each symbol if not already present (preserve order, no duplicates)
            for symbol in symbols_for_site:
                if symbol not in idx[ec]:
                    idx[ec].append(symbol)

    return idx

def aggressive_clean(text: str) -> str:
    """Removes non-printable characters and non-breaking spaces."""
    if not isinstance(text, str):
        return ""
    # Remove non-printable control characters (ASCII 0-31) and non-breaking spaces (CHAR(160))
    cleaned = re.sub(r'[\x00-\x1F\x80-\x9F\xa0]', '', text)
    return cleaned.strip()

def build_metals_db_curated(record: dict, brenda_en: dict, allowed_symbols=None) -> dict:
    """
    Build consolidated list of ionic species (metals, cations, anions) relevant to corrosion.

    NOTE ON TERMINOLOGY: We use "metals" as shorthand for all ionic species/electrolytes affecting water chemistry
    and corrosion. This includes:
    - Metal cations (Fe+2, Cu+, Zn+2, Al+3, etc.)- Non-metal cations (H+, Na+, K+, Ca+2, Mg+2, Ba+2)
    - Anions (Cl-, F-, S-2)- Polyatomic anions (PO4-3, SO4-2, CO3-2, NO3-, NO2-, S2O3-2, SO3-2, MoO4-2)

    RATIONALE: Solid alloys in system components interact with their environment by converting to ionic species or being assimilated as organic complexes, affecting water chemistry
    and corrosion in heating/cooling systems. Cations cannot exist without anions in solution, so we treat all ionic species together as "metals" for simplicity.

    PROCESS:
    1. Extract raw metals list from BRENDA for this EC
    2. Normalize malformed BRENDA ions (Mg2+ → Mg+2, Fe3+ → Fe+3)
    3. Extract metals and anions from PDB for this EC
    4. Apply canonical charge states from DEFAULT_CANONICAL when base element detected
    5. Deduplicate by base element, preserve insertion order (BRENDA first, then PDB)
    6. Filter by allowed_symbols if provided (checks base element only)

    Args:
        record: EC record dict (must have 'ec_number')
        brenda_en: Processed BRENDA data {ec: {'brenda_metals': [...]}}
        allowed_symbols: Optional set/list to filter results (checks base element)

    Returns:
        record with added/updated fields:
        - 'brenda_metals': audit copy of raw BRENDA list
        - 'consolidated_metals': final deduplicated list of ionic species (species if available, else base)
    """
    global _EC2PDB, metal_binding_data

    # Lazy-build EC → PDB symbols index from current metal_binding_data
    if _EC2PDB is None:
        _EC2PDB = _build_ec2pdb_symbols_index(metal_binding_data or {})

    # Canonical charge states for base elements (when BRENDA only gives base element)
    DEFAULT_CANONICAL = {
        'cu': 'Cu+', 'co': 'Co+2', 'cd': 'Cd+2', 'fe': 'Fe+2', 'mn': 'Mn+2', 
        'ni': 'Ni+2', 'zn': 'Zn+2', 'pb': 'Pb+2', 'hg': 'Hg+2', 'al': 'Al+3',
        'cr': 'Cr+3', 'mg': 'Mg+2', 'ca': 'Ca+2', 'ba': 'Ba+2', 'na': 'Na+',
        'k': 'K+', 'h': 'H+',
        # Anions
        'f': 'F-', 'cl': 'Cl-', 's': 'SO4-2', 
        # Polyatomic anions
        'mo': 'MoO4-2', 'po4': 'PO4-3', 'so4': 'SO4-2', 'po3': 'PO3-3', 
        'co3': 'CO3-2', 'no2': 'NO2-', 'no3': 'NO3-', 's2o3': 'S2O3-2', 
        'so3': 'SO3-2', 'v': 'V+5', 'hco3': 'HCO3-'
    }

    metal_mapping = cs.metal_mapping  # Use corrected mapping from config

    # Build allowed base elements set (for filtering)
    allowed_base = None
    if allowed_symbols:
        allowed_base = {_proper_case_symbol(_to_symbol(s)) for s in allowed_symbols if _to_symbol(s)}
        allowed_base = {b for b in allowed_base if b and b.lower() not in EXCLUDE_ELEMENTS}
    # Get EC number and BRENDA data
    ec = normalize_ec_id(record.get('ec_number', '') or '')
    be = brenda_en.get(ec, {}) if brenda_en else {}
    br_raw = []
    for s in (be.get('brenda_metals') or []):
        # 1. Aggressively convert to string, handling non-strings and invisible chars
        try:
            # Coerce to string first
            s_str = aggressive_clean(str(s))
        except Exception:
            continue # Skip if conversion or cleaning fails unexpectedly
        is_nan_like = False
        if not s_str or s_str.lower() in ('nan', 'none', 'na'):
            is_nan_like = True

        try:
            if pd.isna(s): # Check if the original object was NaN
                is_nan_like = True
        except Exception:
            pass # Ignore errors if pd.isna can't handle 's'

        if is_nan_like:
            continue

        # 3. Append the final clean token
        br_raw.append(s_str)

    # !!! ADD THIS PRINT STATEMENT !!!
    print(f"DEBUG: br_raw contents before processing: {br_raw}")


    # Store audit copy of raw BRENDA list
    record['brenda_metals'] = list(dict.fromkeys(br_raw))

    # Get PDB symbols for this EC (already base symbols, proper-cased)
    pdb_syms = _EC2PDB.get(ec, []) or [] 
    # Track ALL unique species (allows Cu+ and Cu+2 to coexist)
    seen_species: set[str] = set()
    order: list[str] = []  # Preserves insertion order: BRENDA first, then PDB

    def observe(token: str):
        """
        Process a single token and update best species for its base element.

        Logic:
        1. Parse token to get (base_symbol, species)
        2. If species detected, store it
        3. If no species but base in DEFAULT_CANONICAL, apply canonical species
        4. Store first-seen base in order list
        """
        if not isinstance(token, str):
            return
        # Skip empty or NaN strings . new
        token = token.strip()
        if not token or token.lower() in ('nan', 'none', '', 'na'):
            return
        try:
            if pd.isna(token):
                return
        except:
            pass

        base, species = _norm_token_any(token, metal_mapping)
        if not base:
            return
        # Proper-case base consistently
        base = _proper_case_symbol(base)
        k = base.lower()
        # Exclude unwanted elements (ALL forms - base and species)
        if k in EXCLUDE_ELEMENTS:
            return
        # 🆕 ALSO check first element of compound (e.g., 'As' in 'AsO3')
        first_element_match = re.match(r'^([A-Z][a-z]?)', base)
        if first_element_match:
            first_elem = first_element_match.group(1).lower()
            if first_elem in EXCLUDE_ELEMENTS:
                print(f"  Excluding '{token}' (compound of excluded element '{first_elem}')")
                return
        # Filter by allowed base elements if provided
        if allowed_base and base not in allowed_base:
            return
        # Determine final species
        if species:
            # Already has charge - keep as-is (Fe+3 stays Fe+3, Cu+2 stays Cu+2)
            final = species
        else:
            # Base only - convert using DEFAULT_CANONICAL (Fe → Fe+2, S → SO4-2)
            final = DEFAULT_CANONICAL.get(k)
            if not final:
                return  # Skip if not in DEFAULT_CANONICAL

        # RIGHT BEFORE the final if statement:
        if final:  # Add this diagnostic
            print(f"  Token '{token}' → base='{base}', species='{species}' → final='{final}'")

        # Add if not duplicate
        final_lower = final.lower()
        if final_lower not in seen_species:
            seen_species.add(final_lower)
            order.append(final)

    # Process BRENDA first (preserves priority), then PDB
    for t in br_raw:
        observe(t)
    for t in pdb_syms:
        if t and t.lower() not in ('nan', 'none', '', 'na'): 
            observe(t)

    record['consolidated_metals'] = order

    return record


# In[43]:


## Sample Uniprot to test
ECcontri_Uniprot= ECcontri_Uniprot.sort_values("idx", ascending=True)
ECcontri_Uniprot_sample = ECcontri_Uniprot.sample(n=15000, random_state=42)
ECcontri_Uniprot_sample= ECcontri_Uniprot_sample.sort_values("idx", ascending=True)


# output EC_RECORDS

# In[44]:


# Takes 30 - 240 min
ec_records = create_metabolism_database(ECcontri_Uniprot, fc_processor, metal_processor, synergy_processor)


# In[45]:


ec_records_df = pd.DataFrame(ec_records)


# In[46]:


sorted_metals = sorted(ec_records_df["consolidated_metals"].explode().astype(str).unique())
print(sorted_metals)


# In[47]:


# Count occurrences of each value in both columns
enzyme_counts = ec_records_df['enzyme_names'].explode().value_counts()
protein_counts = ec_records_df['protein_name'].value_counts()

# Get the number of unique values in each column
enzyme_unique_count = ec_records_df['enzyme_names'].explode().nunique()
protein_unique_count = ec_records_df['protein_name'].nunique()

# Print the results
print("Enzyme Names - Unique Values:", enzyme_unique_count)
print("Enzyme Names - Top 5 Occurrences:\n", enzyme_counts.head())
print("\nProtein Names - Unique Values:", protein_unique_count)
print("Protein Names - Top 5 Occurrences:\n", protein_counts.head())


# Enzyme Names - Unique Values: 7791
# Enzyme Names - Top 5 Occurrences:
#  enzyme_names
# 23s rrna-methyltransferase    13
# trna-methyltransferase        12
# 16s rrna-methyltransferase     9
# malate-dehydrogenase           8
# enoyl-reductase                7
# Name: count, dtype: int64
# 
# Protein Names - Unique Values: 1301
# Protein Names - Top 5 Occurrences:
#  protein_name
#                               63
# 23s rrna-methyltransferase    13
# 16s rrna-methyltransferase     9
# trna-methyltransferase         8
# malate-dehydrogenase           4
# Name: count, dtype: int64

# ### Checking the output
# 'ec_number', 'protein_name',
# 'enzyme_names', # enzymes have list muss be integrated as fallback for protein_name
# __Decisive Scoring Columns__
# 'functional_matches_detailed', # to keep procedence and alignment
# 'functional_categories_detected', # subcategories from fc
# 'functional_child_terms', # child terms fc
# 'brenda_metals', # better than consolidated
# 'synergy_child_terms_found', # child terms found in synergies
# 'synergy_categories_involved', # subcategories involved in synergies
# 'synergy_description', # origin type of synergy
# __Other Columns__
#  'enzyme_class','pathways_db', 'reaction_db', 'hierarchy', # from combined sources 
# 'corrosion_mechanism_categories', # subcategories mechanisms     
# 'corrosion_mechanisms', # child categories
# 'detected_pathways', # subcategories pathways
# ? 'pathway_matches_detailed', # subcategories + child terms
#  'operational_environmental_factors', # child categories
# ? 'operational_matches_detailed' # subcategories + child terms
# __Scoring Columns__
# 'overall_functional_score', # fc score overall
# 'overall_metal_score',  # metal score overall
# 'overall_synergy_score' # overall synergy score
# 'corrosion_relevance_score', # relevance overall
# ______________________
# **I am no seeing the following columns utility**
# 'synergy_type', # useless --type of the synergy
# 'reaction', # wrong extracted
# 'consolidated_metals', # wrong extracted
# 'detected_metal_categories', # no info
# 'functional_categories', # fc strings + scores
# 'metal_categories_detected', # metal strings + scores , useless
# 'metal_matches_detailed', # duplicated info 
# 'functional_score', # scores fc
# fc_cooccurrence_synergy_hit', # useless
# 'metal_score', # scores metal
# 'synergy_score', # synergy score overall
# 'metal_count', # number of metals
# 'legacy_synergy_groups', # useless
# 'metal_binding_info',  'ko','coordination_patterns', 'residue_binding', # maybe interesting in the future

# In[48]:


# Save to JSON with timing, Notice that due to capacity this code was run on a separate script and is loaded here as input
json_path = output_large / "ec_records.json"
print(f"Starting JSON save to {json_path}...")
start_time = time.time()

try:
    with open(json_path, 'w') as f:
        json.dump(ec_records, f)

    end_time = time.time()
    elapsed = end_time - start_time
    size_mb = os.path.getsize(json_path) / 1024 / 1024

    print(f"Successfully saved to {json_path} in {elapsed:.2f} seconds ({size_mb:.2f} MB)")
except Exception as e:
    print(f"Error saving to JSON: {e}")


# In[49]:


# display whole column
pd.set_option("display.max_colwidth", True)
#pd.reset_option('display.max_colwidth')


# In[50]:


# Reading back the parquet as a dataframe
ECcontri_Uniprot_path = output_large / 'ECcontri_Uniprot_pathway.parquet'
ECcontri_Uniprot = pd.read_parquet(ECcontri_Uniprot_path)
print(f"DataFrame loaded from {ECcontri_Uniprot_path} with shape {ECcontri_Uniprot.shape}")
print(f"Memory usage after loading: {ECcontri_Uniprot.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
#32.97


# In[51]:


## opening ec_records dict
json_path = output_large / "ec_records.json"
with open(json_path, 'r') as f:
    ec_records = json.load(f)
ec_records_df = pd.DataFrame(ec_records)
ec_records_df = ec_records_df.rename(columns= {'ec_number': 'ec'})
#print(ec_records_df.head(2))
print(f"Memory usage after loading: {ec_records_df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")


# In[52]:


# Create a lookup from (protein_name, ec) to idx
ec_records_df.rename(columns= {'ec_number': 'ec'}, inplace=True)
idx_map = ECcontri_Uniprot.set_index(['protein_name', 'ec'])['idx'].to_dict()
# Map idx onto ec_records
ec_records_df['idx'] =ec_records_df.apply(
    lambda row: idx_map.get((row['protein_name'], row['ec']), pd.NA),
    axis=1
)
# order by idx
ec_records_df = ec_records_df.set_index('idx').sort_index().reset_index(drop=False)
# fill nans with 'Unknown'
ec_records_df['idx'] = ec_records_df['idx'].fillna('unknown').astype(str)


# In[53]:


# Informative renaming of the columns pendiente cambiar nombres en parse y main no?
ec_records_df = ec_records_df.rename(columns= {'functional_matches_detailed': 'functional_combi', 
        'functional_categories_detected': 'functional_sub', 'functional_child_terms': 'functional_child',  'synergy_child_terms_found': 'synergy_child',
        'synergy_categories_involved': 'synergy_sub','synergy_description': 'synergy_description', 'corrosion_mechanisms': 'mechanisms_child',
         'corrosion_mechanism_categories': 'mechanisms_sub','operational_environmental_factors': 'operational_sub', 'operational_matches_detailed': 'operational_combi'})


# In[54]:


# explode records
ec_records_df["consolidated_metals"].explode().unique()


# ### Cleaning consolidated_metals

# In[55]:


def clean_consolidated_metals(metals_list):
    """
    Post-process consolidated_metals list to fix remaining issues.
    Handles both string 'nan' and pandas NA objects.
    """
    if not isinstance(metals_list, list):
        return []

    # Canonical charge states
    CANONICAL = {
        'fe': 'Fe+2', 's': 'SO4-2', 'cu': 'Cu+', 'co': 'Co+2', 
        'cd': 'Cd+2', 'mn': 'Mn+2', 'ni': 'Ni+2', 'zn': 'Zn+2', 
        'pb': 'Pb+2', 'hg': 'Hg+2', 'al': 'Al+3', 'cr': 'Cr+3', 
        'mg': 'Mg+2', 'ca': 'Ca+2', 'ba': 'Ba+2', 'na': 'Na+', 
        'k': 'K+', 'h': 'H+', 'f': 'F-', 'cl': 'Cl-'
    }

    # Elements to exclude
    EXCLUDE = {'as', 'se', 'sr', 'pt', 'au', 'ag', 'w', 'li', 'u'}

    cleaned = []
    seen = set()

    for item in metals_list:
        # CRITICAL: Check for pandas NA/nan FIRST, before type checking
        try:
            if pd.isna(item):
                continue
        except:
            pass

        # Skip non-strings
        if not isinstance(item, str):
            continue

        item = item.strip()

        # Skip empty or string 'nan'
        if not item or item.lower() in ('nan', 'none', 'na', ''):
            continue

        # Normalize V5+ to V+5
        if item == 'V5+':
            item = 'V+5'

        # Extract first element from compound to check exclusions
        first_elem_match = re.match(r'^([A-Z][a-z]?)', item)
        if first_elem_match:
            first_elem = first_elem_match.group(1).lower()
            if first_elem in EXCLUDE:
                continue

        # Check if it's a base element that needs canonical form
        item_lower = item.lower()
        if item_lower in CANONICAL:
            item = CANONICAL[item_lower]

        # Add if not duplicate (case-insensitive)
        item_check = item.lower()
        if item_check not in seen:
            seen.add(item_check)
            cleaned.append(item)

    return cleaned

# Apply cleaning
ec_records_df['consolidated_metals'] = ec_records_df['consolidated_metals'].apply(clean_consolidated_metals)


# ______________________
# ## 1.3. Flattening ec_records
# _____________________

# In[56]:


def targeted_flatten_ec_records(df):
    """
    Single flat DataFrame.
    - protein_name fallback from enzyme_names (if protein_name invalid)
    - keep list columns as lists (no long tables, no joining)
    - clean list fields (strip, drop '', dedupe preserving order)
    - drop ko/metal_binding_info/coordination/pattern fields
    """
    df = df.copy()

    def _ensure_list(x):
        if x is None or (isinstance(x, float) and pd.isna(x)): return []
        if isinstance(x, (list, tuple, set)): return list(x)
        try:
            if isinstance(x, (np.ndarray, pd.Series)): return list(x.tolist())
        except Exception:
            pass
        return [x]

    def _clean_list(seq):
        seen, out = set(), []
        for item in _ensure_list(seq):
            # Skip None/NaN
            if item is None or (isinstance(item, float) and pd.isna(item)):
                continue

            # Keep numbers as numbers
            if isinstance(item, (int, float)):
                if item not in seen:
                    seen.add(item)
                    out.append(item)
                continue

            # Convert and clean strings
            s = str(item).strip()
            if not s or s.lower() in {'nan', 'none'}:
                continue
            if s not in seen:
                seen.add(s)
                out.append(s)
        return out

    GARBAGE_TERMS = {
    'nan', 'none', 'uncharacterized protein', 'hypothetical protein',
    'transferred to', 'deleted', 'obsolete', 'probable', 'putative', 'possible'
    }

    def _first_valid_name(seq):
        for x in _ensure_list(seq):
            s = str(x).strip().lower()
            if not s:
                continue
            # Skip if exact match to garbage
            if s in GARBAGE_TERMS:
                continue
            # Skip if starts with garbage prefix
            if any(s.startswith(prefix) for prefix in ['transferred to', 'deleted', 'obsolete', 'see ']):
                continue
            # Return the ORIGINAL case version, not lowercased
            return str(x).strip()
        return ''

    def _needs_fallback(name):
        s = (name or '').strip().lower()
        return (len(s) < 3) or ('uncharacterized' in s) or ('hypothetical' in s)

    # protein_name fallback from enzyme_names
    if 'protein_name' not in df.columns:
        df['protein_name'] = ''
    if 'enzyme_names' in df.columns:
        mask = df['protein_name'].fillna('').map(_needs_fallback)
        df.loc[mask, 'protein_name'] = df.loc[mask, 'enzyme_names'].map(_first_valid_name)

    # ---- columns to keep/drop ----
    keep_cols = [
        'ec',                               # strings
        'protein_name',                     # strings
        'idx',                              # strings
        'enzyme_class',                     # strings
        'enzyme_names',                     # list this should be just the first one and discard the others

        #'functional_categories'            # dict sub+child +scoring
        'functional_combi',              # dict sub+child (preserve provenance/alignment), that is just to know where the sub and children come from
        'functional_sub',                   # list subcategories (fc) that was suppose to be just the subcategories terms but it keeps the list of subcategories, so we can explode at least at this level.
        'functional_child',                 # list child terms (fc) to explode the children will be also a good idea, althought not sure if incredible long.

        'synergy_child',                    # list child terms (synergies) same as before
        'synergy_sub',                      # list subcategories (synergies)
        'synergy_description',                 # dict origin/type of synergy

        'consolidated_metals',              # joined list brenda_metals it can not be exploded because makes less sense since metals are ubiquos in these samples
        'pathways_db',                      # list , this column would be removed, so no important
        'reaction_db',                      # list , this column would be removed, so no important
        'hierarchy',                        # list hierarchies, this could be exploded 
        'operational_sub',                  #  list subc terms this are possible to explode too but too long also
        'operational_combi',                # dict sub + child this is another list to explode
        'mechanisms_child',
        'mechanisms_sub',

        'overall_functional_score',         # float
        'overall_metal_score',              # float
        'overall_synergy_score',            # float
        'corrosion_relevance_score',        # float
    ]

    drop_explicit = {
        'ko','metal_binding_info','residue_name','coordination_number', 'corrosion_mechanism_child_terms','detected_pathways', 'pathway_matches_detailed',
        'functional_categories', 'metal_categories_detected','coordination_patterns','residue_binding', 'brenda_metals', 'detected_metal_categories',
        'metal_matches_detailed','metal_count','legacy_synergy_score','legacy_synergy_groups', 'fc_cooccurrence_synergy_hit','synergy_type'
        }
    # drop columns that are explicitly unwanted    
    df.drop(columns=[c for c in drop_explicit if c in df.columns], inplace=True, errors='ignore')

    present_keeps = [c for c in keep_cols if c in df.columns]
    flat = df.loc[:, present_keeps].copy()

    # prepare the lists to be exploded and explode the following columns for granularity
    explode_cols = [
        'functional_sub',    # subcategories (fc)
        'functional_child',  # child terms (fc)
        'mechanisms_sub',    # corrosion mechanisms
        'mechanisms_child',   # corrosion mechanisms
        'hierarchy',         # hierarchies
        'operational_sub',  # subcategories (operational)
        'operational_combi'
    ]
    # Clean as before
    for c in explode_cols:
        if c in flat.columns:
            flat[c] = flat[c].apply(_clean_list)
    # Now explode each column one by one (multiple explosion is not always supported in pandas)
    for c in explode_cols:
        if c in flat.columns:
            flat = flat.explode(c)

    # Cleaning the other list-like columns prefered to keep as lists (not explode)
    list_like_keep = [c for c in [
        'enzyme_names', 'consolidated_metals', 'pathways_db', 'reaction_db', 'synergy_child', 'synergy_sub'
    ] if c in flat.columns]
    for c in list_like_keep:
        flat[c] = flat[c].apply(_clean_list)

    # ensure dict columns are dicts (not NaN)
    for dcol in ['functional_combi','synergy_description','operational_combi']:
        if dcol in flat.columns:
            flat[dcol] = flat[dcol].apply(lambda x: x if isinstance(x, dict) else ({} if x is None or (isinstance(x,float) and pd.isna(x)) else x))

    return flat


# In[57]:


# display whole column
pd.set_option("display.max_colwidth", True)
#pd.reset_option('display.max_colwidth')


# In[58]:


ec_record_flat = targeted_flatten_ec_records(ec_records_df)


# In[59]:


# throw away garbage ec_record_df
del ec_records_df
del idx_map
gc.collect()


# _____________________
# __Enzyme and Protein Ocurrencies__
# ________________

# In[60]:


# Count occurrences of each value in both columns
enzyme_counts = ec_record_flat['enzyme_names'].explode().value_counts()
protein_counts = ec_record_flat['protein_name'].value_counts()
ec_counts = ec_record_flat['ec'].value_counts()
# Get the number of unique values in each column
enzyme_unique_count = ec_record_flat['enzyme_names'].explode().nunique()
protein_unique_count = ec_record_flat['protein_name'].nunique()
ec_unique_count = ec_record_flat['ec'].nunique()

# Print the results
print("EC Numbers - Unique Values:", ec_unique_count)
print("\nProtein Names - Unique Values:", protein_unique_count)
print("Protein Names - Top 5 Occurrences:\n", protein_counts.head())
print("\nEnzyme Names - Unique Values:", enzyme_unique_count)
print("Enzyme Names - Top 5 Occurrences:\n", enzyme_counts.head())


# EC Numbers - Unique Values: 1793
# 
# Protein Names - Unique Values: 1617
# Protein Names - Top 5 Occurrences:
#  protein_name
#                               79
# 23s rrna-methyltransferase    15
# 16s rrna-methyltransferase    11
# trna-methyltransferase         9
# alcohol-dehydrogenase          4
# Name: count, dtype: int64
# 
# Enzyme Names - Unique Values: 9599
# Enzyme Names - Top 5 Occurrences:
#  enzyme_names
# 23s rrna-methyltransferase         15
# 16s rrna-methyltransferase         11
# trna-methyltransferase             10
# hydrogenase                         5
# restriction-modification system     4
# Name: count, dtype: int64

# _________________________
# ## 1.4. Integrating Reaction_Primary
# ________________________
# database coming from section 7 of notebook 6_picrust_functional.ipynb that has the reaction data saved in chunks

# In[61]:


# Reading loading 10 min
output_path = output_large / "chunks_react"
chunk_files = sorted([os.path.join(output_path, f) for f in os.listdir(output_path) if f.startswith('chunk_') and f.endswith('.pkl')])
ECcontri_react = pd.concat((pd.read_pickle(f) for f in chunk_files), ignore_index=True)
# Extract the number from 'site_X'
ECcontri_react['sample_num'] = ECcontri_react['Sites'].str.extract('site_(\d+)').astype('int8')
# Sort by the numerical part
ECcontri_react = ECcontri_react.sort_values(by='sample_num')
# Drop the temporary column
ECcontri_react = ECcontri_react.drop(columns=['sample_num'])
ECcontri_react = ECcontri_react.drop_duplicates(subset=["ec", "Sites", "Genus"])
ECcontri_react = ECcontri_react.copy()
ECcontri_react["Sites"] = ECcontri_react["Sites"].astype("category")
ECcontri_react["Genus"] = ECcontri_react["Genus"].astype("category")
ECcontri_react["rxn"] = ECcontri_react["rxn"].astype("category")
ECcontri_react = ECcontri_react.rename(columns = {"rxn": "react_primary"})


# In[62]:


pre_ECcontri_enriched = pd.merge(ECcontri_Uniprot, ECcontri_react, on=['Sites', 'Genus', 'ec'], how='left')


# In[63]:


# throw away garbage ec_record_df
del ECcontri_react
del ECcontri_Uniprot
gc.collect()


# In[64]:


# List of EC numbers we want
wanted_ecs = ['aspartate kinase', 
    'thioredoxin-disulfide-reductase', 
    'propionate-coa ligase']#['6.2.1.17', '1.8.1.9' , '2.7.2.4']

# Filter rows where ec is in wanted_ecs
columns_to_show = ['idx', 'protein_name', 'pathway_primary', 'react_primary']

filtered_df = pre_ECcontri_enriched[pre_ECcontri_enriched['protein_name'].isin(wanted_ecs)]
print(f"pre_ECcontri_enriched: {filtered_df[columns_to_show]}")


# ______________________________________________
# # 2. Analysis pre-Enriching ECcontri_Uniprot
# ______________________________________________
# 

# Reactions and pathways are entities that came in the results of picrust2 in the form of Picrus_results df it was composed of progamatic pathways(ipath) and desciptive pathways(pathways), they were 366. However the data also have 3000 Reactions and a parsing dictionary between pathways and reactions. So the pathways obtained through this parcing (path_parce 545). Unfortunately the nature of the dbs available in this study are open access and the programatic pathways are no accesible, that is why the descriptors(pathways) were taken going forward. This pathways column, was taken as backbone and the gaps were filled up with the results from path_parce. This column then contains pathways descriptive but also programatic(pathway_primary). In order to make the allocation of the categories of this pathways for the dictionaries on global_terms.py inside pathways and functional categories, a manual allocation was imposible due to the sheer amount of data (545 names). Therefore, several codes were developed to make the distribution accurate (validation notebook in scoring_system), however the results were rudimentary. Ultimately several language models were used to allocate the columns pathways into the different dictionaries, mainly Gemini and Copilot, and ultimately human review.
# ______________________
# ### Auditor Structure 
# ______________________
# the pathway_primary and reaction_primary were used as golden truth and it was keep with no further updating because corresponds to a granular nature of the results from picrust and inputing more data to fillup the gaps would imply diluting the good data and adding noise. The pathway_primary/reaction_primary column was use to validate the retrieved data pathways_db/react_db by objectively measuring the agreement between the two. This helps to improve the scoring system. 

# This section implements a comparison between pathway and reaction annotations from two sources for each EC (Enzyme Commission) number. The primary DataFrame (pre_ECcontri_enriched) provides a fine-grained annotation for individual proteins, sites, and genera, with pathway and reaction assignments in the columns pathway_primary and react_primary. The secondary DataFrame (ec_record_flat) mines the pathways/reactions from brenda and Kegg databases and concatenate them those entities are storaged with all the other metadata on the ec_records dictionary as lists terms per EC. To enable a row-wise comparison, the code first transforms the data to EC-protein level into a long, exploded format, such that each pathway or reaction appears in its own row.

# In[65]:


def counting_df(df, df_name="DataFrame"): 
    print(f"Results for '{df_name}':")
    # count nans on ec and protein_name
    nan_ec_count = df['ec'].isna().sum()
    nan_protein_name_count = df['protein_name'].isna().sum()
    print(f"Number of NaNs in 'ec': {nan_ec_count}")
    print(f"Number of NaNs in 'protein_name': {nan_protein_name_count}")
    # count unique values on ec and protein_name
    unique_ec_count = df['ec'].nunique(dropna=True)
    unique_protein_name_count = df['protein_name'].nunique(dropna=True)
    print(f"Number of unique values in 'ec': {unique_ec_count}")
    print(f"Number of unique values in 'protein_name': {unique_protein_name_count}")
    #count how many pairs ec protein name are empty
    empty_pairs_count = df[df['ec'].isna() & df['protein_name'].isna()].shape[0]
    print(f"Number of rows with both 'ec' and 'protein_name' empty: {empty_pairs_count}")
    # count the nans on idx column
    nan_idx_count = df['idx'].isna().sum()
    print(f"Number of NaNs in 'idx': {nan_idx_count}")
    shape = df.shape
    print(f"'{df_name} shape: {shape}")

# Usage:
counting_df(ec_record_flat, df_name="ec_record_flat")


# In[66]:


counting_df(pre_ECcontri_enriched, df_name="pre_ECcontri_enriched")


# In[67]:


ec_record_flat['idx'] = pd.to_numeric(ec_record_flat['idx'], errors='coerce').astype('Int32')


# In[68]:


ec_pathway_db  = ec_record_flat[['ec', 'protein_name','pathways_db', 'idx']].explode('pathways_db').set_index("idx").sort_index().reset_index(drop=False)
ec_path_contry = pre_ECcontri_enriched[["ec", 'protein_name', "pathway_primary", "idx"]].set_index("idx").sort_index().reset_index(drop=False)


# In[69]:


ec_pathways = pd.merge(ec_path_contry, ec_pathway_db, on=['ec', 'protein_name', 'idx'], how ='left').set_index('idx').sort_index().reset_index(drop=False)
ec_pathways.sample(10).set_index("idx").sort_index()


# In[70]:


# remove the idx column from ec_record_flat without slicing it
ec_record_flat = ec_record_flat.drop(columns=['idx'])
del ec_pathway_db
del ec_path_contry
del ec_pathways
gc.collect()


# Primary Pathway Column (More Accurate)
# 
# Analysing the two columns it can be seen that the pathway_primary which represents BioCyc/MetaCyc pathway annotations, are more acurate and biologically relevant.  Highly curated and specific: Each enzyme is mapped to its most direct, well-characterized pathway
#     Biochemically precise:
#         PWY-5751 for alcohol dehydrogenase represents a specific, defined pathway
#         "Glycolysis III (from glucose)" is the exact variant of glycolysis
#         "Superpathway of L-isoleucine biosynthesis I" is a precise metabolic reconstruction
# 
# Pathways_DB Column (Less Accurate)
# 
# The databases mined pathways come from 2 different sources from Brenda and KEGG pathway annotations, which show several problematic patterns:
# Overly Broad Associations: Alcohol dehydrogenase is linked to several different pathways, many quite remote from its primary function
#     Questionable Connections:
#         ADH in "Chloroalkane degradation" - very specialized, not universal
#         ADH in "Drug metabolism - cytochrome P450" - conflates different enzyme systems
#         ADH in "Naphthalene degradation" - highly organism-specific
#     Generic Catch-all Categories: "Biosynthesis of secondary metabolites" appears for all three enzymes, which is too broad to be meaningful
#     Missing Specificity: While some connections are valid (like GAPDH in glycolysis), the lack of pathway variants/specificity reduces precision
# 
# Valid KEGG Connections:
# 
#     ADH in "Glycolysis/Gluconeogenesis" ✓
#     ADH in "Retinol metabolism" ✓
#     GAPDH in "Glycolysis/Gluconeogenesis" ✓
#     GAPDH in "Carbon fixation by Calvin cycle" ✓
# 
# Conclusion: Because of this comparison the first plan to validate the mined pathways (ec_records) with the picrust_results (ECcontries) is not valid.
# 
# The primary pathway column is significantly more accurate because it provides specific well-defined pathway assignments, avoids over-annotation with tangential associations. It reflects direct biochemical roles rather than peripheral involvement.
# 
# The pathways_db/reaction_db columns suffer from the common KEGG issue of casting too wide a net, leading to many false or marginal associations that dilute the core biochemical function of each enzyme.

# # 3. Enrichment
# The ECcontri dataset will serve as the reference standard. The following columns will be retained:'idx', 'Sites', 'Genus', 'norm_abund_contri', 'protein_name', 'ec', 'Category', 'pathway_primary', 'react_primary'. All rows in ECcontri possess an EC number; however, approximately 272,716 out of 1.5 million entries lack an associated protein_name (NaN). For these cases, protein_name will be supplemented from the corresponding entry in ec_record; if unavailable, the first enzyme name listed in ec_record will be used. For pathway_primary, all ECcontri values will be preserved, and any missing entries will be complemented by pathway_db. The same imputation strategy will be applied to react_primary, utilizing reaction_db as a fallback. ECcontri will be then enriched with supplementary data from ec_record where applicable.

# In[71]:


def enrich_eccontri_data(eccontri_df, ec_records_flat):
    """ Enriches eccontri_df with metadata from ec_records_flat, with fallback logic for pathway/reaction/enzyme names. """
    # 0. Prep: ensure string columns are clean
    eccontri = eccontri_df.copy()
    meta = ec_records_flat.copy()
    eccontri["protein_name"] = eccontri["protein_name"].astype(str).replace('nan', '')
    meta["protein_name"] = meta["protein_name"].astype(str).replace('nan', '')
    meta['enzyme_names'] = meta['enzyme_names'].apply(
        lambda x: x[0] if isinstance(x, (list, tuple)) and len(x) > 0 else x
    )

    # 1. c
    metadata_columns = ['enzyme_class','functional_combi','functional_sub','functional_child',
        'synergy_child','synergy_sub','synergy_description','consolidated_metals',
        'pathways_db','reaction_db','hierarchy','mechanisms_sub','mechanisms_child','operational_sub',
        'operational_combi','overall_functional_score','overall_metal_score',
        'overall_synergy_score','corrosion_relevance_score', 'enzyme_names'
    ]
    present_meta_cols = [c for c in metadata_columns if c in meta.columns]
    #=======================
    # ec_records_flat is flattened, so take first occurrence of each (ec, protein_name)
    # The aggregate fields (functional_combi, synergy_data) are identical across duplicates
    meta_unique = meta.drop_duplicates(subset=['ec', 'protein_name'], keep='first')

    # Merge on (ec, protein_name) - this preserves ALL eccontri rows
    merged = pd.merge(eccontri, meta_unique[['ec', 'protein_name'] + present_meta_cols], 
                     on=['ec', 'protein_name'], how='left', suffixes=('', '_meta'))

    # ============================
    '''merged = pd.merge(eccontri, meta[['idx'] + present_meta_cols], on='idx', how='left', suffixes=('', '_meta'))
    '''
    # 2. Fallback logic for protein_name, pathway_primary, react_primary
    # protein_name: if missing, try protein_name_meta, then enzyme_names
    def fallback_protein(row):
        pn = str(row['protein_name'])
        if not pn or pn in {"uncharacterized protein", "hypothetical protein", "nan", "none"}:
            # Try protein_name_meta first
            pn2 = str(row.get('protein_name_meta', ''))
            if pn2 and pn2 not in {"uncharacterized protein", "hypothetical protein", "nan", "none"}:
                return row['protein_name_meta']
            # Then fallback to enzyme_names
            ens = row.get('enzyme_names', '')
            if isinstance(ens, str) and ens and ens.lower() not in {"uncharacterized protein", "hypothetical protein", "nan", "none", "transferred to", "deleted", "obsolete"}:
                return ens
        return row['protein_name']

    merged['protein_name'] = merged.apply(fallback_protein, axis=1)

    # pathway_primary: if missing, fill from pathways_db (if present and is a list/string)
    def fallback_pathway(row):
        curr = str(row.get('pathway_primary', ''))
        pws = row.get('pathways_db', None)
        if not curr or pd.isna(curr):
            if isinstance(pws, (list, tuple)) and len(pws) > 0:
                return pws[0]
            elif isinstance(pws, str) and pws:
                return pws
        return curr

    if 'pathway_primary' in merged.columns:
        merged['pathway_primary'] = merged.apply(fallback_pathway, axis=1)

    # react_primary: if missing, fill from reaction_db
    def fallback_react(row):
        curr = str(row.get('react_primary', ''))
        rct = row.get('reaction_db', None)
        if not curr or pd.isna(curr):
            if isinstance(rct, (list, tuple)) and len(rct) > 0:
                return rct[0]
            elif isinstance(rct, str) and rct:
                return rct
        return curr

    if 'react_primary' in merged.columns:
        merged['react_primary'] = merged.apply(fallback_react, axis=1)

    # 3. Final columns 
    base_cols = [c for c in ['idx','Sites','Genus','norm_abund_contri','protein_name','ec','Category','pathway_primary','react_primary'] if c in merged.columns]
    final_cols = base_cols + [c for c in present_meta_cols if c not in base_cols]
    enriched_df = merged[final_cols].copy()

    return enriched_df


# In[72]:


ECcontri_Uniprot_enriched = enrich_eccontri_data(pre_ECcontri_enriched, ec_record_flat)


# In[73]:


# drop columns
ECcontri_Uniprot_enriched= ECcontri_Uniprot_enriched.drop(columns=['ec','pathways_db', 'reaction_db', 'enzyme_names']) # used for enrichment purposes and lookups no longer needed
del pre_ECcontri_enriched
del ec_record_flat
gc.collect()


# ## 3.1. Checking data consistency

# In[74]:


def counting_abun_df(df, df_name="DataFrame"): 
    print(f"Results for '{df_name}':")
    # count nans on ec and protein_name
    nan_abundance_count = df['norm_abund_contri'].isna().sum()
    nan_protein_name_count = df['protein_name'].isna().sum()
    print(f"Number of NaNs in 'norm abundance column': {nan_abundance_count}")
    print(f"Number of NaNs in 'protein_name': {nan_protein_name_count}")
    # count unique values on ec and protein_name
    unique_abundance_count = df['norm_abund_contri'].nunique(dropna=True)
    unique_protein_name_count = df['protein_name'].nunique(dropna=True)
    print(f"Number of unique values in 'norm_abund_contri': {unique_abundance_count}")
    print(f"Number of unique values in 'protein_name': {unique_protein_name_count}")
    #count how many pairs ec protein name are empty
    empty_pairs_count = df[df['norm_abund_contri'].isna() & df['protein_name'].isna()].shape[0]
    print(f"Number of rows with both 'norm_abund_contri' and 'protein_name' empty: {empty_pairs_count}")
    # count the nans on idx column
    nan_idx_count = df['idx'].isna().sum()
    print(f"Number of NaNs in 'idx': {nan_idx_count}")
    shape = df.shape
    print(f"'{df_name} shape: {shape}")

# Usage:
counting_abun_df(ECcontri_Uniprot_enriched, df_name="ECcontri_Uniprot_enriched")


# In[75]:


# List of EC numbers we want
wanted_ecs = ['aspartate kinase', 
    'thioredoxin-disulfide-reductase', 
    'propionate-coa ligase']#['6.2.1.17', '1.8.1.9' , '2.7.2.4']

# Filter rows where ec is in wanted_ecs
columns_to_show = ['protein_name', 'enzyme_class', 
    'hierarchy', 'consolidated_metals','functional_combi', 'functional_sub', 'functional_child',
       'synergy_child', 'synergy_sub', 'synergy_description',
        'operational_sub', 'operational_combi', 'mechanisms_child',
       'mechanisms_sub'
]

filtered_df = ECcontri_Uniprot_enriched[ECcontri_Uniprot_enriched['protein_name'].isin(wanted_ecs)]
print(f"ECcontri_Uniprot_enriched: {filtered_df[columns_to_show]}")


# ### Sanitising generic types to pandas types

# In[76]:


# Fix columns with 0-d arrays or numpy generics to pandas native types
def sanitize_for_json(val):
    # Convert 0-d numpy arrays to their item
    if isinstance(val, np.ndarray):
        if val.ndim == 0:
            return val.item()
        else:
            return val.tolist()
    # Convert numpy scalars
    if isinstance(val, (np.generic,)):
        return val.item()
    # Leave lists, dicts, strings, ints, floats alone
    return val

# Apply to all cells
for col in ECcontri_Uniprot_enriched.columns:
    ECcontri_Uniprot_enriched[col] = ECcontri_Uniprot_enriched[col].map(sanitize_for_json)


# In[77]:


counting_abun_df(ECcontri_Uniprot_enriched, df_name="ECcontri_Uniprot_enriched")


# In[78]:


# Count occurrences of each value in both columns
metal_counts = ECcontri_Uniprot_enriched['consolidated_metals'].explode().value_counts()
protein_counts = ECcontri_Uniprot_enriched['protein_name'].value_counts()
reactions_counts =ECcontri_Uniprot_enriched['react_primary'].value_counts()
# Get the number of unique values in each column
metal_unique_count = ECcontri_Uniprot_enriched['consolidated_metals'].explode().nunique()
protein_unique_count = ECcontri_Uniprot_enriched['protein_name'].nunique()
reactions_unique_count = ECcontri_Uniprot_enriched['react_primary'].nunique()

# Print the results
print("Metal - explode value counts:", metal_counts)
print("\nProtein Names - value counts:", protein_counts)
print("Protein Names - Top 5 Occurrences:\n", protein_counts.head())
print("\nmetal consolidated - Unique Values:", metal_unique_count)
print("reactions_unique_count- Top 5 Occurrences:\n", reactions_unique_count)


# ### NOTE
# This notebook has been done in several parts and keep together for thematic continuity. First part: an attempt was made to install picrust2 through several ways on vscode local including using virtual machine, but the capacity was no sufficient. For Colab it was attempted to use anaconda but it failed due to compatibility amongst other pittfails. Ultimately Picrust2 algoritm for predicting functional metabolic pathways was executed using Galaxy Europe platform. Second part:  API calls were performed to retrieve the corresponding protein names from UniProt based on the EC numbers and OTU-abundance data generated by PICRUSt2. Due to the resource-intensive nature of this task, it was conducted on Colab Pro rather than locally. Third part:  A compendium database, ec_records, was created, compiling key information from authoritative databases on proteins associated with corrosion. This task required significant storage capacity and is a resource intense task it was carried out on Kaggle. Four part was done also in kaggle and comprise the flattening of the ec_records and enrichment of the uniprot data to obtain ECcontri_Uniprot_enriched. Fifth part was to prepare the data for the filtering pipeline "analyse_corrosion_protein". This last part was done in vscode local but also in Kaggle, plataform that can hold the whole of the data. Datasets are found in Kaggle, together with the two scripts to make ec_records and eccontri as well as in Github Notebook 6 of this repository. Kaggle is use due to constrains on Github capacity of storage large files. 
# Each part is resource intensive, hence it is recomended to perform one at the time. However, those recomendations should be no necesary, then I recognise that modern machines can do most of the task on brevity and no one would repeate this task wich such scarcity of resources as the present work. 
# Fifh part started on section 9.6, the preparation and continues from this point onward:

# ## 3.2. Preparing data to enter next pipeline
# Columns formating cleaning and dtypes conversion

# In[79]:


#eccontri_path = output_large / 'ECcontri_Uniprot_enriched.json' # I can not longer manage to do it on a reasonable time
#ECcontri_Uniprot_enriched = pd.read_json(eccontri_path, orient='records')
# size of the df
print(f"Loaded enriched DataFrame from {ECcontri_Uniprot_enriched} with shape {ECcontri_Uniprot_enriched.shape}")   
print(f"Memory usage after loading: {ECcontri_Uniprot_enriched.memory_usage(deep=True).sum() / 1024**2:.2f} MB")


# This function prepares the columns to enter the next pipeline by converting the dict into subcategories and parentesis children, the list to semicolons separated list and the special column synergy to simple form. The next function improve datatypes for some columns

# In[80]:


def optimize_dtypes(df):
    """Convert columns to memory-efficient dtypes"""
    # Categorical columns (low cardinality strings) no useful for the present columns
    # Float columns to float32
    float_columns = ['overall_functional_score', 'overall_metal_score',
                     'overall_synergy_score', 'corrosion_relevance_score']
    int_columns = ['idx']  
    cat_columns = ['react_primary', 'enzyme_class','hierarchy']
    int8_columns = ['Category']  
    str_columns = ['protein_name', 'Genus', 'Sites', 'protein_name', 'Genus', 'Sites', 'functional_sub', 'functional_child',
               'hierarchy', 'pathway_primary', 'mechanisms_sub', 'mechanisms_child', 'operational_sub']
    for col in int_columns:
        if col in df.columns:
            df[col] = df[col].astype('Int32')
    for col in int8_columns:
        if col in df.columns:
            df[col] = df[col].astype('Int8')
    for col in str_columns:
        if col in df.columns:
            df[col] = df[col].astype('string')
    for col in float_columns:
        if col in df.columns:
            df[col] = df[col].astype('float32')
    for col in cat_columns:
        if col in df.columns:
            df[col] = df[col].astype('category')    
    return df

ECcontri_Uniprot_Enriched = optimize_dtypes(ECcontri_Uniprot_enriched)


# In[81]:


counting_abun_df(ECcontri_Uniprot_Enriched, df_name="ECcontri_Uniprot_Enriched")


# ### Savind ECcontri Uniprot Enriched

# In[82]:


# Save ECcontri uniprot Enriched to use later
ECcontri_Uniprot_path = output_large / 'ECcontri_Uniprot_Enriched.parquet'
ECcontri_Uniprot_Enriched.to_parquet(ECcontri_Uniprot_path, engine='pyarrow', compression='snappy')
# size df
print(f"Saved enriched DataFrame: {os.path.getsize(ECcontri_Uniprot_path) / (1024**2):.2f} MB")


# In[83]:


# del ECcontri_Uniprot_enriched
gc.collect()


# In[84]:


# starting from here with ECcontri uniprot Enriched dataframe
ECcontri_Uniprot_path = output_large / 'ECcontri_Uniprot_Enriched.parquet'
ECcontri_Uniprot_Enriched = pd.read_parquet(ECcontri_Uniprot_path, engine='pyarrow')
print(f"Read enriched DataFrame: {os.path.getsize(ECcontri_Uniprot_path) / (1024**2):.2f} MB")


# In[85]:


knee_path = output_base / "genus_to_threshold.csv"
with open(knee_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    genus_to_threshold = next(reader) 


# In[86]:


# List of EC numbers we want
wanted_ecs = ['aspartate kinase', 
    'thioredoxin-disulfide-reductase', 
    'propionate-coa ligase']#['6.2.1.17', '1.8.1.9' , '2.7.2.4']

# Filter rows where ec is in wanted_ecs
columns_to_show = ['protein_name', 'enzyme_class', 
    'hierarchy', 'consolidated_metals','functional_combi', 'functional_sub', 'functional_child',
       'synergy_child', 'synergy_sub', 'synergy_description',
        'operational_sub', 'operational_combi', 'mechanisms_child',
       'mechanisms_sub'
]

filtered_df = ECcontri_Uniprot_enriched[ECcontri_Uniprot_enriched['protein_name'].isin(wanted_ecs)]
print(f"ECcontri_Uniprot_enriched: {filtered_df[columns_to_show]}")


# # 4. Corrosion Protein Filtering pairs Bacteria-Protein by significance to the risk category
# 
# The analyze_corrosion_proteins function establishes a systematic framework for identifying protein-genus pairs associated with corrosion across different risk categories. The function prepares the data by converting enzyme records to a searchable dictionary and mapping sites to risk categories, then tracks which sites contain each protein-genus combination to maintain traceability throughout the analysis.
# After initial preparation, the function performs pattern recognition using prevalence, significance, and frequency metrics through the classify_abundance_patterns function. It then integrates biological metadata into the results, classifies proteins into housekeeping and niche-specific categories, and separates increasing patterns from others. Only increasing patterns are used to prioritize markers based on combined pattern and biological significance. These results are further refined using the knee abundance threshold (calculated in section 8.11) to select proteins of biological relevance. Finally, it organizes results into specialized groups for interpretation.
# The analysis uses absolute contribution to abundance (abund_contri) as this metric is comparable across samples and reflects biological influence without relativization. Inverse patterns are preserved separately as they potentially represent bacteria with protective capacities in the ecosystem or indicate bystander organisms whose protein content could provide insights about their ecological role and why they decrease as systems deteriorate.

# In[87]:


def analyze_corrosion_proteins(eccontri_df, per_genus_count=10, genus_to_threshold = genus_to_threshold):
    """
    Comprehensive analysis of protein-genus pairs for corrosion relevance

    Parameters:
    eccontri_df : pandas DataFrame The enriched ECcontri_Uniprot dataframe with EC, protein_name, Genus, abundance, etc.
    ec_records : List of dictionaries with EC, protein_name metadata with corrosion relevance information
    balance_genera : bool, Whether to balance representation across genera (default: False)
    per_genus_count : intNumber of markers to include per genus if balancing (default: 10)
    genus_to_threshold : dict, Optional mapping of genus to Knee_abundance thresholds for filtering (default: None)
    pattern_result : pd.DataFrame, precomputed pattern result DataFrame to bypass compute-intensive steps.

    Returns:
     results :  Dictionary containing various analysis results
    """
    print("Starting corrosion protein analysis...")
    eccontri_df = eccontri_df.copy().reset_index()

    print(f"Analyzing {len(eccontri_df)} data points across {len(eccontri_df['Sites'].unique())} Sites...")
    if 'index' in eccontri_df.columns:
        eccontri_df= eccontri_df.drop(columns="index")

    eccontri_df["idx"] = eccontri_df["idx"].astype('int32')
    eccontri_df = eccontri_df.reset_index(drop=False)

    # Fist insert the Frequency in the original df so that it be granular with transform, the idea is to preserve the info that is droped with protein-genus drop duplicates
    eccontri_df['Frequency'] = eccontri_df.groupby(['Sites', 'Genus', 'protein_name'], observed=True)['idx'].transform('size')
    eccontri_df['Frequency'] = eccontri_df['Frequency'].fillna(0).clip(lower=0).astype('int8')
    print("eccontri_df.columns immediately after creating column Frequency:", eccontri_df.columns.tolist())
    # this was originally an statistical analysis but it was no possible so it was changed the original p_value to significance and prevalence
    pattern_data = classify_abundance_patterns(eccontri_df)
    print("Pattern data created")

    # Merging pattern data with ec_records, which is the colection of databases knowledge on biological meaning
    integrated_results = merging_pattern_record(pattern_data, eccontri_df)
    print("Integrated results created")

    # Clasify ubiquitous, niche genus-protein names by mechanism, pathways. It was made an effort to prefer the universal pathways because all was being assigned to niche
    print("Classifying pathways by specificity...")
    classified_results = classify_pathways_by_specificity(integrated_results)

    # first separation of the data occur here, where the data is divided into positive and inverse patterns
    print("Separating positive and inverse patterns by patterns")
    increasing_results, inverse_results, constant_df = separate_by_pattern(classified_results)

    # Prioritize calculates scores focused only on the positive results,  namely increasing_results
    print("Prioritizing markers based on statistical and biological relevance...")
    prioritized_markers = prioritize_markers(increasing_results)

    # Create base df with only increasing_results and now only the genera which abundance value is above the knee inflexion, but also the per_genus_count to give the change to other genus to show up
    print(f"Balancing genus representation (top {per_genus_count} per genus)...")
    balanced_markers = balance_genus_representation(prioritized_markers, eccontri_df, genus_to_threshold =genus_to_threshold,
                                per_genus_count=per_genus_count)
    print(f"Created balanced dataset with {len(balanced_markers)} markers from {len(prioritized_markers['Genus'].unique())} genera")

    # Create marker groups from balanced markers, only expected categories are 2 and 3 for the risk label category and only increasing pattern at this point
    print("Creating specialized marker groups from balanced markers...")
    marker_groups = create_marker_groups(balanced_markers, top_count=200, threshold_percentile=0.75)

    print("Analysis complete!")

    results_dict = {
        'pattern_data': pattern_data,
        'integrated_results': integrated_results,
        'classified_results': classified_results,
        'increasing_markers': increasing_results,
        'prioritized_markers': prioritized_markers,
        'balanced_markers': balanced_markers,
        'marker_groups': marker_groups,
        'inverse_markers': inverse_results
    }

    return results_dict, marker_groups


# ## 4.1 Determining Abundance Patters
# The classify_abundance_patterns function implements a sophisticated classification system for protein-genus abundance profiles across risk categories, prioritizing patterns with potential corrosion relevance.
# This classification focuses on corrosion-relevant trends by first calculating mean abundance per category for each genus-protein pair. The function begins by capturing the frequency metrics from the original dataframe, where frequency represents the number of occurrences of each unique genus-protein pair across sites and categories.
# The algorithm tracks which categories (1, 2, and 3) have non-zero means for each protein-genus pair and applies detailed pattern recognition logic based on the relationships between these means. These detailed patterns are then grouped into broader categories ("increasing", "decreasing", "constant", or "other") to facilitate downstream analysis focusing on proteins that correlate positively with corrosion risk.
# Beyond pattern identification, the function calculates several quantitative metrics to prioritize patterns:
# 
# Fold changes between categories (with pseudocounts to avoid division by zero)
# Statistical significance through permutation testing with p-value calculation
# Pattern-specific p-value assignment for biologically relevant categories
# 
# Statistical significance assessment is implemented through a permutation-based approach that evaluates whether observed abundance patterns differ significantly from random chance. The algorithm computes a composite score combining presence proportion and abundance ratios, then performs 1000 random permutations of category labels to establish a null distribution. Special p-value assignments are made for patterns of particular biological interest: genera appearing only in category 2 (transitional conditions) receive p = 0.002, those exclusive to category 3 (severe corrosion) receive p = 0.001, and those present only in categories 2 and 3 (excluding normal conditions) receive p = 0.003. This approach prioritizes stress-responsive patterns while maintaining statistical rigor.
# 
# The frequency metric plays an important role in the final significance score calculation, where it's logarithmically scaled to give weight to commonly occurring protein-genus pairs without overly penalizing rare but potentially important cases. This approach balances the importance of frequently observed patterns (which may represent core functions) with unique patterns that might indicate specialized adaptive responses to corrosion conditions.The algorithm identifies distinct patterns
# with 2 transitions between 3 categories, we have 3² = 9 possible patterns:
# 
# "steadily_increasing" (cat1 < cat2 < cat3) - Clear positive correlation with risk  
# "early_plateau" (cat1 < cat2 = cat3) - Initial response that stabilizes  
# "peak_at_transition" (cat1 < cat2 > cat3) - Specific to transitional environments  
# "late_response" (cat1 = cat2 < cat3) - Only activated in high-risk environments  
# "risk_independent" (cat1 = cat2 = cat3) - Housekeeping/core function  
# "late_decline" (cat1 = cat2 > cat3) - Inhibited in severe conditions  
# "stress_recovery" (cat1 > cat2 < cat3) - Recovery pattern after initial stress  
# "adaptation_plateau" (cat1 > cat2 = cat3) - Adjustment to new baseline  
# "steadily_decreasing" (cat1 > cat2 > cat3) - Negative correlation with risk  
# 
# epsilon = 1e-5 → used only to judge whether two means are considered "equal" (floating point tolerance).  (are two means practically equal?)
# pseudocount = 1e-5 → used only to avoid division by zero in fold change calculations.  
# clip(-10, 10) → used only to cap extreme values of log2fc, not to define significance.  
# 
# Fold changes refer to how much the normalized abundance of a protein increases when changing the risk category from a lower to a higher category.
# The threshold follows common biological conventions for fold change interpretation:
#     Fold Change > 1.5 → Moderate increase.
#     Fold Change > 2.0 → Strong increase.
# presence_threshold is separately used for presence detection, meaning whether the abundance is considered meaningfully greater than a small value (e.g., to ignore near-zero noise)
# This convention appears for example in proteomics and transcriptomics studies (Smyth, 2004; Ritchie et al., 2015 [limma package]).
# Reference:
# 
# Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. https://doi.org/10.1093/nar/gkv007

# In[88]:


def classify_abundance_patterns(eccontri_df):
    """
    Comprehensive classification of abundance patterns with biological descriptors,
    fold changes, and other key metrics for prioritization.

    Parameters:
    df: DataFrame with columns for Genus, protein_name, Category, norm_abund_contri, Frequency

    Returns:
    DataFrame with pattern classifications and related metrics
    """
    # Make a copy to avoid modifying the original
    result_df = eccontri_df.copy()
    # Define thresholds
    epsilon = 1e-10  # For comparing means (pattern detection)
    presence_threshold = 1e-10  # For defining presence/absence # Minimum abundance to consider as 'present'
    result_df['Frequency'] = result_df['Frequency'].fillna(0).clip(lower=0).astype('int8') 
    freq_mapping = result_df.groupby(['Sites', 'Genus', 'protein_name'], observed=True)[['Frequency', 'idx']].first().reset_index()

    # Calculate mean abundance per category for each genus-protein pair
    mean_by_cat = result_df.groupby(['Sites','Genus', 'protein_name', 'Category'], observed=True)['norm_abund_contri'].mean().reset_index()
    # Pivot to get means by category
    pattern_data = mean_by_cat.pivot_table(
        index=['Sites', 'Genus', 'protein_name'],
        columns='Category',
        values='norm_abund_contri', observed=True
    ).reset_index()
    # Rename the columns for clarity
    pattern_data.columns.name = None
    if 1 in pattern_data.columns:
        pattern_data.rename(columns={1: 'mean_cat1'}, inplace=True)
    else:
        pattern_data['mean_cat1'] = 0

    if 2 in pattern_data.columns:
        pattern_data.rename(columns={2: 'mean_cat2'}, inplace=True)
    else:
        pattern_data['mean_cat2'] = 0

    if 3 in pattern_data.columns:
        pattern_data.rename(columns={3: 'mean_cat3'}, inplace=True)
    else:
        pattern_data['mean_cat3'] = 0

    # Convert categorical columns to regular types before fillna
    for col in pattern_data.columns:
        if isinstance(pattern_data[col].dtype, pd.CategoricalDtype):
            pattern_data[col] = pattern_data[col].astype('float16' if col.startswith('mean_cat') else 'string')

    pattern_data = pattern_data.fillna(0)
    for cat in [1, 2, 3]:
        if cat not in pattern_data.columns:
            pattern_data[cat] = 0

    # Merge pattern data with frequency data
    pattern_data = pd.merge(pattern_data, freq_mapping, on=['Sites', 'Genus', 'protein_name'], how='left')
    # fill nans and Reindex pattern_data columns after merge
    pattern_data = pattern_data.fillna(0)
    pattern_data = pattern_data[['Sites', 'idx', 'Genus', 'protein_name', 'mean_cat1', 'mean_cat2', 'mean_cat3', 'Frequency']]

    # Define function to determine which categories are present
    def get_category_string(row):
        categories = []
        if row['mean_cat1'] > presence_threshold:
            categories.append('1')
        if row['mean_cat2']  > presence_threshold:
            categories.append('2')
        if row['mean_cat3']  > presence_threshold:
            categories.append('3')
        return ''.join(categories)

    pattern_data['category_str'] = pattern_data.apply(get_category_string, axis=1)

    # Tolerance for floating point comparisons
    epsilon = 1e-5

    # Calculate pattern using descriptive naming
    def determine_pattern(row):
        category_str = row['category_str']
        mean1 = row['mean_cat1']
        mean2 = row['mean_cat2']
        mean3 = row['mean_cat3']

        pattern = "other_pattern"

        if category_str == '123':
            # When present in all categories
            if mean1 < mean2 and mean2 < mean3:
                pattern = "steadily_increasing"
            elif mean1 > mean2 and mean2 > mean3:
                pattern = "steadily_decreasing"
            elif mean1 < mean2 and mean2 > mean3:
                if mean1 < mean3:
                    pattern = "peak_at_transition_high"  # Peak at cat2 but ends higher
                elif abs(mean1 - mean3) < epsilon:
                    pattern = "peak_at_transition_equal"  # Perfect peak at cat2
                else:
                    pattern = "peak_at_transition_low"  # Peak at cat2 but ends lower
            elif mean1 > mean2 and mean2 < mean3:
                if mean1 > mean3:
                    pattern = "stress_recovery_partial"  # Valley at cat2 but ends lower
                elif abs(mean1 - mean3) < epsilon:
                    pattern = "stress_recovery_equal"  # Perfect valley at cat2
                else:
                    pattern = "stress_recovery_full"  # Valley at cat2 but ends higher
            elif abs(mean1 - mean2) < epsilon and mean2 < mean3:
                pattern = "late_response"  # Plateaus then increases
            elif abs(mean1 - mean2) < epsilon and abs(mean2 - mean3) < epsilon:
                pattern = "risk_independent"  # Consistent across all
            elif abs(mean1 - mean2) < epsilon and mean2 > mean3:
                pattern = "late_decline"  # Plateaus then decreases

        elif category_str == '12':
            # When only categories 1 and 2 are present
            if mean1 < mean2:
                pattern = "early_increase_only"
            elif mean1 > mean2:
                pattern = "early_decrease_only"
            else:
                pattern = "early_constant"

        elif category_str == '13':
            # When only categories 1 and 3 are present
            if mean1 < mean3:
                pattern = "extreme_response"  # Responds to normal and severe but not transition
            elif mean1 > mean3:
                pattern = "extreme_inhibition"
            else:
                pattern = "extreme_constant"

        elif category_str == '23':
            # When only categories 2 and 3 are present
            if mean2 < mean3:
                pattern = "late_emergence_increasing"  # Absent in normal, increases with risk
            elif mean2 > mean3:
                pattern = "late_emergence_decreasing"
            else:
                pattern = "late_emergence_constant"

        elif category_str == '1':
            pattern = "normal_exclusive"  # Only in normal condition
        elif category_str == '2':
            pattern = "transition_exclusive"  # Only in transition condition
        elif category_str == '3':
            pattern = "severe_exclusive"  # Only in severe condition

        return pattern

    pattern_data['descriptive_pattern'] = pattern_data.apply(determine_pattern, axis=1)

    # Add broader pattern categories for filtering
    pattern_categories = {
        # Increasing patterns
        'steadily_increasing': 'increasing',
        'peak_at_transition_high': 'increasing',
        'peak_at_transition_equal': 'increasing',
        'late_response': 'increasing',
        'stress_recovery_full': 'increasing',
        'early_increase_only': 'increasing',
        'extreme_response': 'increasing',
        'late_emergence_increasing': 'increasing',
        'transition_exclusive': 'increasing',
        'severe_exclusive': 'increasing',

        # Decreasing patterns
        'steadily_decreasing': 'decreasing',
        'peak_at_transition_low': 'decreasing',
        'stress_recovery_partial': 'decreasing',
        'stress_recovery_equal': 'decreasing',
        'late_decline': 'decreasing',
        'early_decrease_only': 'decreasing',
        'extreme_inhibition': 'decreasing',
        'late_emergence_decreasing': 'decreasing',
        'normal_exclusive': 'decreasing',

        # Other patterns
        'risk_independent': 'constant',
        'early_constant': 'constant',
        'extreme_constant': 'constant',
        'late_emergence_constant': 'constant',
        'other_pattern': 'other'
    }

    pattern_data['pattern_category'] = pattern_data['descriptive_pattern'].map(pattern_categories)

    # Calculate fold changes with pseudocount to avoid division by zero
    pseudocount = 1e-5

    # Cat2 vs Cat1 fold change
    pattern_data['fold_change_2vs1'] = (pattern_data['mean_cat2'] + pseudocount) / (pattern_data['mean_cat1'] + pseudocount)
    pattern_data['log2fc_2vs1'] = np.log2(pattern_data['fold_change_2vs1'])

    # Cat3 vs Cat2 fold change
    pattern_data['fold_change_3vs2'] = (pattern_data['mean_cat3'] + pseudocount) / (pattern_data['mean_cat2'] + pseudocount)
    pattern_data['log2fc_3vs2'] = np.log2(pattern_data['fold_change_3vs2'])

    # Cat3 vs Cat1 fold change
    pattern_data['fold_change_3vs1'] = (pattern_data['mean_cat3'] + pseudocount) / (pattern_data['mean_cat1'] + pseudocount)
    pattern_data['log2fc_3vs1'] = np.log2(pattern_data['fold_change_3vs1'])

    # Cap extreme log2fc values
    for col in ['log2fc_2vs1', 'log2fc_3vs2', 'log2fc_3vs1']:
        pattern_data[col] = pattern_data[col].clip(-10, 10)

    # Calculate maximum log2 fold change (absolute value)
    pattern_data['max_abs_log2fc'] = pattern_data[['log2fc_2vs1', 'log2fc_3vs2', 'log2fc_3vs1']].abs().max(axis=1)

    def run_permutation_genus_only(df, n_perm=1000):
        df = df.copy()

        def compute_score(df, category):
            mask = df['Category'] == category
            if mask.sum() == 0:
                return 0
            ## proportion of non-zero abundances in a category
            pres = (df.loc[mask, 'norm_abund_contri'] > 0).mean()
            ## ratio of mean abundance in category vs outside category
            abun = df.loc[mask, 'norm_abund_contri'].mean() / (df.loc[~mask, 'norm_abund_contri'].mean() + 1e-10)
            # degree of enhancement on this category, log to include extrems
            return np.log1p(pres * abun) 

        # First calculate genus-level p-values
        genus_results = {}
        for genus in df['Genus'].unique():
            genus_data = df[df['Genus'] == genus]

            if genus_data['Category'].nunique() < 1:
                continue

            real_score = max(compute_score(genus_data, c) for c in [1, 2, 3])

            cat_presence = set(genus_data['Category'].unique())
            # Accept only if cat 2 and cat 3 also pass
            if cat_presence == {2}:
                pval = 0.003
            elif cat_presence == {3}:
                pval = 0.001
            elif cat_presence == {2, 3}:
                pval = 0.0002
                print(f"Genus {genus} has categories {cat_presence} - assigning p-value 0.001")
            else:          
                # Permutation test: Randomly shuffles category labels 1000 times, for each finds best score
                # and counts how many times the max score agrees with data
                count = 0
                categories = genus_data['Category'].values
                for _ in range(n_perm):  # Fixed syntax
                    perm_data = genus_data.copy()
                    perm_data['Category'] = np.random.permutation(categories)
                    perm_score = max(compute_score(perm_data, c) for c in [1, 2, 3])
                    if perm_score >= real_score:
                        count += 1

                pval = (count + 1) / (n_perm + 1)
            genus_results[genus] = {'real_score': real_score, 'p_value': pval}

        # Now create results with idx for joining
        results = []
        for _, row in df.iterrows():
            genus = row['Genus']
            if genus in genus_results:
                results.append({
                    'idx': row['idx'],
                    'Genus': genus,
                    'real_score': genus_results[genus]['real_score'],
                    'p_value': genus_results[genus]['p_value']
                })

        return pd.DataFrame(results)

    pval_df = run_permutation_genus_only(eccontri_df, n_perm=1000)

    pval_df = pval_df.drop_duplicates(subset=['Genus', 'idx'])
    pattern_data = pattern_data.merge(pval_df, on=['Genus', 'idx'], how='left')

    # Count distinct patterns
    pattern_counts = pattern_data['descriptive_pattern'].value_counts()
    print("Detailed pattern distribution:")
    print(pattern_counts)

    # Count by category
    category_counts = pattern_data['pattern_category'].value_counts()
    print("\nPattern category distribution:")
    print(category_counts)

    # Calculate a combined significance score
    def calculate_significance(row):
        # Pattern category base score
        pattern_score = {
                'increasing': 3.0,
                'decreasing': 1.0,
                'constant': 0.5
            }.get(row['pattern_category'], 0.0)

        # Descriptive pattern bonus
        unique_bonus = {
                'severe_exclusive': 2.0,
                'transition_exclusive': 1.0,
                'normal_exclusive': 0.0
            }.get(row['descriptive_pattern'], 0.0)

        # Fold change, capped
        fc_score = min(abs(row['log2fc_3vs1']), 3.0)
        pattern_data['Frequency'] = pattern_data['Frequency'].fillna(0).clip(lower=0).astype('int8') #changed
        # Frequency (log-scaled)
        freq_score = np.log2(row['Frequency'] + 1)

        # p-value (inverse-log scaled)
        p_value_score = -np.log10(np.clip(row['p_value'], 1e-300, 1))

        # Final score (weighted)
        return (pattern_score * 0.1) + (unique_bonus * 0.2) +  (fc_score * 0.25) + (freq_score * 0.1) + \
                (p_value_score * 0.35)

    pattern_data['significance_score'] = pattern_data.apply(calculate_significance, axis=1)
    # Return the pattern data
    return pattern_data


# ## 4.2. Integrating pattern data with ec_records columns
# The merging_pattern_record function serves as a critical integration step that combines statistical pattern analysis with rich biological metadata. This function takes the pattern analysis results (containing information about abundance trends across risk categories) and merges them with extensive protein metadata from the original enriched dataset.
# The function first ensures data consistency by sorting both dataframes by index, then performs a left join that preserves all pattern data while incorporating valuable biological context from the original dataset. This includes enzyme names and classes, metabolic pathways, functional hierarchies, metal interactions, corrosion mechanisms, and relevance scores.
# This integration step is essential as it bridges statistical observations with biological meaning, allowing subsequent analyses to interpret abundance patterns in the context of protein function, metal interactions, and corrosion relevance. The comprehensive metadata incorporated here was derived from multiple authoritative databases processed in section 9, providing a rich foundation for biological interpretation of the statistical patterns.

# In[89]:


def merging_pattern_record(pattern_df, eccontri_df):
    """
    Merges pattern results with eccontri_df while preserving the two-step merging
    process for important domain-specific reasons. Removes redundant individual
    prevalence columns, keeping only the unified prevalence.

    Parameters:    pattern_result : pandas DataFrame result from perform_abundance_analysis
    eccontri_df : pandas DataFrame original dataframe with protein metadata

    Returns:  integrated_results : pandas DataFrame, clean merged results with pattern data and metadata
    """
    # Sort both df
    eccontri_df = eccontri_df.copy()
    eccontri_df = eccontri_df.sort_values(by=["idx"])
    pattern_df = pattern_df.copy()
    pattern_df = pattern_df.sort_values(by=["idx"])

    # merge for metadata 'idx', 'Sites', 'Genus', 'norm_abund_contri', 'protein_name',

    integrated_results = pd.merge(pattern_df, eccontri_df[['idx', 'Sites', 'Genus',  
                    'protein_name', 'pathway_primary', 'react_primary', 'enzyme_class',
                    'functional_sub', 'functional_child', 'Frequency', 
                    'synergy_child', 'synergy_sub', 'synergy_description',
                    'consolidated_metals', 'hierarchy', 'mechanisms_sub',
                    'mechanisms_child', 'operational_sub', 'operational_combi',
                    'overall_functional_score', 'overall_metal_score',
                    'overall_synergy_score', 'corrosion_relevance_score']],            
                 on=['idx', 'Sites', 'Genus', 'protein_name', 'Frequency'], how='left'
    )

    integrated_results = integrated_results.rename(columns= {"pathway_primary": "pathways", "react_primary": "reactions",
                                                             "synergy_child": "synergy_child_list", "synergy_sub": "synergy_sub_list"})
    for col in ['synergy_child_list', 'synergy_sub_list']:
        if col in integrated_results.columns:
            integrated_results[col] = integrated_results[col].apply(
                lambda x: (
                    x if isinstance(x, np.ndarray)
                    else (np.array(x, dtype=object) if isinstance(x, (list, tuple)) and len(x) > 0
                    else np.array([], dtype=object))
                )
            )
    integrated_results = integrated_results.sort_values(by=["idx"]).reset_index(drop=True)

    return integrated_results


# ## 4.3. Classifying Housekeeping, Niche and Mixed Protein depending on pathways, mechanisms and hierarchy
# The classify_pathways_by_specificity function implements a systematic approach to categorize protein-associated pathways based on their ecological and metabolic roles. The function distinguishes between universal pathways (found across most microorganisms), niche-specific pathways (specialized for particular environments), and mixed pathways (containing elements of both).
# The classification relies on a comprehensive dictionary of universal pathways organized into six major functional categories: energy production, carbon storage, DNA/RNA/protein processes, membrane transport, stress response, and biomolecule synthesis. Each category contains specific pathways that represent core metabolic functions essential for microbial life across diverse environments.For each protein entry, the function analyzes pathway text using regex pattern matching to identify universal pathway components. Then extracts remaining text as potential niche-specific components, assigns a classification based on the presence of universal and/or niche-specific elements
# Records specific universal pathways detected and remaining niche-specific content
# 
# The function also generates summary statistics showing the distribution of pathway classifications across the dataset and identifies the most common universal pathways. This classification is crucial for downstream analysis as it helps distinguish proteins involved in basic cellular maintenance from those potentially specialized for corrosion environments.This classification step provides important context for understanding which proteins represent housekeeping functions and which might be more directly involved in corrosion-specific adaptations, allowing for more targeted analysis of potentially relevant biomarkers.

# In[90]:


def classify_pathways_by_specificity(integrated_results):
    """
    Classify pathways as 'universal', 'niche-specific', or 'mixed' based on matching to universal pathways list
    Parameters: integrated_results : DataFrame with pathway information
    Returns:    DataFrame with additional columns for pathway classification
    """
    # Create a copy to avoid modifying the original
    results = integrated_results.copy(deep=False)

    # Define the universal pathways based on the document
    universal_pathways = {
        "metabolism": "General Metabolism",
        "central metabolism": "1. Energy Production",
        "energy production": "1. Energy Production",
        # Energy Production
        "glycolysis": "1.1. Glycolysis",
        "gluconeogenesis": "1.2. Gluconeogenesis",
        "pentose phosphate": "1.3. Pentose Phosphate Pathway",
        "tca cycle": "1.4. Krebs/TCA Cycle",
        "krebs cycle": "1.4. Krebs/TCA Cycle",
        "citric acid cycle": "1.4. Krebs/TCA Cycle",
        "carbon cycle": "1.4. Krebs/TCA Cycle",
        "electron transport chain": "1.5. Electron Transport Chain",
        "etc": "1.5. Electron Transport Chain",
        "respiration": "1.5. Electron Transport Chain",
        "fermentation": "1.6. Fermentation",
        "atp synthase": "1.7. ATP Synthase",
        "atp synthesis": "1.7. ATP Synthase",
        # Carbon Storage
        "fatty acid synthesis": "2.1. Fatty Acid Synthesis",
        "fatty acid oxidation": "2.2. Fatty Acid Oxidation",
        "beta-oxidation": "2.2. Fatty Acid Oxidation",
        "amino acid degradation": "2.3. Amino Acid Degradation",
        "glycogenesis": "2.4. Carbohydrate Storage",
        "glycogenolysis": "2.4. Carbohydrate Storage",
        "triacylglycerol synthesis": "2.5. Triacylglycerol Synthesis",

        # DNA/RNA/Protein
        "dna replication": "3.1. DNA Replication",
        "dna polymerase": "3.1. DNA Replication",
        "helicase": "3.1. DNA Replication",
        "transcription": "3.2. Transcription",
        "rna polymerase": "3.2. Transcription",
        "translation": "3.3. Translation",
        "ribosomal": "3.3. Translation",
        "ribosome": "3.3. Translation",
        "ribosomes": "3.3. Translation",
        "protein synthesis": "3.3. Translation",
        "protein folding": "3.4. Protein Folding and Chaperones",
        "chaperone": "3.4. Protein Folding and Chaperones",
        "proteasome": "3.5. Proteasome System",
        "protease": "3.5. Proteasome System",
        "amino acid biosynthesis": "3.6. Amino Acid Biosynthesis",

        # Membrane Transport
        "abc transporter": "4.1. ABC Transporters",
        "pilus": "4.2. Pilus and Flagella Formation",
        "flagella": "4.2. Pilus and Flagella Formation",
        "peptidoglycan": "4.3. Cell Wall Maintenance",
        "s-layer": "4.3. Cell Wall Maintenance",
        "lipid membrane": "4.4. Lipid Membrane Synthesis",
        "glycerophospholipid": "4.4. Lipid Membrane Synthesis",

        # Stress Response
        "oxidative stress": "5.1. Oxidative Stress Response",
        "superoxide dismutase": "5.1. Oxidative Stress Response",
        "catalase": "5.1. Oxidative Stress Response",
        "peroxidase": "5.1. Oxidative Stress Response",
        "heat shock": "5.2. Heat Shock Proteins",
        "cold shock": "5.3. Cold Shock Proteins",
        # Biomolecule Synthesis
        "nucleotide biosynthesis": "6.2. Nucleotide Biosynthesis",
    }
    # Add key pathway identifiers that should be prioritized as universals
    core_universals = {
            "glycolysis": "1.1. Glycolysis",
            "Pyruvate metabolism": "1.1. Glycolysis",
            "gluconeogenesis": "1.2. Gluconeogenesis",
            "pentose phosphate": "1.3. Pentose Phosphate Pathway",
            "tca": "1.4. Krebs/TCA Cycle",
            "krebs": "1.4. Krebs/TCA Cycle",
            "citric acid": "1.4. Krebs/TCA Cycle",
            "electron transport": "1.5. Electron Transport Chain",
            "oxidative phosphorylation": "1.5. Electron Transport Chain",
            "dna replication": "3.1. DNA Replication",
            "Chromosome and associated proteins": "3.1. DNA Replication",
            "transcription": "3.2. Transcription", 
            "Messenger RNA biogenesis": "3.2. Transcription",
            "Transfer RNA biogenesis": "3.2. Transcription", 
            "translation": "3.3. Translation"
        }
    # Define mappings from KEGG KO IDs to universal pathway categories
    kegg_universal_map = {
        # Central Carbon Metabolism
        "ko00010": "1.1. Glycolysis",                # Glycolysis / Gluconeogenesis
        "ko00030": "1.3. Pentose Phosphate Pathway", # Pentose phosphate pathway
        "ko00020": "1.4. Krebs/TCA Cycle",           # Citrate cycle (TCA cycle)
        "ko00190": "1.5. Electron Transport Chain",  # Oxidative phosphorylation

        # Amino Acid Metabolism
        "ko00250": "6.1. Amino Acid Biosynthesis",   # Alanine, aspartate and glutamate metabolism
        "ko00260": "6.1. Amino Acid Biosynthesis",   # Glycine, serine and threonine metabolism
        "ko00270": "6.1. Amino Acid Biosynthesis",   # Cysteine and methionine metabolism
        "ko00280": "6.1. Amino Acid Biosynthesis",   # Valine, leucine and isoleucine degradation
        "ko00290": "6.1. Amino Acid Biosynthesis",   # Valine, leucine and isoleucine biosynthesis
        "ko00520": "4.3. Cell Wall Maintenance",  # Amino sugar and nucleotide sugar metabolism

        # Lipid Metabolism
        "ko00061": "2.1. Fatty Acid Synthesis",      # Fatty acid biosynthesis
        "ko00071": "2.2. Fatty Acid Oxidation",      # Fatty acid degradation

        # Nucleotide Metabolism
        "ko00230": "6.2. Nucleotide Biosynthesis",   # Purine metabolism
        "ko00240": "6.2. Nucleotide Biosynthesis",   # Pyrimidine metabolism

        # Replication and Repair
        "ko03030": "3.1. DNA Replication",           # DNA replicatio       
        # Transcription/Translation
        "ko03010": "3.3. Translation",               # Ribosome
        "ko00970": "3.3. Translation",
        "ko03020": "3.2. Transcription",             # RNA polymerase
        "ko00970": "3.2. Transcription",
        "ko03060": "3.4. Protein Folding and Chaperones", # Protein export
    }
    # Add specific amino acid biosynthesis pathways
    amino_acids = ["isoleucine", "valine", "leucine", "alanine", "arginine",  "asparagine", "aspartate", "cysteine", "glutamate", "glutamine",
                  "glycine", "histidine", "lysine", "methionine", "phenylalanine",  "proline", "serine", "threonine", "tryptophan", "tyrosine"]

    # Define terms that should always be considered niche-specific even if they contain universal keywords
    niche_specific_terms = {
        "sos dna repair", "dna repair", "dna repair and recombination",  "flavonoid biosynthesis", "antibiotic biosynthesis", "antimicrobial resistance",
        "viral proteins"}

    for aa in amino_acids:
        universal_pathways[f"{aa} biosynthesis"] = "6.1. Amino Acid Biosynthesis"

    # Convert universal pathways to lowercase
    universal_pathways = {k.lower(): v for k, v in universal_pathways.items()}
    core_universals = {k.lower(): v for k, v in core_universals.items()}
    niche_specific_terms = {term.lower() for term in niche_specific_terms}

    # Then merge for the full matching
    all_universals = universal_pathways.copy()
    all_universals.update(core_universals)

    #===============

    # Initialize classification columns
    if 'idx' in results.columns:
        results["idx"] = results["idx"].astype('int32')
    results["pathway_classification"] = "niche-specific"  # Default
    results["universal_pathways"] = ""
    results["niche_specific_pathways"] = ""

    # Function to classify a pathway
    def classify_pathway(pathway_text):
        if pd.isna(pathway_text) or pathway_text == "":
            return "unknown", "", ""

        pathway_text_lower = pathway_text.lower()

        # Check if any niche-specific terms are present - these override universal classification
        for term in niche_specific_terms:
            if term in pathway_text_lower:
                return "niche-specific", "", pathway_text

        matched_universals = set()
        matched_core = False

        # 1. Handle KEGG IDs first with strong priority
        kegg_ids = re.findall(r'\[br:ko(\d+)\]', pathway_text_lower)
        for ko_id in kegg_ids:
            ko_key = f"ko{ko_id}"
            if ko_key in kegg_universal_map:
                matched_universals.add(kegg_universal_map[ko_key])
                matched_core = True  # KEGG matches in the map are considered core

        # 2. Check for exact or phrase matches with core universals
        for key, value in core_universals.items():
            if re.search(r'\b' + re.escape(key) + r'\b', pathway_text_lower):
                matched_universals.add(value)
                matched_core = True

        # 3. If no core match, check for broader universal matches
        if not matched_core:
            for key, value in all_universals.items():
                if re.search(r'\b' + re.escape(key) + r'\b', pathway_text_lower):
                    matched_universals.add(value)

        # 4. Identify niche-specific content by removing universal terms
        niche_text = pathway_text_lower
        for uni_term in [k for k, v in all_universals.items() if v in matched_universals]:
            niche_text = re.sub(r'\b' + re.escape(uni_term) + r'\b', '', niche_text)

        # Clean up niche text
        niche_text = re.sub(r'\[br:ko\d+\]', '', niche_text)  # Remove KEGG IDs
        niche_text = re.sub(r'enzymes with ec numbers', '', niche_text)  # Remove problematic phrase
        niche_text = re.sub(r'\s+', ' ', niche_text).strip()
        niche_text = re.sub(r'[,;]\s*[,;]', ',', niche_text)
        niche_text = re.sub(r'^[,;]\s*|\s*[,;]$', '', niche_text)

        # Split into terms for more accurate classification
        niche_terms = [term for term in niche_text.split() if len(term) > 3]

        # 5. Classify based on matched patterns
        if matched_universals:
            if matched_core or len(niche_terms) < 5:  # Core match or few niche terms
                classification = "universal"
            else:
                classification = "mixed"
        elif niche_text.strip():
            classification = "niche-specific"
        else:
            classification = "unknown"

        return classification, ", ".join(sorted(list(matched_universals))), niche_text.strip()

    # Apply classification to pathways column if it exists
    if 'pathways' in results.columns:
        # Ensure pathways are strings
        results['pathways'] = results['pathways'].astype(str).fillna('')

        classifications = results['pathways'].apply(classify_pathway)
        results["pathway_classification"] = classifications.apply(lambda x: x[0])
        results["universal_pathways"] = classifications.apply(lambda x: x[1])
        results["niche_specific_pathways"] = classifications.apply(lambda x: x[2])

        # Count each classification
        class_counts = results['pathway_classification'].value_counts()
        print("Pathway classification results:")
        for cls, count in class_counts.items():
            print(f"  - {cls}: {count} ({count/len(results):.1%})")

    return results

def generate_specificity_report(classified_results, output_file=None):
    """
    Generate a report of pathway classifications with examples

    Parameters:  classified_results : DataFrame with pathway classification
                 output_file : str, optional path to save Excel report
    Returns:     Dict with summary statistics
    """
    #def generate_specificity_report(classified_results, output_file=None):
    """
    Generate a summary report of pathway classifications with useful statistics
    without creating a separate DataFrame that caused the array length error.

    Parameters:
        classified_results : DataFrame with pathway classification
        output_file : str, optional path to save summary (not used in this simplified version)

    Returns:
        classified_results : The same DataFrame with summary printed to screen
    """
    # Calculate basic statistics
    total_pathways = len(classified_results)

    # Count classifications
    classification_counts = classified_results['pathway_classification'].value_counts().to_dict()

    # Calculate percentages for each classification
    percentages = {
        cls: (count / total_pathways * 100)
        for cls, count in classification_counts.items()
    }

    # Get universal pathway counts if available
    universal_pathway_counts = {}
    if 'universal_pathways' in classified_results.columns:
        all_universal = []
        for pathways in classified_results['universal_pathways'].dropna():
            if pathways:
                all_universal.extend([p.strip() for p in pathways.split(',')])

        from collections import Counter
        universal_pathway_counts = dict(Counter(all_universal))

    # Get examples of niche-specific pathways
    niche_examples = []
    if 'niche_specific_pathways' in classified_results.columns:
        niche_examples = classified_results.loc[
            classified_results['pathway_classification'] == 'niche-specific',
            'niche_specific_pathways'
        ].dropna().unique()[:20]  # Top 20 examples

    # Print comprehensive summary to screen
    print("\n=== PATHWAY SPECIFICITY ANALYSIS ===")
    print(f"Total pathways analyzed: {total_pathways}")

    print("\nClassification distribution:")
    for cls, count in classification_counts.items():
        print(f"  - {cls}: {count} ({percentages[cls]:.1f}%)")

    # Get all universal pathway entries
    universal_entries = classified_results[classified_results['pathway_classification'].isin(['universal', 'mixed'])]['universal_pathways'].dropna()

    # Split multi-pathway entries and count individual pathways
    all_pathways = []
    for entry in universal_entries:
        pathways = [p.strip() for p in entry.split(', ')]
        all_pathways.extend(pathways)

    # Count occurrences of each unique pathway
    from collections import Counter
    pathway_counts = Counter(all_pathways)

    # Print sorted results
    print("Universal pathway distribution (without duplication):")
    for pathway, count in sorted(pathway_counts.items(), key=lambda x: x[1], reverse=True):
        if pathway:  # Skip empty entries
            print(f"  - {pathway}: {count} occurrences")

    # Return summary statistics as a dictionary as documented
    summary = {
        'total_pathways': total_pathways,
        'classification_counts': classification_counts,
        'percentages': percentages,
        'universal_pathway_counts': universal_pathway_counts,
        'niche_examples': list(niche_examples)
    }

    # Return the original DataFrame 
    return classified_results


# ## 4.4. Patterns Separtation
# The separate_by_pattern function organizes protein-genus pairs into three categories based on their abundance patterns across risk levels:
# 
# Increasing patterns - Proteins showing higher abundance with increasing corrosion risk, suggesting involvement in corrosion processes.
# Inverse patterns - Proteins decreasing as corrosion risk increases. These may represent protective mechanisms being overwhelmed or protective bacteria that actively inhibit corrosion in healthy systems (Category 1) but decline as systems deteriorate, allowing corrosion-promoting organisms to flourish.
# Constant/Mixed patterns - Proteins showing consistent abundance or complex patterns without clear directional trends.
# 
# The function handles any unaccounted patterns, sorts by significance score, and provides summary statistics on distribution and significance across categories, enabling focused analysis on patterns most relevant to corrosion mechanisms.

# In[91]:


def separate_by_pattern(classified_results):
    """
    Separates classified_results into groups based on their pattern
    of presence across risk categories, and sorts by significance.

    Returns:
      increasing_df : DataFrame with patterns suggesting increasing risk
      inverse_df : DataFrame with patterns suggesting decreasing risk
      constant_df : DataFrame with other patterns
    """
    # Create a copy to avoid rename warning
    results = classified_results.copy()

    # Rename to be able to use the name pattern
    if 'pattern_category' in results.columns and 'pattern' not in results.columns:
        results = results.rename(columns={'pattern_category': 'pattern'})

    # Define patterns associated with increasing risk
    increasing_patterns = ['increasing']

    # Define patterns associated with decreasing risk
    inverse_patterns = ['decreasing']

    # Define patterns that are ambiguous or mixed
    mixed_patterns = ['constant', 'other']

    # Separate based on patterns
    increasing_df = results[results['pattern'].isin(increasing_patterns)].copy()
    inverse_df = results[results['pattern'].isin(inverse_patterns)].copy()
    constant_df = results[results['pattern'].isin(mixed_patterns)].copy()

    # Check if there are any patterns not accounted for
    unaccounted_patterns = set(results['pattern'].unique()) - set(increasing_patterns + inverse_patterns + mixed_patterns)
    if unaccounted_patterns:
        print(f"Warning: Found {len(unaccounted_patterns)} unaccounted patterns: {unaccounted_patterns}")
        # Add rows with unaccounted patterns to constant_df
        unaccounted_df = results[results['pattern'].isin(unaccounted_patterns)]
        constant_df = pd.concat([constant_df, unaccounted_df])

    # Sort each dataframe by significance_score if it exists
    if 'significance_score' in results.columns:
        increasing_df = increasing_df.sort_values('significance_score', ascending=False)
        inverse_df = inverse_df.sort_values('significance_score', ascending=False)
        constant_df = constant_df.sort_values('significance_score', ascending=False)

    print(f"Separated {len(increasing_df)} increasing patterns, {len(inverse_df)} inverse patterns, and {len(constant_df)} mixed/other patterns.")

    # Calculate average significance score for each group if available
    if 'significance_score' in results.columns:
        inc_avg = increasing_df['significance_score'].mean() if len(increasing_df) > 0 else 0
        inv_avg = inverse_df['significance_score'].mean() if len(inverse_df) > 0 else 0
        con_avg = constant_df['significance_score'].mean() if len(constant_df) > 0 else 0
        print(f"Average significance scores: Increasing={inc_avg:.2f}, Inverse={inv_avg:.2f}, Constant={con_avg:.2f}")

    return increasing_df, inverse_df, constant_df


# In[92]:


#increasing_df , inverse_df, constant_df = separate_by_pattern(classified_results)


# ## 4.5. Prioritize Markers
# The prioritize_markers function implements a sophisticated scoring system to rank protein-genus pairs based on their potential relevance to corrosion processes. The function combines statistical pattern significance with rich biological context, using a multi-dimensional approach that evaluates candidates from multiple perspectives. The scoring system incorporates several key components:
# 
# Pattern significance (up to 3 points) - Leverages previously calculated significance scores that consider specificity and frequency of proteins across risk categories. Functional categories (up to 4.5 points) - Evaluates proteins against comprehensive functional categories relevant to corrosion, including iron/sulfur redox, acid production, biofilm formation, and other mechanisms, with additional points for synergistic combinations.Pathways and hierarchy relevance (up to 3 points each) - Uses weighted keyword groups to identify corrosion-relevant pathways and functional hierarchies while avoiding redundant scoring. Metal interaction synergies (up to 2 points) - Identifies key metal combinations known to accelerate corrosion, such as iron-sulfur interactions, metal complexation, and Fe-S clusters. Pathway classification tier (up to 2 points) - Assigns higher scores to niche-specific pathways likely specialized for corrosion environments. Effect size (up to 2 points) - Considers the magnitude of abundance changes across risk categories. Corrosion relevance (up to 3 points) - Incorporates pre-existing relevance scores from metadata.
# 
# The function produces a combined score capped at 10 points, creating a balanced ranking that identifies proteins most likely to play meaningful roles in corrosion processes. Additionally, it generates explanatory columns documenting which factors contributed to each protein's score, enhancing interpretability and providing biological context for the statistical rankings.

# In[93]:


def prioritize_markers(increasing_results):
    """
    Prioritize markers based on statistical and biological significance,
    incorporating Frequency and keeping only specific patterns.
    """
    # Copy to avoid modifying original
    results = increasing_results.copy(deep=True)

    # Initialize combined score
    results['combined_score'] = 0.0
    def minmax_norm(series):
        min_ = series.min()
        max_ = series.max()
        if max_ > min_:
            return (series - min_) / (max_ - min_)
        else:
            return pd.Series(0.0, index=series.index)

   # 1. Pattern component (20%) significance + pValue + frequency
    if 'significance_score' in results.columns:
        results['combined_score'] = minmax_norm(results['significance_score']).fillna(0.0) 
    else:
        results['combined_score'] = 0  # Initialize combined_score

    # Add p-value component (smaller p-value = higher score). Prefer existing p_value_score/frequency_score; only compute when missing
    if 'p_value_score' in results.columns:
        results['p_value_score'] = pd.to_numeric(results['p_value_score'], errors='coerce').fillna(0.0)
    elif 'p_value' in results.columns:
        results['p_value_score'] = -np.log10(np.clip(pd.to_numeric(results['p_value'], errors='coerce').fillna(1.0), 1e-300, 1))
    else:
        results['p_value_score'] = 0.0

    if 'frequency_score' in results.columns:
        results['frequency_score'] = pd.to_numeric(results['frequency_score'], errors='coerce').fillna(0.0)
    elif 'Frequency' in results.columns:
        freq = pd.to_numeric(results['Frequency'], errors='coerce').fillna(0.0)
        results['frequency_score'] = (freq - freq.min()) / (freq.max() - freq.min()) if freq.max() != freq.min() else 0.0
    else:
        results['frequency_score'] = 0.0

    results['combined_score'] += minmax_norm(results['p_value_score'])
    results['combined_score'] += minmax_norm(results['frequency_score'])
    results['combined_score'] =  results['combined_score'] * 0.2

    # 2. DATABASE DERIVED SCORES (already computed upstream) (60%)
    try:
        if 'database_combined_score' not in results.columns:
            results['database_combined_score'] = 0.0
        # Components to consider
        comps = ['overall_functional_score','overall_metal_score',
                'overall_synergy_score','corrosion_relevance_score']
        # ensure numeric
        for c in comps:
            if c not in results.columns: results[c] = 0.0
        results[comps] = results[comps].apply(pd.to_numeric, errors='coerce').fillna(0.0)

        database_score = (
            results['overall_functional_score'] * 0.6 +   # differential importance
            results['overall_metal_score']      * 0.4 +   # metals common, still relevant
            results['overall_synergy_score']    * 0.7 +   # multi-pathway advantage
            results['corrosion_relevance_score']* 0.8     # study focus
        )
        results['database_combined_score'] = database_score
        results['combined_score'] += minmax_norm(results['database_combined_score'])*0.6
    except Exception as e:
        print(f"Warning: Error in database scoring: {e}")

    # 3. PATHWAYS COMPONENT (20%) - Focusing ONLY on 'niche-specific' pathways
    if 'pathway_classification' in results.columns:
        # Initialize the score column to 0.0
        results['niche_specific_score'] = 0.0
        niche_mask = results['pathway_classification'] == 'niche-specific'
        results.loc[niche_mask, 'niche_specific_score'] = 2.0
        # Normalize and incorporate the score into the 'combined_score'
        results['combined_score'] += minmax_norm(results['niche_specific_score']) * 0.2
    # 4. FUNCTIONAL SUBCATEGORY already counted on the db part # check again for synergies with synergy detector
    # --- Apply the synergy detector ---
    results_synergy = results.apply(
        lambda row: synergy_detector._detect_row_synergies(row),
        axis=1
    )

    # Convert list of dicts to DataFrame
    synergy_df = pd.DataFrame(results_synergy.tolist(), index=results.index)

    # Merge synergy results with original DataFrame
    results = pd.concat([results, synergy_df], axis=1)

    # Ensure numeric and normalized score
    results['synergy_combi_score'] = pd.to_numeric(results['synergy_combi_score'], errors='coerce').fillna(0.0)
    results['combined_score'] += minmax_norm(results['synergy_combi_score']) * 0.4

    # Sort by updated score
    results = results.sort_values(by='combined_score', ascending=False).reset_index(drop=True)
    # 5. EFFECT SIZE COMPONENT (10%) magnitude of differential expression
    if 'max_abs_log2fc' in results.columns:
        fc_magnitude = pd.to_numeric(results['max_abs_log2fc'], errors='coerce')
        results['fc_score'] = 0.0
        results.loc[fc_magnitude  >= 4,                       'fc_score'] = 2.0
        results.loc[(fc_magnitude >= 2) & (fc_magnitude < 4), 'fc_score'] = 1.5
        results.loc[(fc_magnitude >= 1) & (fc_magnitude < 2), 'fc_score'] = 1.0
        results.loc[(fc_magnitude > 0)  & (fc_magnitude < 1), 'fc_score'] = 0.5
        results['combined_score'] += minmax_norm(results['fc_score'])*0.1

    # 6. Create explanation for combined score
    def create_explanation(row):
        """Generate human-readable explanation of why marker was prioritized"""
        # Get top scoring components
        score_components = {
            'Pattern significance': row.get('significance_score', 0),
            'P-value': row.get('p_value_score', 0),
            'Frequency': row.get('frequency_score', 0),
            'Database': row.get('database_combined_score', 0),
            'Pathways': row.get('pathways_score', 0),
            'Functional mechanisms': row.get('functional_mechanisms_score', 0),
            'Effect size': row.get('fc_score', 0),
            'niche_specific_score': row.get('niche_specific_score', 0),
            'synergy_combi_score': row.get('synergy_combi_score', 0)
        }

        # Sort by score and get top 3
        top_components = sorted(score_components.items(), key=lambda x: x[1], reverse=True)[:3]

        # Only include components with non-zero scores
        explanation = [f"{name} ({score:.1f})" for name, score in top_components if score > 0]

        # Add functional subcategories if available
        if 'functional_sub' in row and pd.notna(row['functional_sub']) and str(row['functional_sub']).strip():
            explanation.append(f"Functions: {row['functional_sub']}")

        return " | ".join(explanation)

    results['explanation'] = results.apply(create_explanation, axis=1)
    # Progress / distribution (avoid warnings if some rows are NaN by summation elsewhere)
	#s = results['combined_score'].dropna()

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Prioritization complete")
    print(f"Score distribution: Mean={results['combined_score'].mean():.2f}, "
          f"Median={results['combined_score'].median():.2f}, "
          f"Min={results['combined_score'].min():.2f}, "
          f"Max={results['combined_score'].max():.2f}")
    print(results[['combined_score','database_combined_score','significance_score','p_value','Frequency']].describe())

    return results


# ## 4.6. Balancing Genus Representation
# 
# The balance_genus_representation function addresses the challenge of over-representation by dominant genera in the protein dataset. Without appropriate balancing, a few abundant genera can dominate the results, obscuring potentially important contributions from less abundant but functionally significant microorganisms.The function implements a two-tier filtering approach:
# 
# * Abundance Threshold Filtering - Uses genus-specific knee-point thresholds to filter out low-abundance proteins that are likely not metabolically significant. These thresholds identify natural breakpoints in abundance distributions for each genus, ensuring that only proteins present at functionally relevant levels are considered.
# * Per-Genus Selection - After abundance filtering, selects the top proteins per genus based on their combined score, maintaining a consistent representation limit across all genera.By selecting top protein per genus, the aim is to prevent bias toward dominant genera by enforcing representation limits whiles maintaining biological diversity by ensuring all genera contribute to the final dataset. 
# 
# The function also incorporates a fallback mechanism: if the abundance threshold filtering removes all proteins for a particular genus, it defaults to selecting the top proteins by combined score, ensuring all genera remain represented in the final dataset.
# By balancing representation across genera while maintaining a focus on functionally significant proteins, this approach creates a more diverse and biologically meaningful set of protein-genus pairs for downstream analysis and visualization.

# In[94]:


def balance_genus_representation(prioritized_markers, eccontri_df, genus_to_threshold, per_genus_count=10):
    """
    Balance representation by selecting 1. only rows with norm_abund_contri ≥ Knee_Abundance are eligible(biologically significantly drives corrosion). 
    2. top proteins per genus based on combined score to limit dominance and preserves breadth across genera
    and applying abundance thresholds from knee point analysis.
    parameters: prioritized_markers are filtered data markers pairs and it is the core data
                eccontri_df original data need to source the column "norm_abund_contri" for comparison
                genus_to_threshold is a dictionary that holds the knee abundance calculate on section 8.11
    """
    # Initialize empty DataFrame for balanced markers
    balanced_markers = pd.DataFrame()

    prioritized_markers['idx'] = prioritized_markers['idx'].astype('int32')
    eccontri_df['idx'] = eccontri_df['idx'].astype('int32')

    # Filter eccontri_df to include only rows present in prioritized_markers
    relevant_idxs = prioritized_markers['idx'].unique()
    filtered_eccontri = eccontri_df[eccontri_df['idx'].isin(relevant_idxs)]

    # Enrich prioritized_markers with abundance values from eccontri_df
    abundance_data = filtered_eccontri[['idx', 'norm_abund_contri']].drop_duplicates()
    # Merge to add abundance data to prioritized_markers
    enriched_markers = pd.merge(
        prioritized_markers,
        abundance_data,
        on='idx',
        how='left'
    )
    # Add Knee Abundance thresholds using genus_to_threshold dictionary
    enriched_markers['Knee_Abundance'] = enriched_markers['Genus'].map(genus_to_threshold)
    enriched_markers['Knee_Abundance'] = enriched_markers['Knee_Abundance'].astype('float16')
    # Ensure no missing values in Knee_Abundance; raise error if any are found
    if enriched_markers['Knee_Abundance'].isnull().any():
        missing_genera = enriched_markers[enriched_markers['Knee_Abundance'].isnull()]['Genus'].unique()
        raise ValueError(f"The following genera are missing from the threshold dictionary: {missing_genera}")
    # Process each genus individually
    for genus, group in enriched_markers.groupby('Genus'):

        # Filter by Knee Abundance threshold
        filtered_group = group[group['norm_abund_contri'].astype('float16') >= group['Knee_Abundance'].astype('float16')]
           # If filtering removes everything, fall back to top N entries by combined score
        if len(filtered_group) == 0:  # FIX: Indented inside loop
            top_group = group.sort_values(
                ['combined_score', 'p_value'],
                ascending=[False, True],
                kind='mergesort'
            ).head(per_genus_count)
        else:  # FIX: Indented inside loop
            top_group = filtered_group.sort_values(
                ['combined_score', 'p_value'],
                ascending=[False, True],
                kind='mergesort'
            ).head(per_genus_count)
        # Append the selected entries to balanced_markers
        balanced_markers = pd.concat([balanced_markers, top_group])
    # Keep relevant columns
    balanced_markers['Category'] = balanced_markers['Sites'].map(category_dict)

    # 2) Quick verification (diagnostic only)
    _missing = balanced_markers['Category'].isna()
    _missing_count = int(_missing.sum())
    if _missing_count > 0:
        print(f"Warning: {_missing_count} entries in balanced markers have missing Category information.")

    balanced_markers  = balanced_markers.drop(columns=['category_str', 'Frequency', 'descriptive_pattern', 'pattern',  'Knee_Abundance'])

    # Return the final balanced markers sorted by combined score
    return balanced_markers.sort_values('combined_score', ascending=False)


# In[95]:


#balanced_markers = balance_genus_representation(prioritized_markers, ECcontri_Uniprot_Enriched, genus_to_threshold, per_genus_count=10)


# ## 4.7 Creating Marker Groups
# 
# The `create_marker_groups` function organizes filtered protein-genus pairs into specialized subsets that facilitate targeted analysis of different aspects of corrosion processes. Starting with the balanced marker dataset, this function generates multiple biologically meaningful groupings based on statistical patterns, functional characteristics, and corrosion relevance, they are:
# 
# | Group Name                  | Short Description                                                                           |
# |-----------------------------|--------------------------------------------------------------------------------------------|
# | top_markers                 | High-confidence subset of top N markers by combined score                                  |
# | significant_markers         | Statistically significant markers across risk categories                                   |
# | high_metals_relevance       | Markers with high relevance for metal interactions                                         |
# | high_functional_relevance   | Markers with high functional category relevance                                            |
# | high_synergy_relevance      | Markers with top synergy relevance                                                         |
# | high_pathway_relevance      | Markers most relevant to known pathways                                                    |
# | high_niche_relevance         | Markers with high tier/niche specificity                                                   |
# | mechanism_all               | Markers associated with any corrosion mechanism                                            |
# | consolidated_metals         | Markers interacting with metallic species                                                  |
# | pathways_all                | Markers with any annotated metabolic/functional pathway                                    |
# | functional_categories       | Markers assigned to any functional category (child or present)                             |
# | inorganic_acid_complexes    | Markers predicted to form Fe-inorganic acid complexes                                      |
# | organic_acid_complexes      | Markers linked to Fe-organic acid/humic/fulvic complexation                               |
# | biofilm_complexes           | Markers linked to biofilm-mediated complexation                                            |
# | unified_synergy             | Markers with top synergy scores across dimensions                                          |
# | high_biological_relevance   | Markers with highest combined biological relevance                                         |
# | operational_all             | Markers annotated with operational/field environmental factors                            |
# | corrosion_critical          | Markers with highest scores and critical metal/functional relevance (potential drivers)    |

# In[96]:


def create_marker_groups(balanced_markers, top_count=200, threshold_percentile=0.75):
    """ Create specialized marker groups based on various criteria.
    Parameters: balanced_markers: DataFrame with markers balanced by knee_abundance and genus representation
        top_count: Number of top markers to include (default: 200)
        threshold_percentile: Percentile for determining high-relevance markers (default: 0.75)
    Returns:   Dictionary of marker groups by different criteria """
    balanced_markers = balanced_markers.copy(deep=True)

    # Clean None/NaN values in list columns
    list_cols = ['consolidated_metals', 'synergy_child_list', 'synergy_sub_list']
    for col in list_cols:
        if col in balanced_markers.columns:
            def deep_clean_list(x):
                # None/NaN as the list itself
                if x is None or (isinstance(x, float) and pd.isna(x)):
                    return []

                # If it's a list, clean all items inside
                if isinstance(x, list):
                    cleaned = []
                    for item in x:
                        # Skip None/NaN items
                        if item is None:
                            continue
                        try:
                            if pd.isna(item):
                                continue
                        except:
                            pass
                        # Skip string 'nan'
                        if isinstance(item, str) and item.strip().lower() in ('nan', 'none', ''):
                            continue
                        cleaned.append(item)
                    return cleaned

                # Not a list - return as-is (shouldn't happen after preprocessing)
                return x

            balanced_markers[col] = balanced_markers[col].apply(deep_clean_list)


    # Minimal helper regex used to split semicolon/comma lists in strings
    _SPLIT_RE = re.compile(r'[;,]')

    def create_combi(row: dict, sub_col: str, child_col: str) -> Optional[Dict[str, List[str]]]:
        """
        Create a small combi dict {sub: [child, ...]} for a row.
        - Preserves single-string children as single-element lists.
        - Accepts child as: None, str, list/tuple/set, dict (will flatten values).
        - If child is a string containing ';' or ',' it will be split into multiple children.
        - Returns None when no valid child values exist.
        Usage (same signature as the existing helper):
            balanced_markers['functional_combi'] = balanced_markers.apply(
                lambda r: create_combi(r, 'functional_sub', 'functional_child'), axis=1)
        """
        sub = row.get(sub_col)
        if sub is None or pd.isna(sub) or str(sub).strip() in ('<NA>', 'nan', 'None', ''):
            return None

        child = row.get(child_col)
        if child is None or pd.isna(child) or str(child).strip() in ('<NA>', 'nan', 'None', ''):
            return None

        children: List[str] = []

        # list-like children
        if isinstance(child, (list, tuple, set)):
            for x in child:
                if x is None:
                    continue
                sx = str(x).strip()
                if sx:
                    children.append(sx)

        # dict child -> flatten values
        elif isinstance(child, dict):
            for v in child.values():
                if v is None:
                    continue
                if isinstance(v, (list, tuple, set)):
                    for x in v:
                        if x is None:
                            continue
                        sx = str(x).strip()
                        if sx:
                            children.append(sx)
                else:
                    sv = str(v).strip()
                    if sv:
                        children.append(sv)

        # scalar / string child
        else:
            s = str(child).strip()
            if not s:
                return None
            if ';' in s or ',' in s:
                parts = [p.strip() for p in _SPLIT_RE.split(s) if p.strip()]
                children.extend(parts)
            else:
                children.append(s)

        if not children:
            return None

        # Check if sub is actually valid (not string '<NA>')
        sub_str = str(sub).strip()
        if sub_str.lower() in ('<na>', 'nan', 'none', ''):
            return None

        return {sub_str: children}

    # Apply it to mechanism and functional
    if {'mechanisms_sub', 'mechanisms_child'} <= set(balanced_markers.columns):
        balanced_markers['mechanisms_combi'] = balanced_markers.apply(
            lambda row: create_combi(row, 'mechanisms_sub', 'mechanisms_child'), axis=1)
    else:
        balanced_markers['mechanisms_combi'] = None

    if {'functional_sub', 'functional_child'} <= set(balanced_markers.columns):
        balanced_markers['functional_combi'] = balanced_markers.apply(
            lambda row: create_combi(row, 'functional_sub', 'functional_child'), axis=1)
    else:
        balanced_markers['functional_combi'] = None

    # parse the combi columns:
    for col in ['operational_combi', 'functional_combi', 'mechanisms_combi']:
        if col in balanced_markers.columns:                    
            def parse_dict_column(val: Any, col_name: Optional[str] = None) -> Optional[Dict[str, List[str]]]:
                """
                Parse and clean various representations into a dict{sub:[child,...]} or None.
                Accepts:  - dict -> cleaned
                - list/tuple of dicts -> merged
                - string representations (literal_eval of dict/list) -> parsed & cleaned
                - strings with ';' or ',' -> treated as scalar children where appropriate
                Cleans:- removes keys with None or empty children
                - normalizes scalar children to 1-element lists
                - splits semicolon/comma-separated strings into multiple children
                Returns:- dict with non-empty lists OR None when no valid content found.
                Note: col_name is optional and only used for a safe diagnostic message on parse errors.
                """
                if val is None or (isinstance(val, float) and getattr(val, 'nan', False)):
                    return None

                d: Optional[Dict] = None

                # already a dict
                if isinstance(val, dict):
                    d = val

                # merge list/tuple of dicts
                elif isinstance(val, (list, tuple)):
                    merged: Dict[str, Any] = {}
                    for item in val:
                        if isinstance(item, dict):
                            merged.update(item)
                    d = merged if merged else None

                # string: attempt literal_eval, else fail safe
                elif isinstance(val, str):
                    s = val.strip()
                    if s in ('{}', '', 'nan', 'None'):
                        return None
                    try:
                        parsed = ast.literal_eval(s)
                        if isinstance(parsed, dict):
                            d = parsed
                        elif isinstance(parsed, (list, tuple)):
                            merged: Dict[str, Any] = {}
                            for item in parsed:
                                if isinstance(item, dict):
                                    merged.update(item)
                            d = merged if merged else None
                        else:
                            return None
                    except Exception:
                        snippet = (repr(val)[:120] + '...') if val is not None else 'None'
                        key_label = f" for column '{col_name}'" if col_name else ''
                        print(f"Warning: Could not parse{key_label} value: {snippet}")
                        return None

                else:
                    return None

                if not isinstance(d, dict):
                    return None

                cleaned: Dict[str, List[str]] = {}
                for k, v in d.items():
                    if v is None:
                        continue

                    # list-like values
                    if isinstance(v, (list, tuple, set)):
                        kids = [str(x).strip() for x in v if x is not None and str(x).strip()]
                        if kids:
                            cleaned[str(k).strip()] = kids
                        continue

                    # nested dict -> flatten its values
                    if isinstance(v, dict):
                        kids: List[str] = []
                        for vv in v.values():
                            if isinstance(vv, (list, tuple, set)):
                                for x in vv:
                                    if x is None:
                                        continue
                                    sx = str(x).strip()
                                    if sx:
                                        kids.append(sx)
                            elif isinstance(vv, str) and vv.strip():
                                kids.append(vv.strip())
                        if kids:
                            cleaned[str(k).strip()] = kids
                        continue

                    # scalar/string value
                    if isinstance(v, str):
                        s = v.strip()
                        if not s:
                            continue
                        if ';' in s or ',' in s:
                            parts = [p.strip() for p in _SPLIT_RE.split(s) if p.strip()]
                            if parts:
                                cleaned[str(k).strip()] = parts
                        else:
                            cleaned[str(k).strip()] = [s]
                        continue

                    # other scalar -> coerce to string
                    s = str(v).strip()
                    if s:
                        cleaned[str(k).strip()] = [s]

                return cleaned if cleaned else None
            balanced_markers[col] = balanced_markers[col].apply(parse_dict_column)
            valid_count = balanced_markers[col].apply(lambda x: isinstance(x, dict) and len(x) > 0).sum()
            print(f"  {col}: {valid_count} valid dicts")

    # Initialize groups dictionary
    groups = {}

    # Check if combined_score exists, otherwise calculate it if needed
    if 'combined_score' not in balanced_markers.columns:
        print("Warning: combined_score column not found. Creating basic score.")
        if 'significance_score' in balanced_markers.columns:
            balanced_markers['combined_score'] = balanced_markers['significance_score']
        else:
            balanced_markers['combined_score'] = 0.0

    # 1. TOP MARKERS
    top_threshold = balanced_markers['combined_score'].quantile(threshold_percentile)
    groups['top_markers'] = balanced_markers[balanced_markers['combined_score'] >= top_threshold].sort_values('combined_score', ascending=False)\
                            .head(top_count)
    # 2. SIGNIFICANCE GROUPS
    if 'significance_score' in balanced_markers.columns:
        sig_threshold = balanced_markers['significance_score'].quantile(threshold_percentile)
        groups['significant_markers'] = balanced_markers[balanced_markers['significance_score'] >= sig_threshold]

    # 3. COMPONENT SCORE GROUPS
    component_score_fields = {
        'overall_metal_score': 'high_metals_relevance',
        'overall_functional_score': 'high_functional_relevance',
        'overall_synergy_score': 'high_synergy_relevance',
        'corrosion_relevance_score': 'high_corrosion_relevance',
        'niche_specific_score': 'high_niche_relevance',
    }
    for score_field, group_name in component_score_fields.items():
        if score_field in balanced_markers.columns:
            s = pd.to_numeric(balanced_markers[score_field], errors='coerce')
            if not s.isnull().all():
                threshold = s.quantile(threshold_percentile)
                if threshold > 0:
                    groups[group_name] = balanced_markers[s >= threshold]
                else:
                    print(f"Warning: Threshold for {score_field} is non-positive ({threshold}). Skipping group {group_name}.")

    # 4. MECHANISMS - One consolidated group for all mechanisms
    if 'mechanisms_combi' in balanced_markers.columns:
        mask = balanced_markers['mechanisms_combi'].apply(
                    lambda v: isinstance(v, dict) and bool(v) )
        if mask.any():
            groups['mechanisms_all'] = balanced_markers[mask].sort_values('combined_score', ascending=False)

    # 5. CONSOLIDATED METALS - Only markers with MULTIPLE metals or HIGH metal score
    if 'consolidated_metals' in balanced_markers.columns and 'overall_metal_score' in balanced_markers.columns:

        # Require either 2+ metals OR high metal score
        metal_count_mask = balanced_markers['consolidated_metals'].apply(
            lambda v: ((isinstance(v, (list, tuple, np.ndarray)) and 
                not (isinstance(v, np.ndarray) and v.ndim == 0) and 
                len(v) >= 2) if v is not None else False))

        metal_score_threshold = balanced_markers['overall_metal_score'].quantile(0.6)  # Top 40%
        metal_score_mask = balanced_markers['overall_metal_score'] >= metal_score_threshold

        metal_mask = metal_count_mask | metal_score_mask
        if metal_mask.any():
            groups['consolidated_metals'] = balanced_markers[metal_mask].sort_values('combined_score', ascending=False)

    # 6. PATHWAYS GROUP (pathways is expected as a single string; guard against non-strings)
    if 'pathways' in balanced_markers.columns and 'niche_specific_score' in balanced_markers.columns:
        # Require both pathway presence AND high niche specificity
        pathway_threshold = balanced_markers['niche_specific_score'].quantile(0.25)  # Top 75%
        path_mask = balanced_markers['pathways'].apply(lambda v: isinstance(v, str) and v.strip() != '') & \
            (balanced_markers['niche_specific_score'] >= pathway_threshold)
        if path_mask.any():
            groups['pathways_all'] = balanced_markers[path_mask].sort_values('combined_score', ascending=False)
    # 7. FUNCTIONAL CATEGORIES all functional cat are single strings. Create single unified group
    if 'functional_combi' in balanced_markers.columns and 'overall_functional_score' in balanced_markers.columns:
        # Require both functional presence AND high functional score
        func_threshold = balanced_markers['overall_functional_score'].quantile(0.6)  # Top 40%
        def _has_nonempty_child(d):
            return isinstance(d, dict) and any(
                isinstance(v, (list, tuple)) and any(isinstance(x, str) and x.strip() for x in v)
                for v in d.values()
            )
        func_mask = balanced_markers['functional_combi'].apply(_has_nonempty_child) & \
                    (balanced_markers['overall_functional_score'] >= func_threshold)
        if func_mask.any():
            groups['functional_all'] = balanced_markers[func_mask].sort_values('combined_score', ascending=False)

    # 8. COMPLEXATION GROUPS
    # 1) Curated inorganic complexation ligands (ions/oxyanions + common acid names)
    ORGANIC_INDICATORS = ('acetate', 'acetic acid', 'oxalate', 'oxalic acid', 'fatty acid', 'butyric acid','propionate', 'lactate', 'formate', 'citrate', 'succinate', 'fumarate', 'malate', 'pyruvate', 
                            'glycolysis', 'tca cycle', 'hydrocarbon degradation', 'iron oxalate', 'iron acetate','organic acid', 'propionic acid', 'lactic acid', 'formic acid', 'citric acid', 
                            'succinic acid', 'fumaric acid', 'malic acid', 'pyruvic acid', 'carboxylic acid')
    INORGANIC_INDICATORS = ('iron oxide', 'manganese oxide', 'ferrihydrite', 'hematite', 'goethite', 'magnetite', 
                                'siderite', 'carbonates', 'phosphates', 'sulfates', 'nitrates', 'nitrites', 'chlorides', 'fluorides', 'oxides', 'ochre', 'mno2', 'aluminum oxide', 
                                'chromium oxide', 'zinc oxide', 'lead oxide')
    BIOFILM_INDICATORS = ('biofilm', 'adhesion', 'colonization', 'biofilm_formation', 'extracellular_matrix', 'exopolysaccharide', 'eps', 'quorum_sensing', 'biofilm maturation', 'metal_chelation', 
                        'siderophore', 'chelator', 'metallophore', 'metal complexation', 'metal sequestration', 'fe-s cluster', 'iron sulfide', 'bimetallic', 'galvanic couple', 'extracellular polymeric substance', 'exopolymer', 'attachment', 'surface attachment', 
                        'matrix production', 'metal_binding', 'complexation', 'iron chelation', 'metal transport', 'metal solubilization', 'mineral dissolution', 'metal immobilization')
    # all lowercase for matching non-metal children
    METALS_INDICATORS = ('al','as','ba','ca','co','cr','cu','cd','cl-','fe','h+','hg','mg','ni','pb','po4-','so3-','so4-','zn',
        'na+','s','se','v5+','f-','mo','no2-','no3-','s2-','co3-','k+','mn','s2o3-','sr','v'
    )
    INORGANIC_INDICATORS = tuple(set(INORGANIC_INDICATORS) | set(METALS_INDICATORS))
    # normalize inputs once (use global normalize_listlike)
    ORGANIC_INDICATORS_LOWER = {s.lower() for s in ORGANIC_INDICATORS}
    INORGANIC_MINERAL_SUBSTRINGS_LOWER = {s.lower() for s in INORGANIC_INDICATORS}
    BIOFILM_INDICATORS_LOWER = {s.lower() for s in BIOFILM_INDICATORS}

    # lightweight tokeniser that preserves + and - inside tokens
    _token_re = re.compile(r"[A-Za-z0-9\+\-]+")
    def _token_set(text):
        return set(_token_re.findall(str(text).lower()))

    # canonical metal set (lower) from global metal_mapping values
    METAL_CANONICAL_LOWER = {str(v).strip().lower() for v in metal_mapping.values() if v}

    # Classify child terms 
    def classify_complexation_terms(row):
        """
        Classify all child terms into inorganic/organic/biofilm complexation categories.
        Only processes rows with corrosion-relevant metals.
        """        
        # 1. Check metal requirement  build flat metal sets (preserve case for display; lower for matching)
        mets_seq = normalize_listlike(row.get('consolidated_metals'))
        met_set_display = {str(m).strip() for m in mets_seq  if m}

        met_set_lower = {m.lower() for m in met_set_display}
        # exact token universe (lowercased) taken from exploded uniques, excluding None/nan
        CORROSION_METALS_LOWER = {'cr+3', 'cu+',  'fe+2', 'fe+3', 'ni+2', 'mn+2', 'zn+2', 'moo4-2'}
        if not (met_set_lower & CORROSION_METALS_LOWER):
            return '', '', ''
        # 🔍 DIAGNOSTIC: Print when we have a corrosion metal
        corr_metals_found = met_set_lower & CORROSION_METALS_LOWER
        print(f"Row with corrosion metals: {corr_metals_found}")

        # 2. Collect ALL child terms from all child columns
        all_children = []

        # Metal terms (already clean list)
        if 'consolidated_metals' in row:
            all_children.extend(list(met_set_display))

        # Synergy children
        #balanced_markers["synergy_child_list"] =  balanced_markers["synergy_child_list"].apply(lambda x: x if isinstance(x, list) else ([] if pd.isna(x) else [x]))
        scl = row.get('synergy_child_list', [])
        if isinstance(scl, (list, tuple)):
            # Assuming clean_list_elements or prior logic has cleaned it
            all_children.extend([str(it).strip() for it in scl if it is not None and str(it).strip()])

        # Collect from combi dictionaries (now properly parsed)
        for col in ['operational_combi', 'functional_combi', 'mechanisms_combi']:
            combi_dict = row.get(col)
            if isinstance(combi_dict, dict):
                for kids in combi_dict.values():
                    if isinstance(kids, (list, tuple, set)):
                        for k in kids:
                            if k is not None:
                                kk = str(k).strip()
                                if kk and kk.lower() not in ('nan', 'none', '<na>', ''):
                                    all_children.append(kk)
       # Functional child (single string)
        val = row.get('functional_child')
        if val is not None and not (isinstance(val, float) and pd.isna(val)):
            s = str(val).strip()
            if s and s.lower() not in ('nan', 'none', '<na>', ''):
                all_children.append(s)

        # Mechanisms child (single string)
        val = row.get('mechanisms_child')
        if val is not None and not (isinstance(val, float) and pd.isna(val)):
            s = str(val).strip()
            if s and s.lower() not in ('nan', 'none', '<na>', ''):
                all_children.append(s)

        # Operational child doesnt exist,we extract them from operational_combi dictionary
        oc = row.get('operational_combi')
        if isinstance(oc, dict):
            for kids in oc.values():
                if isinstance(kids, (list, tuple, set)):
                    for k in kids:
                        if k is None:
                            continue
                        kk = str(k).strip()
                        if kk:
                            all_children.append(kk)
                elif isinstance(kids, dict):
                    # flatten nested dict values as well
                    for vv in kids.values():
                        if isinstance(vv, (list, tuple, set)):
                            for x in vv:
                                if x is None:
                                    continue
                                xx = str(x).strip()
                                if xx:
                                    all_children.append(xx)
                        elif isinstance(vv, str) and vv.strip():
                            all_children.append(vv.strip())
                elif isinstance(kids, str) and kids.strip():
                    all_children.append(kids.strip())
        # 🔍 DIAGNOSTIC: Print collected children
        print(f"  All children collected: {all_children[:10]}") 

        # 3. Classify each child into the 3 complexation categories
        inorganic_terms = []
        organic_terms = []
        biofilm_terms = []

        # seed explicit metals once (preserve case in output)
        metal_children = list(met_set_display)
        inorganic_terms.extend(metal_children)

        # classify non-metal children; match on lowercased        # child classification (token-aware)
        for child in all_children:
            if child is None:
                continue
            cs = str(child).strip()
            if not cs:
                continue
            child_low = cs.lower()
            toks = _token_set(child_low)

            # 1) exact-metal token match (high confidence)
            if toks & METAL_CANONICAL_LOWER:
                inorganic_terms.append(cs)
                continue

            # 2) inorganic minerals (longer substrings only)
            if any(min_sub in child_low for min_sub in INORGANIC_MINERAL_SUBSTRINGS_LOWER):
                inorganic_terms.append(cs)

            # 3) organic indicators: token overlap or reasonably long substring
            if toks & ORGANIC_INDICATORS_LOWER or any(len(ind) > 3 and ind in child_low for ind in ORGANIC_INDICATORS_LOWER):
                organic_terms.append(cs)

            # 4) biofilm indicators: token overlap or longer substring
            if toks & BIOFILM_INDICATORS_LOWER or any(len(ind) > 3 and ind in child_low for ind in BIOFILM_INDICATORS_LOWER):
                biofilm_terms.append(cs)

        # 4. Return semicolon-separated unique terms
        return (
            '; '.join(sorted(set(inorganic_terms))) if inorganic_terms else '',
            '; '.join(sorted(set(organic_terms))) if organic_terms else '',
            '; '.join(sorted(set(biofilm_terms))) if biofilm_terms else ''
        )
    # Apply to dataframe
    complexation_results = balanced_markers.apply(classify_complexation_terms, axis=1)
    balanced_markers['inorganic_complex'] = complexation_results.apply(lambda x: x[0])
    balanced_markers['organic_complex'] = complexation_results.apply(lambda x: x[1])
    balanced_markers['biofilm_complex'] = complexation_results.apply(lambda x: x[2])

    # Create groups from non-empty complexation columns
    inorganic_markers = balanced_markers[balanced_markers['inorganic_complex'].str.len() > 0]
    if len(inorganic_markers) > 0:
        groups['inorganic_acid_complexes'] = inorganic_markers.sort_values('combined_score', ascending=False)

    organic_markers = balanced_markers[balanced_markers['organic_complex'].str.len() > 0]
    if len(organic_markers) > 0:
        groups['organic_acid_complexes'] = organic_markers.sort_values('combined_score', ascending=False)

    biofilm_markers = balanced_markers[balanced_markers['biofilm_complex'].str.len() > 0]
    if len(biofilm_markers) > 0:
        groups['biofilm_complexes'] = biofilm_markers.sort_values('combined_score', ascending=False)
    # 9. OPERATIONAL ENVIRONMENTAL FACTORS GROUP
    if 'operational_combi' in balanced_markers.columns:
        ope_mask = balanced_markers['operational_combi'].apply(
            lambda x: isinstance(x, dict) and len(x) > 0)
        if ope_mask.any():
            groups['operational_all'] = balanced_markers[ope_mask].sort_values('combined_score', ascending=False)
            print(f"  - operational_all: {len(groups['operational_all'])} markers")
    # 10. Synergy group different from synergy child and synergy sub lists. It was made from priority functional, operational and metal subcategories on the actual filtered data 
    if 'synergy_combi' in balanced_markers.columns and 'overall_synergy_score' in balanced_markers.columns:
        balanced_markers["synergy_combi"] =  balanced_markers["synergy_combi"].apply(
        lambda x: x if isinstance(x, list) else ([] if pd.isna(x) else [x]))
        # Require both synergy presence AND high synergy score
        synergy_threshold = balanced_markers['overall_synergy_score'].quantile(0.5)  # Top 50%

        # Also check that synergy_combi has multiple elements (real synergy, not single category)
        def has_multi_element_synergy(combi):
            if not isinstance(combi, list): # this is a list:[organic_acid_metabolism, Fe, Ni, S, direct_effect]
                return False
            # Filter out empty strings and None
            valid_elements = [x for x in combi if x and str(x).strip() not in ['', 'nan', 'None']]
            return len(valid_elements) >= 2

        synergy_mask = (balanced_markers['synergy_combi'].apply(has_multi_element_synergy)) & \
                    (balanced_markers['overall_synergy_score'] >= synergy_threshold)
        if synergy_mask.any():
            groups['synergy_all'] = balanced_markers[synergy_mask].sort_values('combined_score', ascending=False)
    # 12. HIGH BIOLOGICAL RELEVANCE, based on pathway score
    try:
        bio_score = np.zeros(len(balanced_markers))
        bio_cols = [c for c in ['overall_functional_score', 'overall_synergy_score', 'niche_specific_score'] if c in balanced_markers.columns]
        for col in bio_cols:
            bio_score += pd.to_numeric(balanced_markers[col], errors='coerce').fillna(0)
        if bio_score.sum() > 0:
            bio_threshold = np.percentile(bio_score[bio_score > 0], threshold_percentile * 100)
            groups['high_biological_relevance'] = balanced_markers[bio_score >= bio_threshold]
    except Exception as e:
        print(f"Warning: Could not create biological relevance group: {e}")

    # 13. CRITICAL CORROSION MARKERS -- 
    try:
        score_mask = pd.Series(False, index=balanced_markers.index)
        metal_mask = pd.Series(False, index=balanced_markers.index)
        # Base: top 20% combined score
        if 'consolidated_metals' in balanced_markers.columns:
            score_mask = balanced_markers['combined_score'] >= balanced_markers['combined_score'].quantile(0.8)
            # Metal/ionic species filter# column rows are symbols separated by semicolons ;
            CORROSION_METALS = {"Cr", "Cu", "Fe", "Ni", "Mn", "Zn", "Mo"}
            metal_mask = balanced_markers["consolidated_metals"].apply(
                lambda seq: (
                    bool({metal_mapping.get(str(m).strip().lower(), str(m).strip()) 
                        for m in normalize_listlike(seq)} & CORROSION_METALS)
                    if seq is not None and not (isinstance(seq, float) and pd.isna(seq))
                    else False
                )
            )

        # --- Pathway mask --- single str column not to be striped
        if 'niche_specific_pathways' in balanced_markers.columns:
            path_mask = balanced_markers['niche_specific_pathways'].apply(
                lambda v: isinstance(v, str) and v.strip().lower() not in ('', 'nan', 'none')
            )
        else:
            path_mask = pd.Series(False, index=balanced_markers.index)

        if 'synergy_child_list' in balanced_markers.columns:
            synergy_mask = balanced_markers['synergy_child_list'].apply(
                    lambda v: ((isinstance(v, (list, tuple, np.ndarray)) and 
         not (isinstance(v, np.ndarray) and v.ndim == 0) and 
         len(v) > 0) if v is not None else False)
            )

        else:
            synergy_mask = pd.Series(False, index=balanced_markers.index)

        # --- Combine masks ---
        final_mask = score_mask & ((metal_mask & path_mask) | (path_mask & synergy_mask) | (metal_mask & synergy_mask))
        groups['corrosion_critical' ] = balanced_markers[final_mask].sort_values('combined_score', ascending=False)

    except Exception as e:
        print(f"Error in corrosion_critical: {type(e).__name__}: {e}")

    # 14. Group sizes
    group_sizes = {name: len(df) for name, df in groups.items()}
    print("Marker groups created:")
    for name, size in sorted(group_sizes.items(), key=lambda x: x[1], reverse=True):
        print(f"  - {name}: {size} markers")

    # Print Fe complexation analysis summary
    complexation_groups = ['inorganic_acid_complexes', 'organic_acid_complexes', 'biofilm_complexes']
    present_complexation = [g for g in complexation_groups if g in groups]
    if present_complexation:
        print(f"\nFe Complexation Analysis (based on physicochemical TOC study):")
        for group_name in present_complexation:
            group_df = groups[group_name]
            print(f"  - {group_name}: {len(group_df)} markers")
        if 'biofilm_complexes' in groups and len(present_complexation) > 1:
            biofilm_ids = set(groups['biofilm_complexes'].index)
            for acid_group in ['inorganic_acid_complexes', 'organic_acid_complexes']:
                if acid_group in groups:
                    acid_ids = set(groups[acid_group].index)
                    overlap = len(biofilm_ids.intersection(acid_ids))
                    if overlap > 0:
                        print(f"    * Biofilm-{acid_group.replace('_complexes', '')} overlap: {overlap} markers")

    return groups


# ## 4.8. Generating Marker Report
# The generate_streamlined_report function transforms complex analytical results into an accessible, multi-sheet report that highlights the most relevant protein-genus pairs associated with corrosion processes. The function extracts prioritized markers (using balanced markers if available) and organizes them into four specialized sheets:
# 
# Top Corrosion Markers: A focused presentation of the highest-scoring markers with core information including genus, protein name, enzyme details, corrosion mechanisms, metal involvement, and pattern characteristics.
# Visualization Data: A streamlined dataset optimized for creating abundance pattern visualizations across risk categories.
# Detailed Scores: A comprehensive breakdown of scoring components (prevalence, specificity, frequency, pathway relevance, etc.) to provide transparency into the prioritization process.
# Mechanisms Focus: An expanded view that separates multi-mechanism proteins into individual entries, facilitating mechanism-specific analysis.
# 
# The function also generates a summary sheet with metadata about the analysis, including total markers analyzed, balancing status, and pattern distribution. When provided with an output path, the report is saved as an Excel file with appropriately formatted column names for improved readability. This structured reporting approach enables researchers to quickly identify the most promising corrosion biomarkers while maintaining access to the underlying analytical details.

# In[97]:


def consolidate_analysis_results(analysis_results, marker_groups):
    """
    Consolidate all analysis results and marker groups into a single dictionary
    for easy parquet storage.
    Parameters:   analysis_results: Dictionary containing analysis dataframes
                  marker_groups: Dictionary containing marker groups
    Returns:      consolidated_results: Dictionary with all dataframes in one place
    """
    # Create a new consolidated dictionary
    consolidated_results = {}

    # Add main analysis dataframes
    key_dataframes = [
        'pattern_data',
        'integrated_results',
        'classified_results',
        'increasing_markers',
        'inverse_markers',
        'prioritized_markers',
        'balanced_markers'
    ]

    for key in key_dataframes:
        if key in analysis_results and analysis_results[key] is not None:
            consolidated_results[key] = analysis_results[key]
            print(f"Added {key}: {len(analysis_results[key])} rows")

    # Add marker groups with detailed logging
    print("\n=== Marker Groups ===")
    for group_name, group_df in marker_groups.items():
        consolidated_results[f'group_{group_name}'] = group_df
        print(f"Added group_{group_name}: {len(group_df)} rows, {len(group_df.columns)} columns")

    # Check specifically for complexation groups
    complexation_groups = ['inorganic_acid_complexes', 'organic_acid_complexes', 'biofilm_complexes']
    print("\n=== Complexation Groups Check ===")
    for cg in complexation_groups:
        if f'group_{cg}' in consolidated_results:
            print(f"✓ {cg}: {len(consolidated_results[f'group_{cg}'])} markers")
        else:
            print(f"✗ {cg}: NOT FOUND in marker_groups")

    print(f"\nSuccessfully consolidated {len(consolidated_results)} dataframes")
    notify_complete()
    return consolidated_results


# In[98]:


## calling main function
analysis_results, marker_groups = analyze_corrosion_proteins(ECcontri_Uniprot_Enriched, genus_to_threshold=genus_to_threshold, per_genus_count=10)

# Making the report #140
corrosion_report = consolidate_analysis_results(analysis_results, marker_groups)


# In[99]:


import psutil
print(psutil.virtual_memory())


# In[100]:


def print_group_shapes(corrosion_report):
    """
    Print the shapes of all DataFrames in the consolidated corrosion report.

    Parameters:
        corrosion_report: Dictionary containing all analysis results and marker groups
    """
    print("DataFrame shapes in the consolidated corrosion report:")
    for key, df in corrosion_report.items():
        if isinstance(df, pd.DataFrame) and df is not df.empty:
            sites = df['Sites'].nunique() if 'Sites' in df.columns else 'N/A'
            genus = df['Genus'].nunique() if 'Genus' in df.columns else 'N/A'
            print(f" {key:<35} Sites: {sites: < 5} Genus: {genus: < 5} Rows: {len(df): < 6}  Columns: {len(df.columns)}")
print_group_shapes(corrosion_report)


# In[101]:


# Main analysis frames
pattern_data           = corrosion_report["pattern_data"]            # Sites:  70   Genus:  85   Rows:  1168944  Columns: 21
integrated_results     = corrosion_report["integrated_results"]      # Sites:  70   Genus:  85   Rows:  1168944  Columns: 39
classified_results     = corrosion_report["classified_results"]      # Sites:  70   Genus:  85   Rows:  1168944  Columns: 42
increasing_markers     = corrosion_report["increasing_markers"]      # Sites:  53   Genus:  84   Rows:  849769  Columns: 42
inverse_markers        = corrosion_report["inverse_markers"]         # Sites:  17   Genus:  80   Rows:  319175  Columns: 42
prioritized_markers    = corrosion_report["prioritized_markers"]     # Sites:  53   Genus:  84   Rows:  849769  Columns: 51
balanced_markers       = corrosion_report["balanced_markers"]        # Sites:  53   Genus:  84   Rows:  779    Columns: 49

# Marker groups (accessing keys stored as 'group_<name>' in the report)
top_markers                  = corrosion_report["group_top_markers"]                   # Sites:  34   Genus:  42   Rows:  200    Columns: 51
significant_markers          = corrosion_report["group_significant_markers"]           # Sites:  21   Genus:  60   Rows:  467    Columns: 51
high_metals_relevance        = corrosion_report["group_high_metals_relevance"]       # Sites:  47   Genus:  54   Rows:  246    Columns: 51
high_functional_relevance    = corrosion_report["group_high_functional_relevance"]   # Sites:  48   Genus:  62   Rows:  263    Columns: 51
high_synergy_relevance       = corrosion_report["group_high_synergy_relevance"]      # Sites:  49   Genus:  69   Rows:  423    Columns: 51
high_corrosion_relevance     = corrosion_report["group_high_corrosion_relevance"]    # Sites:  46   Genus:  52   Rows:  199    Columns: 51
high_niche_relevance         = corrosion_report["group_high_niche_relevance"]        # Sites:  53   Genus:  84   Rows:  769    Columns: 51
mechanisms_all               = corrosion_report["group_mechanisms_all"]            # Sites:  53   Genus:  76   Rows:  677    Columns: 51
consolidated_metals          = corrosion_report["group_consolidated_metals"]       # Sites:  53   Genus:  77   Rows:  753    Columns: 51
pathways_all                 = corrosion_report["group_pathways_all"]             # Sites:  53   Genus:  84   Rows:  769    Columns: 51
functional_all               = corrosion_report["group_functional_all"]           # Sites:  53   Genus:  72   Rows:  466    Columns: 51

# Complexation / environmental groups (keep left side names without 'group_')
inorganic_acid_complexes     = corrosion_report["group_inorganic_acid_complexes"]   # Sites:  53   Genus:  73   Rows:  599    Columns: 54
organic_acid_complexes       = corrosion_report["group_organic_acid_complexes"]     # Sites:  30   Genus:  42   Rows:  132    Columns: 54
biofilm_complexes            = corrosion_report["group_biofilm_complexes"]          # Sites:  11   Genus:  7    Rows:  13     Columns: 54

# Operational, synergy and biological relevance groups
operational_all              = corrosion_report["group_operational_all"]          # Sites:  53   Genus:  77   Rows:  748    Columns: 54
synergy_all                  = corrosion_report["group_synergy_all"]              # Sites:  49   Genus:  69   Rows:  423    Columns: 54
high_biological_relevance    = corrosion_report["group_high_biological_relevance"]  # Sites:  48   Genus:  59   Rows:  244    Columns: 54
corrosion_critical           = corrosion_report["group_corrosion_critical"]        # Sites:  25   Genus:  15   Rows:  50     Columns: 54


# Checking for correct columns in grops 

# In[102]:


# Pipeline validation utility
# Usage:
#   report = check_pipeline_success(corrosion_report, verbose=True, n_examples=10)
#
# Assumes: pandas as pd, numpy as np, normalize_listlike, metal_mapping, category_dict already imported.

def check_pipeline_success(corrosion_report: dict, verbose: bool = True, n_examples: int = 10):
    """
    Run a set of defensive checks on the output dictionary produced by the pipeline (corrosion_report).
    Prints colored diagnostics and returns a structured dict with pass/fail booleans and details.

    What it checks (non-destructively):
      - expected keys/groups exist and are DataFrames
      - balanced_markers exists and has required columns
      - consolidated_metals is list-like for rows (counts non-list-like examples)
      - synergy_child_list / synergy_sub_list / synergy_combi are list-like
      - functional_child / mechanisms_child are strings (or empty)
      - combi columns presence & type-distribution (dict vs list_of_dicts)
      - combined_score and key numeric columns are finite numeric
      - Category column present and mapped (non-null fraction)
      - quick counts for top groups and sensible sizes (non-zero)
    Returns:
      dict with boolean 'ok' plus per-check details for programmatic assertions.
    """
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    RESET = '\033[0m'

    result = {'ok': True, 'details': {}}

    # Basic structure checks
    expected_groups = [
        'pattern_data', 'integrated_results', 'classified_results',
        'increasing_markers', 'inverse_markers', 'prioritized_markers',
        'balanced_markers', 'group_top_markers', 'group_significant_markers',
        'group_high_metals_relevance', 'group_high_functional_relevance',
       'group_mechanisms_all', 'group_consolidated_metals',
        'group_pathways_all', 'group_functional_all', 'group_inorganic_acid_complexes',
        'group_organic_acid_complexes', 'group_biofilm_complexes', 'group_operational_all',
        'group_synergy_all', 'group_high_biological_relevance', 'group_corrosion_critical'
    ]
    present = set(corrosion_report.keys())
    missing = [k for k in expected_groups if k not in present]
    result['details']['missing_groups'] = missing
    if missing:
        result['ok'] = False
        if verbose:
            print(f"{YELLOW}Missing expected groups: {missing}{RESET}")
    else:
        if verbose:
            print(f"{GREEN}All expected top-level groups present (checked list).{RESET}")

    # Balanced markers must exist and be a DataFrame
    if 'balanced_markers' not in corrosion_report:
        result['details']['balanced_markers'] = 'MISSING'
        result['ok'] = False
        if verbose:
            print(f"{RED}balanced_markers not present in report - pipeline failed early{RESET}")
        return result

    balanced = corrosion_report['balanced_markers']
    if not hasattr(balanced, 'shape') or not hasattr(balanced, 'columns'):
        result['details']['balanced_markers'] = 'NOT_DF'
        result['ok'] = False
        if verbose:
            print(f"{RED}balanced_markers is not a DataFrame-like object{RESET}")
        return result

    # Required columns we expect downstream (from schema)
    required_cols = [
        'Sites', 'idx', 'Genus', 'protein_name', 'mean_cat1', 'mean_cat2', 'mean_cat3',
        'combined_score', 'consolidated_metals', 'synergy_child_list', 'synergy_sub_list',
        'functional_child', 'mechanisms_child', 'operational_combi', 'functional_combi',
        'mechanisms_combi', 'pathways', 'niche_specific_pathways', 'Category'
    ]
    missing_cols = [c for c in required_cols if c not in balanced.columns]
    result['details']['missing_columns'] = missing_cols
    if missing_cols:
        result['ok'] = False
        if verbose:
            print(f"{YELLOW}Missing expected balanced_markers columns (examples): {missing_cols}{RESET}")
    else:
        if verbose:
            print(f"{GREEN}All required columns present in balanced_markers.{RESET}")

    # Helper: safe checks for list-like or string
    def _is_listlike(x):
        try:
            if x is None:
                return False
            # strings should NOT be listlike here
            if isinstance(x, str):
                return False
            return isinstance(x, (list, tuple, set, np.ndarray, pd.Series))
        except Exception:
            return False

    # Check consolidated_metals entries list-like
    non_list_mask = balanced['consolidated_metals'].apply(lambda x: not _is_listlike(x))
    non_list_count = int(non_list_mask.sum())
    result['details']['consolidated_metals_non_list_count'] = non_list_count
    if non_list_count > 0:
        result['ok'] = False
        if verbose:
            examples = balanced.loc[non_list_mask, 'consolidated_metals'].head(n_examples).tolist()
            print(f"{YELLOW}consolidated_metals has {non_list_count} non-list-like entries. Examples:{RESET}")
            for ex in examples:
                print("   ", repr(ex))
    else:
        if verbose:
            print(f"{GREEN}consolidated_metals appears list-like for all rows.{RESET}")

    # Check synergy lists are list-like (if present)
    for col in ('synergy_child_list', 'synergy_sub_list', 'synergy_combi'):
        if col in balanced.columns:
            non_list = balanced[col].apply(lambda x: not _is_listlike(x) and not (isinstance(x, str) and x.strip() == ''))
            cnt = int(non_list.sum())
            result['details'][f'{col}_non_list_count'] = cnt
            if cnt and verbose:
                print(f"{YELLOW}{col}: {cnt} rows are not list-like (or empty-string). Example values:{RESET}")
                for i, v in balanced.loc[non_list, col].head(n_examples).items():
                    print("   ", i, repr(v))
            elif verbose:
                print(f"{GREEN}{col} looks list-like or empty in all sampled rows.{RESET}")

    # Check functional_child and mechanisms_child are strings or empty
    for col in ('functional_child', 'mechanisms_child'):
        if col in balanced.columns:
            non_str_mask = balanced[col].apply(lambda x: not isinstance(x, str) and not pd.isna(x))
            cnt = int(non_str_mask.sum())
            result['details'][f'{col}_non_str_count'] = cnt
            if cnt:
                result['ok'] = False
                if verbose:
                    print(f"{YELLOW}{col}: {cnt} non-string values. Examples:{RESET}")
                    for i, v in balanced.loc[non_str_mask, col].head(n_examples).items():
                        print("   ", i, repr(v))
            else:
                if verbose:
                    print(f"{GREEN}{col} is string-like (or NaN) across rows.{RESET}")

    # Combi columns: report type distribution (dict vs list_of_dicts vs None vs str)
    combi_cols = [c for c in ('functional_combi', 'mechanisms_combi', 'operational_combi') if c in balanced.columns]
    combi_summary = {}
    for col in combi_cols:
        def _class(x):
            if x is None or (isinstance(x, float) and pd.isna(x)):
                return 'None'
            if isinstance(x, dict):
                return 'dict'
            if isinstance(x, (list, tuple)):
                if any(isinstance(i, dict) for i in x):
                    return 'list_of_dicts'
                return 'list_other'
            if isinstance(x, str):
                return 'str'
            return type(x).__name__
        types = balanced[col].map(_class).value_counts().to_dict()
        combi_summary[col] = types
        if verbose:
            print(f"{BLUE}{col} types: {types}{RESET}")
    result['details']['combi_types'] = combi_summary

    # Numeric columns sanity checks (finite & not-all-NaN)
    numeric_cols = ['combined_score', 'database_combined_score', 'overall_functional_score', 'overall_metal_score', 'overall_synergy_score', 'p_value', 'p_value_score', 'frequency_score', 'fc_score']
    numeric_issues = {}
    for col in numeric_cols:
        if col in balanced.columns:
            ser = pd.to_numeric(balanced[col], errors='coerce')
            n_nan = int(ser.isna().sum())
            n_inf = int(np.isinf(ser).sum())
            numeric_issues[col] = {'n_nan': n_nan, 'n_inf': n_inf, 'n_total': int(len(ser))}
            if (n_nan > 0 or n_inf > 0) and verbose:
                print(f"{YELLOW}Numeric sanity: {col} -> NaN={n_nan}, Inf={n_inf}{RESET}")
            if (n_nan > 0 or n_inf > 0):
                result['ok'] = False
    result['details']['numeric_issues'] = numeric_issues

    # Category mapping check (presence + non-null fraction)
    if 'Category' in balanced.columns:
        cat_nonnull = balanced['Category'].notna().sum()
        cat_total = len(balanced)
        result['details']['category_mapped_count'] = int(cat_nonnull)
        result['details']['category_total'] = int(cat_total)
        frac = cat_nonnull / cat_total if cat_total else 0.0
        if verbose:
            print(f"{BLUE}Category mapped: {cat_nonnull}/{cat_total} rows ({frac:.1%}){RESET}")
        # If many rows unmapped warn
        if frac < 0.5:
            result['ok'] = False
            if verbose:
                print(f"{YELLOW}Less than 50% of rows have Category mapped. Consider re-mapping Sites->{RESET}")
    else:
        result['details']['category_mapped_count'] = 0
        result['details']['category_total'] = int(len(balanced))
        result['ok'] = False
        if verbose:
            print(f"{RED}Category column missing from balanced_markers. Map Sites->Category before plotting.{RESET}")

    # Quick group-size sanity for several expected result groups (if present)
    group_checks = {}
    for g in ['group_top_markers', 'group_significant_markers', 'group_consolidated_metals', 'group_inorganic_acid_complexes', 'group_corrosion_critical']:
        if g in corrosion_report:
            try:
                dfg = corrosion_report[g]
                group_checks[g] = {'n_rows': int(len(dfg)), 'n_cols': int(dfg.shape[1]) if hasattr(dfg, 'shape') else None}
                if len(dfg) == 0:
                    result['ok'] = False
                    if verbose:
                        print(f"{YELLOW}Group {g} present but empty (0 rows).{RESET}")
            except Exception as e:
                group_checks[g] = {'error': str(e)}
                result['ok'] = False
        else:
            group_checks[g] = {'missing': True}
            result['ok'] = False
            if verbose:
                print(f"{YELLOW}Expected plotting group missing: {g}{RESET}")
    result['details']['group_checks'] = group_checks

    # Final summary
    if result['ok']:
        if verbose:
            print(f"{GREEN}Pipeline checks PASSED: the corrosion_report looks consistent for plotting.{RESET}")
    else:
        if verbose:
            print(f"{RED}Pipeline checks FAILED: see details in returned dict result['details'] for actionable items.{RESET}")

    return result


# In[103]:


print(top_markers[ ['functional_combi', 'mechanisms_combi', 'operational_combi', 'synergy_combi']].sample(5))


# In[104]:


result = check_pipeline_success(corrosion_report, verbose=True, n_examples=10)
notify_complete()


# In[105]:


specificity_report = generate_specificity_report(classified_results, output_file=None)
print(specificity_report)


# Let us remember that this is the refractoring of the original design, previously the databases and ec_records with the metadata information(pathways, hierarchy, mechanisms, functional categories, etc) were extracted directly from the public databases and as such the number of pathways were general and less granular. In this new implementation, the refractoring was done parsing the granular pathways and reactions coming from the data analysis using picrust2 and hence more granular pathways are accounted for. That implementation is in addition to the improved "corrosion system" for mining data and evaluating its relation to corrorsion. In the updated analysis, the proportion of mixed pathways has drastically decreased, while the proportion of niche-specific pathways has substantially increased. This shift is likely due to improvements in the pathway classification process and the increased granularity of the new dataset. More pathways can now be confidently classified as niche-specific, rather than being grouped into the “mixed” or “unknown” categories.
# 
# This is a summary of the pathways classification done in section 10.3 Classifying Housekeeping, Niche and Mixed Protein depending on pathways, mechanisms and hierarchy. The pathways by specificity resulted on :
# ## 4.9 Classification Summary on Prioritized Markers
# 
# Let us remember that this is the refractoring of the original design, previously the databases and ec_records with the metadata information(pathways, hierarchy, mechanisms, functional categories, etc) were extracted directly from the public databases and as such the number of pathways were general and less granular. In this new implementation, the refractoring was done parsing the granular pathways and reactions coming from the data analysis using picrust2 and hence more granular pathways are accounted for. That implementation is in addition to the improved "corrosion system" for mining data and evaluating its relation to corrorsion. In the updated analysis, the proportion of mixed pathways has drastically decreased, while the proportion of niche-specific pathways has substantially increased. This shift is likely due to improvements in the pathway classification process and the increased granularity of the new dataset. More pathways can now be confidently classified as niche-specific, rather than being grouped into the “mixed” or “unknown” categories.
# 
# This is a summary of the pathways classification done in section 10.3 Classifying Housekeeping, Niche and Mixed Protein depending on pathways, mechanisms and hierarchy. 
# For reference, the previous summary of pathway specificity (section 10.3) was:
#   - niche-specific: 1,097,851 (73.6%)
#   - mixed: 330,939 (22.2%)
#   - unknown: 62,495 (4.2%)
# 
# In the updated analysis, based on the improved methodology, the classification now is as follows:
# 
# === PATHWAY SPECIFICITY ANALYSIS ===
# Total pathways analyzed: 1390007
# Classification distribution:
#   - niche-specific: 1008798 (86.3%)
#   - universal: 159014 (13.6%)
#   - mixed: 1132 (0.1%)
# 
# Note that the genera here study was preselected to be the 85 most influencial genera out 887 total organisms types.
# 

# The predominance of niche-specific pathways suggests that corrosion environments select for specialized metabolic functions rather than just general housekeeping processes. This makes biological sense given the unique challenges of metal-rich, potentially oxygen-limited environments where corrosion occurs.
# 
# The universal pathways identified represent core cellular functions essential for any microorganism's survival (translation, glycolysis, TCA cycle, etc.), which explains their presence across the dataset. The distribution of pathway and reaction abundance (see sections 7.3 and 7.4) showed at a general level that abundance varied according to the risk level; this may differ at a more granular level, which is currently under investigation.
# 
# Many of the niche-specific examples directly relate to processes that could influence corrosion, such as:
# 
#     Methane metabolism (potentially creating anaerobic conditions)
#     Nitrogen metabolism (can affect local pH and redox conditions)
#     Metal-related metabolic pathways (directly interacting with metal surfaces)
#     Various secondary metabolite biosynthesis pathways (which may produce compounds that interact with metals)
# 
# This analysis corroborates the pre-selected 85 most influential genera (of the total analyzed), helping focus on the most relevant players in corrosion processes.
# Universal pathway breakdown (top examples, new data):
# 
# Universal pathway distribution (without duplication):
#   - 6.1. Amino Acid Biosynthesis: 80203 occurrences
#   - 4.3. Cell Wall Maintenance: 28404 occurrences
#   - 1.4. Krebs/TCA Cycle: 13331 occurrences
#   - 3.6. Amino Acid Biosynthesis: 10167 occurrences
#   - 1.1. Glycolysis: 8900 occurrences
#   - 1.6. Fermentation: 7819 occurrences
# 

# ### Saving the Results

# In[106]:


# Create directory for parquet files
parquet_dir = os.path.join(output_large, "Markers") 
os.makedirs(parquet_dir, exist_ok=True) 

# Save each DataFrame in the corrosion report as a parquet file 
for name, df in corrosion_report.items():
    parquet_path = os.path.join(parquet_dir, f"{name}.parquet")
    df.to_parquet(parquet_path, engine='pyarrow', compression='snappy')     
    print(f"Saved {name}.parquet - Size: {os.path.getsize(parquet_path) / (1024**2):.2f} MB")


# ## 5. Summary of the Results  
# This analysis pipeline implements a comprehensive approach to identify and prioritize protein-genus pairs associated with microbial-influenced corrosion across varying risk environments. By integrating statistical pattern recognition with biological metadata, the pipeline systematically progresses from identifying abundance patterns to classifying pathways, prioritizing markers, and organizing results into specialized groups. The methodology prioritizes biological relevance through metrics like prevalence, specificity, and frequency, moving beyond traditional statistical tests to better capture the complex relationships between microbial proteins and corrosion processes. Each output dataset serves a specific purpose in refining our understanding of the microbial mechanisms underlying corrosion, ultimately producing a curated set of high-confidence biomarkers with demonstrated relevance to corrosion severity. The following sections detail each component of the results and their contributions to the overall analysis framework.
# Analysis Files Description  
# * pattern_data: Previously included traditional statistical hypothesis testing with Kruskal-Wallis H statistics and p-values, however the data showed that this test were no adequated to the data. Therefore, the methodology has been refined to focus on more direct category-specific metrics: prevalence (proportion of samples in a category where a protein appears), specificity (how unique a protein is to certain categories), and frequency (how often a protein appears). This change provides more direct biological interpretability and relevance to corrosion environments.The pattern analysis has been refined to focus specifically on three key patterns: 'increasing_abundance' (proteins that increase with corrosion severity), 'only_cat3' (proteins exclusive to high corrosion environments), and 'only_cat2' (proteins exclusive to medium corrosion environments). This focused approach allows for more targeted identification of corrosion-relevant markers.
# 
# * integrated_results combines pattner data with functional metadata from the eccontri_df. This dataset enriches statistical findings with biological context including enzyme names, EC numbers, metabolic pathways, functional hierarchies, metal interactions, and corrosion-related mechanisms. It serves as the foundation for subsequent analysis by linking statistical significance to biological function.
# 
# * classified_results Annotates integrated results with pathway classifications using the comprehensive taxonomy of universal vs. niche-specific pathways. Each entry is labeled as "universal" (common across all bacteria/archaea), "niche-specific" (specialized for particular environments), or "mixed" (containing elements of both). This classification helps distinguish between housekeeping functions and potentially corrosion-specific adaptations. classification_summary is a quantitative breakdown of pathway classifications with counts and percentages for each category, providing insight into the distribution of universal vs. specialized pathways in corrosion environments. Universal_frequency is a detailed frequency analysis of specific universal pathways, showing which core metabolic functions are most prevalent in the dataset.
# 
# * Correlation Separated Results has tree types of results increasing_results which is a subset of integrated results showing positive correlation with corrosion risk.inverse_results which is a subset showing inverse correlation with corrosion risk, These are preserved for analysis of potentially protective or competitive proteins that might inhibit corrosion and analyse on a different section. Lastly constant_results a subset showing no significant correlation with risk categories.  
# 
# * prioritized_markers comprising markers ranked according to a comprehensive scoring system that integrates: statistical pattern significance, biological relevance (mechanisms, pathways, hierarchies), prevalence, specificity, frequency, and corrosion-specific factors (metal interactions, synergistic mechanisms). The scoring components include metals_score, functional_mechanisms_score, pathways_score, hierarchy_score, tier_score, and corrosion_relevance.
# 
# * balanced_markers is an optional filtered subset of prioritized markers that ensures balanced representation across genera, preventing over-representation of abundant genera while maintaining significant proteins.
# * marker_groups are categorized groups from the top 100 markers at 75% confidence threshold, organized by biological function and relevance.  
# 
# * Marker_groups: The marker groups organize proteins based on multiple biological and statistical criteria, providing different perspectives on corrosion-relevant functions.  
# 
# Main analysis dataframes: Core pipeline data from initial pattern detection through final balanced selection
# Top marker groups: Highest-scoring proteins and statistically significant patterns
# Relevance groups: Proteins with high scores in specific biological domains
# Mechanism-specific groups: Proteins associated with particular corrosion mechanisms
# Metal-related groups: Proteins interacting with specific metals relevant to corrosion
# Pathway classification groups: Proteins categorized by ecological specialization
# Effect size groups: Proteins showing dramatic abundance changes between corrosion states
# 
# This pipeline provides a comprehensive analysis framework that progresses from statistical identification to biological classification and prioritization, creating both comprehensive and focused views of corrosion-relevant microbial proteins

# # 6. High Confidence Corrosion Microorganisms
# ## 6.1. Retrieving Selected Significant Groups
# 
# From notebook 3_Feature_selection the file finalist.xlsx contain the groups worked and that were statistically significant in relation to the risk label. This groups posses interest since the relationship to the label could show better understanding in contrast with the different groups of known bacteria, core taxa, checked bacteria and the mixed groups between the former.
# The idea is to understand if the core taxa which make up a large influence on the comunities on the water and cooling systems are also influencing corrosion. The known_bacteria_list corresponds to the selected group found on notebook 4, where we conducted an exhaustive api literature review and get the most referenced bacteria to be in this group. In section ## 3.1. Classifying Bacteria by their Source DataFrame, we extracted the list for each of these groups: "known_bacteria", "pure_checked", "pure_core" and "checked_core".
# 
# We can verify the corrosion label "Category" by correlating the known_bacteria_list with the physicochemical parameters and the agreement if the known_bacteria_list have major abundance or correlate on someway with no only the label but the physicochemical conditions

# In[107]:


# Read to Excel with one sheet per category, this output base is located inside data_picrust which is the propietary folder of notebooks 6 and 7 
bacteria_clas_path = output_base / "bacteria_classification.xlsx"
usual_bacteria = pd.read_excel(bacteria_clas_path, sheet_name='usual_bacteria', engine ='openpyxl')
components_bacteria = pd.read_excel(bacteria_clas_path, sheet_name='components', engine ='openpyxl')
# Extracting the genus lists for each group:
usual_list = usual_bacteria["Genus"].tolist()
components_list = components_bacteria["Genus"].tolist()


# ## 6.2. Usual MIC bacteria

# In[108]:


# identify which genus are on my original data belogning to the known bacteria list, but genus are pair to protein_names
usual_eccontri = ECcontri_Uniprot_Enriched[ECcontri_Uniprot_Enriched['Genus'].astype(str).isin(usual_list)]
# Filter for components list
component_eccontri = ECcontri_Uniprot_Enriched[ECcontri_Uniprot_Enriched['Genus'].astype(str).isin(components_list)]


# In[109]:


## calling it with balancing the genera complete_results def balance_genus_representation(prioritized_markers, eccontri_df, genus_to_threshold, per_genus_count=10)
analysis_results_us, marker_groups_us = analyze_corrosion_proteins(usual_eccontri, genus_to_threshold=genus_to_threshold, per_genus_count=10)

# Making the report
corrosion_report_us = consolidate_analysis_results(analysis_results_us, marker_groups_us)


# In[110]:


# Create directory for parquet files
parquet_dir_us = os.path.join(output_large, "Usuals")
os.makedirs(parquet_dir_us, exist_ok=True)

# Save each dataframe as a separate parquet file
for name, df in corrosion_report_us.items():
    parquet_path = os.path.join(parquet_dir_us, f"{name}.parquet")
    df.to_parquet(parquet_path, engine='pyarrow', compression='snappy')     
    print(f"Saved {name}.parquet - Size: {os.path.getsize(parquet_path) / (1024**2):.2f} MB")


# ## 6.3. Component MIC bacteria
# The group components corresponds to all the bacteria that is part of the usual bacteria and the check bacteria that for this study has shown correlation.

# In[111]:


## calling it with balancing the generacomplete_results def balance_genus_representation(prioritized_markers, eccontri_df, genus_to_threshold, per_genus_count=10)
analysis_results_component, marker_groups_component = analyze_corrosion_proteins(component_eccontri, genus_to_threshold=genus_to_threshold, per_genus_count=10)

# Making the report
corrosion_report_component = consolidate_analysis_results(analysis_results_component, marker_groups_component)


# In[112]:


# Create directory for parquet files
parquet_dir_component = os.path.join(output_large, "Components")
os.makedirs(parquet_dir_component, exist_ok=True)

# Save each dataframe as a separate parquet file
for name, df in corrosion_report_component.items():
    parquet_path = os.path.join(parquet_dir_component, f"{name}.parquet")
    df.to_parquet(parquet_path, engine='pyarrow', compression='snappy')     
    print(f"Saved {name}.parquet - Size: {os.path.getsize(parquet_path) / (1024**2):.2f} MB")

