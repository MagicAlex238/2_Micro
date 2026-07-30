#!/usr/bin/env python
# coding: utf-8

# # Introduction
# 
# Methodology: Taxonomic Harmonization and Cross-Batch Integration
# I actually have done the harmonisation myself, however I realise it has to be documanted on line as you did, I only did color code, I wrote the method but It has to be double checked. 
# Taxonomic Standardization and Structural Restructuring
# To enable robust cross-batch integration and statistical comparison across multi-year sequencing projects, a standardized taxonomic nomenclature was enforced. Raw taxonomic lineages were systematically cleaned to remove low-confidence annotations beyond the genus level, discarding all species-level designations. Standard nomenclature formats were applied uniformly across all datasets using consistent capitalization and joining multi-word terms with underscores.
# 
# To optimize algorithmic sorting and indexing, the conventional hierarchical prefix notation was inverted. Lineages failing to resolve to lower taxonomic ranks were refactored from a prefix-dependent format (e.g., unclassified_FamilyName) to a suffix-dependent format (e.g., FamilyName_unclassified). This restructuring allowed for systematic, alphabetical grouping and matrix sorting strictly by hierarchy, moving sequentially from Kingdom down through Phylum, Class, Order, and Family. Missing phylum-level assignments were subsequently harmonized and remapped to modern phylum designations.
# 
# A Priori Baseline Selection and Environmental Filtering
# Cross-batch harmonization was guided by a deeply curated, high-confidence baseline taxonomy dataset comprising approximately 800 recognized genera. This baseline dataset, established via rigorous biological curation, was designated as the target structural framework. Given that samples spanned disparate computational workflows, laboratories, and collection timelines, the core microbial communities were expected to mirror the specialized ecological niches of the engineering systems under study.
# 
# To neutralize batch effects and eliminate false diversity expansions driven by technical artifacts or laboratory-specific misclassifications, unmapped or newly introduced taxa were managed using a conservative ecological filtering protocol. Unidentified or newly introduced operational taxonomic units (OTUs) were systematically collapsed into known baseline genera under the following restrictive conditions:
# 
# The target organism shared a direct, identical Family assignment with an established baseline representative.
# 
# The biological and ecological traits of the candidate genus aligned strictly with the specific environmental profiles of closed-loop industrial heating and cooling water infrastructure.
# 
# Taxa whose primary metadata tied them exclusively to disparate environmental systems—such as marine, deep-ocean, or petroleum reservoirs—were classified as cross-batch noise and omitted from direct genus-level mapping.
# 
# Abundance Thresholding, Merging, and Mass Conservation
# When the incoming sequencing matrices introduced distinct taxonomic splits (e.g., separating an established baseline genus into multiple newly resolved variants), data integrity was maintained via abundance pooling. The abundance counts of these newly split taxa were collapsed back into the dominant, high-confidence baseline representative by summing their respective cell values, ensuring total absolute mass conservation within the sample profiles.
# 
# To protect the dataset from inflation by minor sequencing artifacts while retaining significant biological signals, unclassified lineages at the family level were preserved as valid analytical placeholders only if their individual relative abundance exceeded an ecological significance threshold of greater than 0.01%. High-abundance unclassified fractions critical to specific matrix profiles (such as the dominant unspecific_Bacteria_meta2 faction, representing 66% of Sample 1 total reads) were strictly locked and retained to prevent mathematical distortion of the remaining relative abundance profiles. All manual fusions, data corrections, and abundance shifts were permanently color-flagged within the master matrix to ensure complete traceability. Orange for the Genus representing the 2 or more genus and pink for the cell which get absorbed by the orange. Lila for the genus likely to become the representing genus. 
# 
#                      MIC PROJECT  
#                         │  
#           ┌─────────────┴─────────────┐  
#           │                           │  
#     Biotot_noncured ~800 GENERA    
#    (preserved, never modified)       NEW 150/300 GENERA  
#           │                           │  
#           └─────────────┬─────────────┘  
#                         │    
#                     MultiTax   
#                         │  
#                         ▼  
#                     GTDB R232  
#                         │  
#                         ▼  
#                 Current taxonomy 800-genus harmonised list  
#                         │  
#             ┌───────────┴───────────┐  
#             │                       │  
#     Update old reference     Classify new data  
#             │                       │  
#             ▼                       ▼  
#     GID identity preserved       Match / new / unknown  

# In[2]:


from pathlib import Path
import pandas as pd
from multitax import GtdbTx
import inspect
import re
import openpyxl
#gtdb = GtdbTx(version="232")
#print(gtdb)


# In[3]:


print(inspect.signature(GtdbTx))


# In[4]:


# Load  current harmonised 800-genus reference
biotot_path = Path("data/Biotot.xlsx")

biotot = pd.read_excel(biotot_path, sheet_name="Biotot", engine = "openpyxl")
biotot["GID"] = biotot["GID"].astype("Int16")
# print(biotot.head())


# In[5]:


def clean_genus_string(genus):
    """General cleanup: trim whitespace, collapse repeated/mixed
    whitespace and underscores into single underscores, strip
    leading/trailing underscores."""
    if pd.isna(genus):
        return genus
    genus = str(genus).strip()
    genus = re.sub(r"\s+", "_", genus)      # any whitespace -> underscore
    genus = re.sub(r"_+", "_", genus)       # collapse multiple underscores into one
    genus = genus.strip("_")                # drop leading/trailing underscores
    return genus

biotot["Genus"] = biotot["Genus"].apply(clean_genus_string)


# # Query GTDB R232 with MultiTax
# 
# Create a gtdb_taxonomy dataframe.
# Query data from the GTDB database, to get the taxonomic ranks for each genus

# In[6]:


genus_queries = biotot[["Genus", "GID"]].copy()


# In[7]:


# Query data from the GTDB database, to get the taxonomic ranks for each genus
def clean_gtdb_query(genus):
    """
    Convert the project's Genus label into the best GTDB query name,
    while leaving the original Genus column untouched.
    """
    if pd.isna(genus):
        return None

    genus = str(genus).strip()

    # Remove known placeholder suffixes
    genus = re.sub(
        r"_(unclassified|uncultured|unknown)$", "", genus, flags=re.IGNORECASE)

    # Remove strain/sequence-style numeric suffix:
    # Selenomonas_3 -> Selenomonas
    genus = re.sub(r"_\d+$", "", genus)

    return genus


# In[8]:


def get_query_rank(original_genus):
    '''ranking it so that it queries depending on the rank of the genus, if it is a family or genus'''
    if pd.isna(original_genus):
        return None

    original_genus = str(original_genus).strip()

    # Family-level unresolved names # Family-level placeholder:# Rhodocyclaceae_unclassified -> Rhodocyclaceae
    # This will then be queried as a family, not as a genus.
    if original_genus.endswith("aceae"):
        return "f__"
    if original_genus.endswith("_unclassified") or original_genus.endswith("_uncultured"):
        base = re.sub(
            r"_(unclassified|uncultured)$","", original_genus,
            flags=re.IGNORECASE
        )

        if base.endswith("aceae"):
            return "f__"

    # Numeric variants
    if re.search(r"_\d+$", original_genus):
        return "g__"

    # Normal genus
    return "g__"


# In[9]:


# the querry will have the Genus and keep the ID column, so to not lose the ID, id already in df
genus_queries["GTDB_query_genus"] = (genus_queries["Genus"].astype(str).map(clean_gtdb_query))

genus_queries["GTDB_query_rank"] = (genus_queries["Genus"].map(get_query_rank))

genus_queries["GTDB_query"] = (genus_queries["GTDB_query_rank"] + genus_queries["GTDB_query_genus"])

genus_queries["Query_transformation"] = (genus_queries["Genus"] != genus_queries["GTDB_query_genus"])


# In[10]:


rank_prefix = {"d__": "Kingdom", "p__": "Phylum", "c__": "Class", "o__": "Order", "f__": "Family", "g__": "Genus"}

def lookup_lineage(gtdb_query):
    """Takes a full query string (e.g. 'g__Escherichia' or 'f__Aerococcaceae').
    Returns a dict of rank -> value if a real lineage was found, else None."""
    try:
        lineage_nodes = gtdb.lineage(gtdb_query)
    except Exception:
        return None
    parsed = {}
    for node in lineage_nodes:
        for prefix, rankname in rank_prefix.items():
            if node.startswith(prefix):
                parsed[rankname] = node[len(prefix):]
                break
    return parsed if parsed else None
def search_gtdb_nodes(substring):
    """Case-insensitive substring search across all GTDB node keys.
    Used to find suffixed variants (Calothrix_A) or check whether
    something exists under a different name."""
    substring = str(substring).lower()
    return [n for n in gtdb._nodes.keys() if substring in n.lower()]

def lookup_lineage_with_suffix_retry(base_name, rank_prefix_char="g__"):
    """Try exact match first. If that fails, search for GTDB's suffixed
    variants (_A, _B, _A1, etc.) of the same base name and use the first one found."""
    base_name = str(base_name) 
    parsed = lookup_lineage(f"{rank_prefix_char}{base_name}")
    if parsed:
        return parsed, base_name

    candidates = search_gtdb_nodes(base_name)
    suffixed = [c for c in candidates if re.match(rf"^{rank_prefix_char}{re.escape(base_name)}_[A-Z]\d*$", c)]
    for candidate_node in suffixed:
        parsed = lookup_lineage(candidate_node)
        if parsed:
            resolved_name = candidate_node[len(rank_prefix_char):]
            return parsed, resolved_name

    return None, None


# In[11]:


records = []
for _, r in genus_queries.iterrows():
    original_genus = r["Genus"]
    row = {"GID": r["GID"], "Genus": original_genus}

    parsed, resolved_name = lookup_lineage_with_suffix_retry(original_genus)
    if parsed:
        row.update(parsed)
        row["GTDB_status"] = "direct" if resolved_name == original_genus else f"resolved_via_suffix({resolved_name})"
    else:
        stripped = clean_gtdb_query(original_genus)
        if stripped and stripped != original_genus:
            parsed, resolved_name = lookup_lineage_with_suffix_retry(stripped)
            if parsed:
                row.update(parsed)
                row["GTDB_status"] = f"resolved_via_base_genus({resolved_name})"
            else:
                row["GTDB_status"] = "not_in_gtdb"
        else:
            row["GTDB_status"] = "not_in_gtdb"

    records.append(row)  # always runs, every genus gets a row regardless of outcome

table_b = pd.DataFrame(records)


def strip_gtdb_suffix(name):
    """Burkholderiaceae_B -> Burkholderiaceae. Only strips a single trailing
    underscore + capital letter, GTDB's specific convention ."""
    if pd.isna(name):
        return name
    return re.sub(r"_[A-Z]$", "", str(name))

for rank in ["Kingdom", "Phylum", "Class", "Order", "Family"]:
    table_b[rank] = table_b[rank].apply(strip_gtdb_suffix)


# In[ ]:


table_b.sample(10)


# In[ ]:


merged = biotot[["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "GID"]].merge(
    table_b, on="GID", how="left", suffixes=("_original", "_gtdb"))


# In[ ]:


RANKS_TO_COMPARE = [
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family"
]

def make_comment(row):

    # No GTDB result
    if row["GTDB_status"] == "not_in_gtdb":
        return "No GTDB lineage retrieved"

    if row["GTDB_status"] == "skipped_placeholder":
        return "Placeholder/composite name; not directly checked against GTDB"

    diffs = []

    for rank in RANKS_TO_COMPARE:

        original = row.get(f"{rank}_original")
        gtdb_val = row.get(f"{rank}_gtdb")

        if pd.isna(gtdb_val):
            continue

        if pd.isna(original):
            diffs.append(f"{rank}: missing -> {gtdb_val}")

        elif str(original).strip() != str(gtdb_val).strip():
            diffs.append(f"{rank}: {original} -> {gtdb_val}")

    return "; ".join(diffs) if diffs else "Exact lineage match"

merged["Comment"] = merged.apply(make_comment, axis=1)
mismatches = merged[merged["Comment"].str.startswith(("Phylum", "Class", "Order", "Family"))]


# In[ ]:


merged.sample(10)


# In[ ]:


not_found = merged[merged["GTDB_status"] == "not_in_gtdb"].copy()

def categorize(genus):
    genus = str(genus)
    if any(m in genus.lower() for m in ("uncultured", "clone", "enrichment_culture", "bacterium_")):
        return "environmental_clone_or_placeholder"
    if genus.startswith("Candidatus_"):
        return "candidatus_name"
    if "_" in genus and re.match(r"^[A-Z][a-z]+_[A-Z][a-z]+$", genus):
        return "composite_two_genus_name"
    if re.match(r"^[A-Z][a-z]+$", genus):
        return "clean_genus_name"  # worth checking individually -- this is the group most likely a real gap or rename
    return "other"

not_found["category"] = not_found["Genus_original"].apply(categorize)
not_found["category"].value_counts()


# New data Parcing
# Found under revising some files and matching it to the late publication on 2024.

# In[ ]:


df = pd.read_csv("data/Linage_Family_from_Taxonomietabelle.csv", header=None)


# In[ ]:


# Combine Row 0 (TubeID) and Row 1 (FF... ID) to keep track of both sample identifiers
header_cols = []
for col in range(df.shape[1]):
    val1 = str(df.iloc[0, col]).strip() if pd.notna(df.iloc[0, col]) else ""
    val2 = str(df.iloc[1, col]).strip() if pd.notna(df.iloc[1, col]) else ""
    if val1 and val2:
        header_cols.append(f"{val1}_{val2}")
    elif val1:
        header_cols.append(val1)
    elif val2:
        header_cols.append(val2)
    else:
        header_cols.append(f"Col_{col}")
# Slice the dataframe to drop the first two header rows, then apply new column names
df_clean = df.iloc[2:].copy()
df_clean.columns = header_cols
df_clean.reset_index(drop=True, inplace=True)


# In[ ]:


# 3. Split the OTU ID from the Kingdom string in the first column
# "4362609.0 k__Bacteria" -> OTUID: "4362609.0", Kingdom: "k__Bacteria"
split_first_col = df_clean.iloc[:, 0].str.split(" ", n=1, expand=True)
df_clean.insert(0, "OTU_ID", split_first_col[0])
df_clean.iloc[:, 1] = split_first_col[1]  # Overwrite the original column with just the Kingdom string


# In[ ]:


# Rename taxonomy columns so they are clean and recognizable
tax_ranks = ["Kingdom", "Phylum", "Class", "Order", "Family"]
for i, rank in enumerate(tax_ranks):
    df_clean.rename(columns={df_clean.columns[i + 1]: rank}, inplace=True)
# removing the prefixes "k__", "p__"
tax_columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]

# Strip the prefixes (e.g., "k__", "p__") only from these columns
for col in tax_columns:
    if col in df_clean.columns:
        # This replaces any letter followed by '__' at the start of the string
        df_clean[col] = df_clean[col].astype(str).str.replace(r'^[a-z]__', '', regex=True)


# The lineage from Kingdom to Family was on a page that was separated from the Genus, so off course was difficult to match them because the reorganisation of the lineage by the experts has been continued a long the years, so it was just merged.

# In[ ]:


df_genus = pd.read_csv('data/lineage_genus_from_Taxonomietabelle.csv')
# merge family df with genus df
df_merged = pd.merge(df_clean, df_genus, on='Family') #, how=
# order the columns Otu, Kingdom, Phylum, Class, Order, Family Genus, rest
df_merged = df_merged[["OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"] + [col for col in df_merged.columns if col not in ["OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]]] 


# Now the new genus from this list are going to be harmonised with the original curated data

# In[ ]:


# saving the df to a excel file
df_merged.to_excel("data/merged_taxonomy_from_taxonomytabelle.xlsx", index=False)

