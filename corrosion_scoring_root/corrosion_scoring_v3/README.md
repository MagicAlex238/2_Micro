# Corrosion Scoring System v3.0

A refactored Python module for analyzing corrosion relevance in scientific texts through text mining, scoring, and synergy detection.

## Overview

This system provides a clean, modular architecture for:
- **Text Mining**: Extracting corrosion-related terms and categories from scientific texts
- **Scoring**: Calculating relevance scores based on functional categories, metals, and synergies
- **Synergy Detection**: Identifying co-occurrence patterns that indicate corrosion mechanisms

## Architecture

The system is built with separation of concerns in a nested module structure:

```
corrosion_scoring_root/
├── __init__.py                     # Main entry point with backward compatibility, I dont have
|
└── corrosion_scoring_V2/
|    |
|    ├── __init__.py                 # Term dictionaries and mappings
|    ├── global_terms.py             # Term dictionaries and mappings      
|    ├── scoring_system.py           # mining and scoring functions
|    ├── utils_ec.py                # ec name normalisation
|    ├── name_utils.py              # protein_name normalisation
|    └── term_processor.py           # sdg Term normalization and matching
|
└── corrosion_scoring_v3/
    ├── __init__.py
    ├── config.py                  # Configuration settings and parameters
    ├── exceptions.py              # Custom exception classes
    ├── text_miner.py              # Feature extraction (no scoring)
    ├── score_calculator.py        # Score computation from features
    ├── synergy_detector.py        # Synergy detection and analysis
    ├── term_processor.py          # sdg Term normalization and matching
    ├── validators.py              # Input validation utilities
    ├── utils_ec.py                # ec name normalisation
    ├── name_utils.py              # protein_name normalisation
    └── global_terms.py            # Term dictionaries and mappings
```

## Key Features

### 🔍 Text Mining
- **Functional Category Detection**: 16 categories including iron metabolism, biofilm formation, etc.
- **Metal Term Extraction**: Standardization and consolidation from text and BRENDA data
- **Mechanism Detection**: Corrosion mechanism pathway identification
- **Operational Factor Detection**: Environmental and operational conditions

### 📊 Scoring System
- **Weighted Scoring**: Configurable weights for functional categories, metals, and synergies
- **Priority-Based Processing**: Term conflicts resolved using priority order
- **Synergy Scoring**: Co-occurrence pattern analysis with priority synergy pairs

### 🔗 Synergy Detection
- **Co-occurrence Analysis**: Multi-category term pattern detection
- **Priority Synergies**: High-value pairs (e.g., TOC-chelation, Iron-Organic Acid)
- **Configurable Thresholds**: Minimum terms and categories for synergy detection

## Global Terms Structure

The system uses several dictionaries for term classification:

### Metal Terms Dictionary
```python
metal_terms_dict = {
    'iron': {
        'terms': ['Fe2+', 'Fe3+', 'iron', 'ferrous', 'ferric', 'heme', 'iron-sulfur', 'rust', 'ochre', 'iron oxide', 'siderophore', 'ferritin'],
        'score': 1.0
    },
    'copper': {
        'terms': ['Cu2+', 'Cu+', 'copper', 'cupric', 'cuprous', 'azurin', 'plastocyanin'],
        'score': 1.0
    }
    # ... additional metal categories
}
```

### Functional Categories Dictionary
```python
functional_categories_dict = {
    'iron_metabolism': {
        'terms': ['iron uptake', 'siderophore', 'ferritin', 'transferrin', 'iron transport'],
        'score': 2.0
    },
    'biofilm_formation': {
        'terms': ['biofilm', 'extracellular matrix', 'quorum sensing', 'adhesion'],
        'score': 1.8
    }
    # ... 14 additional categories
}
```

### Mechanisms Dictionary
```python
mechanisms_dict = {
    'o2_consumption': [
        '2-aminophenol degradation',
        'PWY-5088',
        'PWY-6953',
        'aerobic respiration',
        'oxygen reduction'
    ],
    'sulfur_oxidation': [
        'sulfur oxidation pathway',
        'thiosulfate oxidation',
        'elemental sulfur oxidation'
    ]
    # ... additional mechanism categories
}
```

## Quick Start

```python
from corrosion_scoring_v3 import calculate_overall_scores, ScoringConfig
from corrosion_scoring_v3.term_processor import TermProcessor
from corrosion_scoring_v3.global_terms import (
    functional_categories_dict, 
    metal_terms_dict, 
    corrosion_synergies_dict
)

# Initialize processors (required)
processors = {
    'fc_processor': TermProcessor(functional_categories_dict),
    'metal_processor': TermProcessor(metal_terms_dict),
    'synergy_processor': TermProcessor(corrosion_synergies_dict)
}

# Analyze text
text = "Iron metabolism and biofilm formation with Fe2+ and organic acid chelation"
brenda_metals = ['Fe2+', 'Cu2+']  # Optional metals from BRENDA database

results = calculate_overall_scores(text, processors, brenda_metals=brenda_metals)

print(f"Corrosion Relevance Score: {results['corrosion_relevance_score']}")
print(f"Functional Categories: {results['functional_categories_detected']}")
print(f"Detected Metals: {results['consolidated_metals']}")
print(f"Synergy Score: {results['synergy_score']}")
```

## Advanced Usage

### Custom Configuration

```python
config = ScoringConfig(
    metal_score_weight=0.5,
    functional_score_weight=1.5,
    synergy_score_weight=2.0,
    min_synergy_categories=2,
    min_synergy_terms=2,
    max_text_length=15000
)

results = calculate_overall_scores(text, processors, config=config)
```

### Direct Component Usage

```python
from corrosion_scoring_v3 import TextMiner, ScoreCalculator, SynergyDetector

# Initialize components separately
text_miner = TextMiner(processors)
score_calculator = ScoreCalculator(config)
synergy_detector = SynergyDetector(config)

# Extract features
features = text_miner.extract_all_features(text, brenda_metals)

# Calculate scores
scores = score_calculator.calculate_scores(features)

# Detect synergies
synergies = synergy_detector.detect_synergies(features)
```

## Priority Synergies

The system recognizes key corrosion synergies with weighted scoring:

| Synergy Pair | Score | Description |
|-------------|-------|-------------|
| organic_acid_metabolism + metal_binding_chelation | 3.0 | TOC-chelation Synergy |
| iron_metabolism + organic_acid_metabolism | 2.8 | Acid-enhanced Fe corrosion |
| biofilm_formation + metal_binding_chelation | 2.7 | Biofilm-chelate corrosion |
| o2_consumption + iron_metabolism | 2.5 | Aerobic Fe corrosion |
| sulfur_metabolism + iron_metabolism | 2.4 | SRB-mediated corrosion |

## Output Structure

The `calculate_overall_scores` function returns a comprehensive dictionary:

```python
{
    # Feature extraction results
    'functional_categories_detected': ['iron_metabolism', 'biofilm_formation'],
    'functional_child_terms': ['iron uptake', 'biofilm', 'siderophore'],
    'functional_matches_detailed': {
        'iron_metabolism': ['iron uptake', 'siderophore'],
        'biofilm_formation': ['biofilm']
    },
    
    # Metal analysis
    'consolidated_metals': ['Fe', 'Cu'],
    'detected_metal_categories': ['iron', 'copper'],
    'brenda_metals': ['Fe2+', 'Cu2+'],
    'metal_matches_detailed': {
        'iron': ['Fe2+', 'iron'],
        'copper': ['Cu2+', 'copper']
    },
    
    # Mechanism and pathway detection
    'corrosion_mechanisms': ['o2_consumption'],
    'detected_pathways': ['aerobic_respiration'],
    'operational_environmental_factors': ['ph_conditions'],
    
    # Scoring results
    'functional_score': 4.2,
    'metal_score': 2.0,
    'synergy_score': 2.8,
    'overall_functional_score': 6.3,
    'overall_metal_score': 1.0,
    'overall_synergy_score': 5.6,
    'corrosion_relevance_score': 12.9,
    
    # Synergy analysis
    'fc_cooccurrence_synergy_hit': True,
    'synergy_categories_involved': ['iron_metabolism', 'organic_acid_metabolism'],
    'synergy_description': 'Iron-Organic Acid Synergy (acid-enhanced Fe corrosion)',
    'synergy_child_terms_found': ['iron uptake', 'organic acid']
}
```

## Error Handling

The system includes comprehensive error handling with custom exceptions:

```python
from corrosion_scoring_v3.exceptions import (
    ScoringError, 
    TextMiningError, 
    SynergyDetectionError,
    ValidationError
)

try:
    results = calculate_overall_scores(text, processors)
except TextMiningError as e:
    print(f"Text mining failed: {e}")
except ScoringError as e:
    print(f"Scoring failed: {e}")
except ValidationError as e:
    print(f"Input validation failed: {e}")
```

## Functional Categories

The system recognizes 16 functional categories in priority order:

1. **iron_metabolism** - Iron uptake, transport, and processing
2. **sulfur_metabolism** - Sulfur oxidation/reduction pathways
3. **organic_acid_metabolism** - Organic acid production and processing
4. **biofilm_formation** - Biofilm development and maintenance
5. **o2_consumption** - Oxygen consumption processes
6. **metal_binding_chelation** - Metal chelation and binding
7. **nitrogen_metabolism** - Nitrogen cycle processes
8. **manganese_processes** - Manganese-related metabolism
9. **h2_consumption** - Hydrogen consumption
10. **methanogenesis** - Methane production
11. **carbon_metabolism** - Carbon cycle processes
12. **fumarate_formation** - Fumarate-related pathways
13. **phosphorus_metabolism** - Phosphorus cycle
14. **direct_eet** - Direct extracellular electron transfer
15. **indirect_eet** - Indirect extracellular electron transfer
16. **halogen_related** - Halogen-related processes

## Requirements

- Python 3.7+
- Standard library only (no external dependencies)
- Pre-configured processor instances with global term dictionaries

## Migration from v2.0

The v3.0 system maintains backward compatibility through the `calculate_overall_scores` function while providing improved modularity and error handling. Key improvements:

- Separated text mining from scoring calculations
- Enhanced synergy detection with priority pairs
- Comprehensive input validation
- Modular architecture for easier testing and maintenance

## Contributing

1. Ensure all tests pass
2. Follow the existing code style and error handling patterns
3. Add appropriate validation for new features
4. Update documentation for any API changes

## License

MIT License
____________________________________
# v3 SCORING SYSTEM TEXT MINING - SCORING 
____________________________________

cd /home/beatriz/MIC/2_Micro/corrosion_scoring_root/corrosion_scoring_v3
{  echo "# ========== __init__.py =========="
    cat __init__.py
    echo -e "\n\n# ========== config.py =========="
    cat config.py
    echo -e "\n\n# ========== exceptions.py =========="
    cat exceptions.py
    echo -e "\n\n# ========== global_terms.py ========== too large" # it is too large
    #cat global_terms.py
    echo -e "\n\n# ========== name_utils.py  =========="
    cat name_utils.py
    echo -e "\n\n# ========== utils_ec.py  =========="
    cat utils_ec.py
    echo -e "\n\n# ========== score_calculator.py =========="
    cat score_calculator.py
    echo -e "\n\n# ========== synergy_detector.py ========="
    cat synergy_detector.py
    echo -e "\n\n# ========== term_processor.py =========="
    cat term_processor.py
    echo -e "\n\n# ========== text_miner.py =========="
    cat text_miner.py
    echo -e "\n\n# ========== validators.py =========="
    cat validators.py
} > /home/beatriz/MIC/2_Micro/corrosion_scoring_root/corrosion_scoring_v3/.v3_summary_combined_files_cs.py
´´´

´´´
# Convert to Python script
jupyter nbconvert --to python 6_picrust_functional.ipynb
####_________________-
MD
'''
### Symplifying Protein_name values

The full pipeline was displayed with too many repeated names on the protein_name side, that is the reason why it was decided to simplify the protein_name keeping the biological meaning without removing too much. We remove words such as "probable", "putative", "possible", "uncharacterized", "hypothetical", etc.
→ These terms indicate annotation uncertainty but are not part of the function name. Additionaly EC was remove which are already available on a different column, if required. And since the original protein_name values contain multiple names that are suppose to be synonim and annotated according to the enviroment, it was necesary to take a decision to be able to leverage this information without overhelming the plot because of the long name. The selected name out of the whole variety is the first as widerly accepted or cannonical name. The spaces extra were cleaned whitespace at start/end and between words and it was lowercased. Several runs to the whole pipeline including the notebook 7 were necesary to really detect the names and the removals necesary. The first core name up to the first parenthesis or bracket only if what follows is clearly additional information (e.g., enzyme classification, isoforms) and keep all text if the parentheses are truly part of the main name.
Lastly the depest names were remove as redundant. These are selected examples 2,3 of the protein_name before
|idx| protein_name|
|--|--|
|2|   alcohol-dehydrogenase; aldehyde-reductase; adh; alcohol-dehydrogenase (nad); aliphatic alcohol-dehydrogenase; ethanol-dehydrogenase; nad-dependent alcohol-dehydrogenase; nad-specific aromatic alcohol-dehydrogenase; nadh-alcohol-dehydrogenase; nadh-aldehyde-dehydrogenase; primary alcohol-dehydrogenase; yeast alcohol-dehydrogenase|
|3 |  alcohol-dehydrogenase; aldehyde-reductase; adh; alcohol-dehydrogenase (nad); aliphatic alcohol-dehydrogenase; ethanol-dehydrogenase; nad-dependent alcohol-dehydrogenase; nad-specific aromatic alcohol-dehydrogenase; nadh-alcohol-dehydrogenase; nadh-aldehyde-dehydrogenase; primary alcohol-dehydrogenase; yeast alcohol-dehydrogenase|
|1491283|   cobaltochelatase; hydrogenobyrinic acid a,c-diamide cobaltochelatase; cobnst; cobncobst; hydrogenobyrinic-acid-a,c-diamide:cobalt cobalt-ligase (adp-forming)|
|1491285 |    cobaltochelatase; hydrogenobyrinic acid a,c-diamide cobaltochelatase; cobnst; cobncobst; hydrogenobyrinic-acid-a,c-diamide:cobalt cobalt-ligase dp-forming)|

Name: protein_name, Length: 1491288, dtype: object

And these are the simplified versions with which we continue:
| idx | protein_name|
|--|--|
|2  |alcohol-dehydrogenase|
|3    |alcohol-dehydrogenase|
|1491283|  cobaltochelatase|
|1491285| cobaltochelatase|

protein_name simplification as 
UniProt, NCBI RefSeq standards (The UniProt Consortium, 2021)
Text cleaning standard in biomedical data preparation 
(Müller et al., 2016, Introduction to Biomedical Text Mining).
Logical matching best practice in name cleaning (APA: UniProt Consortium, 2021).

References : 

Müller, H. M., Kenny, E. E., & Sternberg, P. W. (2016). Textpresso: An ontology-based information retrieval and extraction system for biological literature. PLoS Biology, 2(11), e309. https://doi.org/10.1371/journal.pbio.0020309

Source (for regex best practices):
Friedl, J. E. F. (2006). Mastering Regular Expressions (3rd ed.). O'Reilly Media.

McKinney, W. (2012). Python for Data Analysis. O'Reilly Media. (p.77-78)

Pruitt, K. D., Tatusova, T., & Maglott, D. R. (2007). NCBI reference sequences (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins. Nucleic Acids Research, 35(Database issue), D61–D65. https://doi.org/10.1093/nar/gkl842

The UniProt Consortium. (2021). UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Research, 49(D1), D480–D489. https://doi.org/10.1093/nar/gkaa1100

'''