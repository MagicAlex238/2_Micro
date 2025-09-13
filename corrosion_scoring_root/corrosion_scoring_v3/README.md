# Corrosion Scoring System v3.0

A refactored Python module for analyzing corrosion relevance in scientific texts through text mining, scoring, and synergy detection.

## Overview

This system provides a clean, modular architecture for:
- **Text Mining**: Extracting corrosion-related terms and categories from scientific texts
- **Scoring**: Calculating relevance scores based on functional categories, metals, and synergies
- **Synergy Detection**: Identifying co-occurrence patterns that indicate corrosion mechanisms

## Architecture

The system is built with separation of concerns:

```
corrosion_scoring_v3/
├── __init__.py              # Main entry point with backward compatibility
├── config.py                # Configuration settings and parameters
├── exceptions.py            # Custom exception classes
├── text_miner.py           # Feature extraction (no scoring)
├── score_calculator.py     # Score computation from features
├── synergy_detector.py     # Synergy detection and analysis
├── term_processor.py       # Term normalization and matching
├── validators.py           # Input validation utilities
└── global_terms.py         # Term dictionaries and mappings
```

## Key Features

### 🔍 Text Mining
- Functional category detection (16 categories including iron metabolism, biofilm formation, etc.)
- Metal term extraction and standardization
- Mechanism and pathway identification
- Operational factor detection

### 📊 Scoring System
- Weighted scoring for functional categories
- Metal presence scoring
- Priority-based synergy scoring
- Configurable weights and thresholds

### 🔄 Synergy Detection
- Co-occurrence pattern analysis
- Priority synergy pairs (e.g., TOC-chelation, Iron-Organic Acid)
- Multi-pathway synergy detection
- Configurable synergy parameters

## Quick Start

```python
from corrosion_scoring_v3 import calculate_overall_scores, ScoringConfig

# Initialize processors 
processors = {
    'fc_processor': functional_processor,
    'metal_processor': metal_processor,
    'synergy_processor': synergy_processor
}

# Analyze text
text = "mined text with corrosion trigger words is on global terms..."
results = calculate_overall_scores(text, processors)

print(f"Corrosion Relevance Score: {results['corrosion_relevance_score']}")
print(f"Functional Categories: {results['functional_categories_detected']}")
print(f"Detected Metals: {results['consolidated_metals']}")
```

## Configuration

Customize the scoring system:

```python
config = ScoringConfig(
    metal_score_weight=0.5,
    functional_score_weight=1.5,
    synergy_score_weight=2.0,
    min_synergy_categories=2,
    min_synergy_terms=2
)

results = calculate_overall_scores(text, processors, config)
```

## Priority Synergies

The system recognizes key corrosion synergies:

| Synergy Pair | Score | Description |
|-------------|-------|-------------|
| TOC + Chelation | 3.0 | TOC-chelation Synergy |
| Iron + Organic Acid | 2.8 | Acid-enhanced Fe corrosion |
| Biofilm + Chelation | 2.7 | Biofilm-chelate corrosion |
| O2 + Iron | 2.5 | Aerobic Fe corrosion |

## Error Handling

The system includes comprehensive error handling:

```python
from corrosion_scoring_v3.exceptions import ScoringError, TextMiningError

try:
    results = calculate_overall_scores(text, processors)
except TextMiningError as e:
    print(f"Text mining failed: {e}")
except ScoringError as e:
    print(f"Scoring failed: {e}")
```

## Requirements

- Python 3.7+
- Standard library only (no external dependencies)
- existing processor implementations

## Migration from v2.0

The v3.0 system maintains backward compatibility through the `calculate_overall_scores` function.

## Contributing

1. Ensure all tests pass
2. Follow the existing code style
3. Add appropriate error handling
4. Update documentation

## License

MIT

```BASH
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
} > /home/beatriz/MIC/2_Micro/corrosion_scoring_root/corrosion_scoring_v3/v3_combined_files_cs.py
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