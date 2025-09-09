**README**

This is a scoring system module package made to reuse the scoring system and the global terms so that scoring remains homogeneous and there is no need to recode the terms each time, reducing the risk of inconsistency.

**Project Structure:**
2_Micro
|
corrosion_scoring_root/
├── corrosion_scoring/         # Python package (importable code)  
│   ├── __init__.py            # Package initialization  
│   ├── scoring_system.py      # Functions for scoring found terms  
│   ├── global_terms.py        # Realistic terms from bioinformatics databases (ec_records)  
│   └── term_processor.py      # Functions to avoid duplication and prioritize columns  
│
├── Explanatory/                   
│   └── validation_tools.ipynb # Notebook for validation and term realism checks  
│
├── tests/                      
│   ├── test_scoring_system.py  # Test functions for scoring  
│   └── test_global_terms.py    # Test functions for term lists  
│
[text](.)├── pyproject.toml              # Project settings and version control  
├── README.md                   # This file  
└── .gitignore                  # Files/folders to ignore in version control  

```BASH

cd /home/beatriz/MIC/2_Micro/corrosion_scoring_root/corrosion_scoring
{
    echo "# ========== __init__.py =========="
    cat __init__.py
    #echo -e "\n\n# ========== global_terms.py ==========" # it is too large
    #cat global_terms.py
    echo -e "\n\n# ========== name_utils.py =========="
    cat name_utils.py
    echo -e "\n\n# ========== scoring_system.py =========="
    cat scoring_system.py
    echo -e "\n\n# ========== term_processor.py =========="
    cat term_processor.py
    echo -e "\n\n# ========== utils_ec.py =========="
    cat utils_ec.py
} > /home/beatriz/MIC/2_Micro/corrosion_scoring_root/corrosion_scoring_combined.py
´´´

####_________________-
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