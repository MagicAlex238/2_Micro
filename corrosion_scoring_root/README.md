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