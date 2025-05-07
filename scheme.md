```mermaid
graph TD;
    A[Data Preparation & Filtering] -->|↓| B[Feature Selection];
    B -->|↓| C[Literature Analysis];
    C -->|↓| D[Phylogenetic Analysis];
    D -->|↓| E[Functional Analysis];
    E -->|↓| F[Visualization];

    A --> A1[Data categorization by corrosion degree];
    A --> A2[Distribution analysis];
    A --> A3[Genus relevance filtering];
    A --> A4[Influence frequency calculation];

    B --> B1[Principal Component Analysis];
    B --> B2[Random Forest feature importance];
    B --> B3[Multivariate statistical analysis];
    B --> B4[Generation of checked_genera, usual_suspects, core_taxa];

    C --> C1[MIC literature research];
    C --> C2[Known vs. candidate bacteria differentiation];
    C --> C3[API calls for taxonomic information];
    C --> C4[Functional gene analysis (dsrAB, aprAB genes)];

    D --> D1[Sequence retrieval for bacteria taxa];
    D --> D2[Alignment, consensus trees];
    D --> D3[Bootstrapping];
    D --> D4[Phylogenetic tree visualization];

    E --> E1[Preparing data for Galaxy Europe];
    E --> E2[Analyzing picrust2 General Results];
    E --> E3[Making database, enrichment];
    E --> E4[Filtering protein-genus markers];

    F --> F1[Visualizing the results];
    F --> F2[Protein visualization];
    F --> F3[Integrated data presentation];
