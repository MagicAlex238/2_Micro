@startuml
digraph TD {
    A [label="Data Preparation & Filtering"]
    B [label="Feature Selection"]
    C [label="Literature Analysis"]
    D [label="Phylogenetic Analysis"]
    E [label="Functional Analysis"]
    F [label="Visualization"]
    
    A -> B
    B -> C
    C -> D
    D -> E
    E -> F
    
    A -> A1 [label="Data categorization by corrosion degree"]
    A -> A2 [label="Distribution analysis"]
    A -> A3 [label="Genus relevance filtering"]
    A -> A4 [label="Influence frequency calculation"]

    B -> B1 [label="Principal Component Analysis"]
    B -> B2 [label="Random Forest feature importance"]
    B -> B3 [label="Multivariate statistical analysis"]
    B -> B4 [label="Generation of checked_genera, usual_suspects, core_taxa"]

    C -> C1 [label="MIC literature research"]
    C -> C2 [label="Known vs. candidate bacteria differentiation"]
    C -> C3 [label="API calls for taxonomic information"]
    C -> C4 [label="Functional gene analysis (dsrAB, aprAB genes)"]

    D -> D1 [label="Sequence retrieval for bacteria taxa"]
    D -> D2 [label="Alignment, consensus trees"]
    D -> D3 [label="Bootstrapping"]
    D -> D4 [label="Phylogenetic tree visualization"]

    E -> E1 [label="Preparing data for Galaxy Europe"]
    E -> E2 [label="Analyzing picrust2 General Results"]
    E -> E3 [label="Making database, enrichment"]
    E -> E4 [label="Filtering protein-genus markers"]

    F -> F1 [label="Visualizing the results"]
    F -> F2 [label="Protein visualization"]
    F -> F3 [label="Integrated data presentation"]
}
@enduml