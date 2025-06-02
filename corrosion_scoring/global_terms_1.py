#===============================================
# 1. Common metal ions to look for Types of metal ions and compounds
metal_terms = {
    'iron': ['Fe2+', 'Fe3+', 'iron', 'ferrous', 'ferric', 'heme', 'iron-sulfur', 'rust', 'ochre', 'iron oxide', 'iron precipitation', 'siderophore', 'ferritin'],
    'manganese': ['Mn2+', 'manganese', 'mn', 'manganous', 'manganic', 'manganese oxidation', 'manganese oxide', 'MnO2'],
    'copper': ['Cu+', 'Cu2+', 'copper', 'cupric', 'cuprous', 'copper oxide', 'copper corrosion'],
    'nickel': ['Ni2+', 'nickel', 'nickelous', 'nickel oxidation', 'nickel reduction'],
    'cobalt': ['Co2+', 'cobalt', 'cobaltous', 'cobalamin', 'vitamin B12'],
    'magnesium': ['Mg2+', 'magnesium', 'magnesium oxide'],
    'calcium': ['Ca2+', 'calcium', 'calcium carbonate', 'calcite', 'calcium precipitation'],
    'Mo': ['Mo', 'molybdenum', 'molybdopterin', 'molybdenum cofactor'],
    'V5+': ['V5+', 'vanadium', 'vanadate', 'vanadyl'],
    'Al3+': ['Al3+', 'aluminum', 'aluminate', 'aluminum oxide'],
    'Cr3+': ['Cr3+', 'chromium', 'chromate', 'dichromate', 'chromium oxide'],
    'zinc': ['Zn2+', 'zinc', 'zinc finger', 'zinc oxide'],
    'sodium': ['Na+', 'sodium', 'NaCl', 'sodium transport', 'sodium gradient'],
    'potassium': ['K+', 'potassium', 'KCl', 'potassium transport', 'potassium channel'],
    'selenium': ['selenium', 'Se', 'selenocysteine', 'selenoprotein', 'selenite'],
    'barium': ['Ba2+', 'barium', 'barium sulfate', 'barite'],
    'phosphate': ['HPO4-2', 'PO4-3', 'phosphate', 'phosphates'],
    'nitrate': ['NO3-', 'nitrate', 'nitrates'],
    'nitrite': ['NO2-', 'nitrite', 'nitrites'],
    'chloride': ['Cl-', 'chloride', 'chlorine', 'hypochlorite', 'chlorate'],
    'sulfate': ['SO4-2', 'sulfate', 'sulfates'],
    'sulfide': ['S', 'sulfide', 'sulfides', 'H2S', 'hydrogen sulfide', 'pyrite', 'pyrrhotite'],
    'thiosulfate': ['S2O3-2', 'thiosulfate'],
    'oxygen': ['O2', 'oxygen', 'oxidase', 'superoxide', 'peroxide', 'hydroxyl radical'],
    'hydrogen': ['H2', 'hydrogen', 'hydrogenase', 'hydrogen uptake', 'hydrogen evolution'],
    'organics': ['methane', 'CH4', 'methane', 'methanogenic', 'methanogenesis', 'formate','formate', 'formic acid', 'HCOO-', 'acetate', 'acetate', 'acetic acid', 'CH3COO-', 
                    'propionate','propionate', 'propionic acid', 'butyrate','butyrate', 'butyric acid','lactate', 'lactic acid', 'mercaptans', 'mercaptan', 'thiol', 'methanethiol',
                    'ethanethiol', 'H2S', 'alcohol', 'ethanol', 'methanol', 'propanol', 'alcohol']
}

# 2. Comprehensive corrosion mechanisms classification Specific biochemical processes causing corrosion
corrosion_mechanisms = {
    'direct_eet': ['cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase', 'electron conduit',  'direct electron transfer', 'extracellular electron transfer'],
    'indirect_eet': ['shuttle', 'mediator', 'redox mediator', 'electron shuttle', 'flavin', 'quinone', 'humic substance'],
    'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 'lactate metabolism', 'formate production', 'proton generation', 'low pH', 'carbonic acid', 'citric acid', 'gluconic acid'],
    'h2_consumption': ['hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism'],
    'o2_consumption': ['oxidase', 'oxygen reduction', 'aerobic respiration', 'oxygen reduc', 'oxygen consum',  'oxygen scavenging', 'oxygen stress', 'oxidative phosphorylation', 'respiratory burst', 'oxygen depletion'],
    'biofilm_formation': ['polysaccharide', 'adhesin', 'biofilm', 'EPS', 'extracellular polymeric substance', 'curli', 'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 'attachment', 'surface', 'adherence', 'biofilm maturation'],
    'nitrogen_metabolism': ['denitrification', 'nitrification', 'nitrate reduction', 'nitrite reduction', 'nitrate respiration', 'nitrite respiration', 'nitrous oxide reduction', 'ammonia oxidation', 'anammox', 'nitrogen fixation', 'ammonification'],  
    'manganese_metabolism': ['manganese oxidation', 'manganese reduction', 'Mn-oxide formation', 'Mn-oxide dissolution', 'Mn cycling',  'birnessite formation', 'pyrolusite formation', 'Mn precipitation'],  
    'sulfur_metabolism': ['sulfate reduc', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB',     'sulfur disproportionation', 'sulfate-reducing bacteria', 'sulfur respiration'],
    'metal_transformation': ['iron reduction', 'manganese oxidation', 'metal oxide', 'ochre formation',  'iron oxide deposits', 'iron precipitation', 'rust formation', 'metal deposition', 'metal solubilization', 'mineral dissolution', 'mineral precipitation'],
    'iron_metabolism': ['iron reduc', 'ferric reduc', 'iron oxid', 'ferrous oxid', 'iron uptake', 'iron transport', 'iron storage', 'iron homeostasis', 'ferritin', 'bacterioferritin', 'heme biosynthesis'],
    'metal_chelation': ['siderophore', 'metal binding', 'chelator', 'metallophore', 'iron complex',   'metal transport', 'chelation', 'metal complexation', 'metal sequestration',  'chelate formation', 'metal ligand'], 
    'carbon_metabolism': ['carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis',   'TCA cycle', 'carbon flux', 'carbon assimilation'],
    'ph_modulation': ['acid tolerance', 'alkaline tolerance', 'proton pump', 'pH homeostasis', 'pH stress',  'pH regulation', 'acid resistance']
}

# 3. Expanded pathway categories Metabolic pathways relevant to corrosion
pathway_categories = {
    'hydrogen_metabolism': [
        'h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 
        'h2', 'H2 oxidation', 'H2ase', 'hydrogen production',
        'FeFe-hydrogenase', 'NiFe-hydrogenase', 'hydrogen evolution',
        'hydrogen cycling', 'H2 sensing', 'proton reduction'
    ],
    'oxygen_metabolism': [
        'o2_consumption', 'aerobic_respiration', 'oxygen reduction', 
        'oxygen consumption', 'cytochrome oxidase', 'terminal oxidase',
        'oxygen reductase', 'superoxide dismutase', 'catalase',
        'oxidative stress', 'oxygen sensor', 'oxygen tolerance'
    ],
    'nitrogen_metabolism': [
        'denitrification', 'nitrate_reduction', 'nitrite_reduction',
        'nitrogen metabolism', 'nitrate', 'nitrite',
        'ammonia oxidation', 'nitrification', 'nitrogen fixation',
        'ammonification', 'nitrous oxide reduction', 'nitric oxide reduction',
        'dissimilatory nitrate reduction', 'nitrite reductase', 'nitrate reductase'
    ],
    'manganese_processes': [
        'manganese_reduction', 'mn_redox', 'manganese oxidation',
        'manganese oxide', 'pyrolusite', 'birnessite',
        'manganese cycling', 'manganese mineral', 'manganese transport'
    ], 
    'iron_sulfur_redox': [
        'iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 
        'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 
        'thiosulfate', 'sulfur oxidation', 'SRB',
        'Fe-S cluster', 'iron-sulfur cluster', 'ferredoxin',
        'rubredoxin', 'ferritin', 'bacterioferritin',   'PWY-7221',  # Iron reduction
            'HEME-BIOSYNTHESIS-II',  # Iron-containing compounds
            'P125-PWY'  # Metal resistance
    ],
    'ocre_formation': [
        'ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 
        'iron oxide deposits', 'iron precipitation', 'rust formation',
        'ferrihydrite', 'goethite', 'magnetite', 'hematite',
        'mineral precipitation', 'iron mineral', 'PWY-7219',  # Iron oxidation
    ],
    'sulfur_metabolism': [
        'sulfur', 'sulfate', 'sulfide', 'thiosulfate', 'sulfite', 'sulfonate',
        'sulfate reduction', 'sulfur oxidation', 'sulfur respiration',
        'SRB', 'dsrAB', 'APS reductase', 'sulfide:quinone oxidoreductase',
        'sulfur disproportionation', 'dissimilatory sulfate reduction',
        'PWY-6932',  # Sulfate reduction
        'SO4ASSIM-PWY',  # Sulfate assimilation
        'SULFATE-CYS-PWY'  # Sulfate to cysteine
    ],
    'electron_transfer': [
        'cytochrome', 'electron transport', 'oxidoreductase', 'redox',
        'electron transfer', 'direct EET', 'indirect EET', 'nanowire',
        'conductive pili', 'electron shuttle', 'mtrABC', 'omcS',
        'c-type cytochrome', 'multi-heme cytochrome', 'flavin'
    ],
    'organic_acid_metabolism': [
        'acetate', 'acetic acid', 'acetyl', 'acetate metabolism', 'acetate production',
        'oxalate', 'oxalic acid', 'oxalate metabolism', 'oxalate production',
        'organic acid', 'fatty acid', 'butyric acid', 'butyrate', 'propionate', 'propionic acid',
        'carboxylic acid', 'lactate', 'lactic acid', 'formate', 'formic acid',
        'citrate', 'citric acid', 'succinate', 'succinic acid', 'fumarate', 'fumaric acid',
        'malate', 'malic acid', 'pyruvate', 'pyruvic acid',  'CENTFERM-PWY',  # Central fermentation pathways
        'FERMENTATION-PWY',  # Mixed acid fermentation
        'GLYCOLYSIS',  # Glucose fermentation
        'PWY-5100',  # Pyruvate fermentation
        'GALACTUROCAT-PWY'  # Galacturonate degradation
    ],
    'metal_organic_interaction': [
        'siderophore', 'metal binding', 'metal chelation', 'iron chelation', 'iron complex', 'metal transport', 'metallophore', 'metal organic',
        'organometallic', 'iron uptake', 'metal uptake', 'metalloprotein',  'iron-sulfur cluster', 'metal coordination', 'metal sequestration',
        'ferric reductase', 'ferrous oxidase', 'metal homeostasis', 'metal bioavailability', 'metallophore production', 'metal-ligand complex', 'chelator secretion',    ],
    
    'biofilm_formation': [
        'biofilm', 'exopolysaccharide', 'EPS production', 'extracellular matrix', 'adhesion', 'attachment', 'colonization', 'surface attachment',
        'polysaccharide biosynthesis', 'extracellular polymeric substance',  'cell aggregation', 'quorum sensing', 'biofilm maturation',
        'biofilm regulation', 'biofilm dispersion', 'cell-cell adhesion',  'matrix production', 'pellicle', 'floc formation', 'COLANSYN-PWY',  # Colanic acid (biofilm)
            'EXOPOLYSACC-PWY',  # Exopolysaccharide
            'GLUCOSE1PMETAB-PWY'  # UDP-glucose synthesis
    ],
    'carbon_metabolism': [
        'carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'pentose phosphate pathway', 'gluconeogenesis', 'carbon assimilation', 'carbon flux',
        'Calvin cycle', 'reductive acetyl-CoA pathway', 'carbon monoxide dehydrogenase'
    ],
    'ph_modulation': [
        'acid', 'alkaline', 'proton pump', 'pH homeostasis',   'pH stress', 'acid tolerance', 'alkaline tolerance',  'proton motive force', 'pH regulation', 'acidic environment',
        'alkaline environment', 'acid resistance', 'proton antiporter' ],
    'temp_response': [
        'heat shock', 'cold shock', 'temperature response', 'thermophilic', 'psychrophilic', 'mesophilic', 'thermal adaptation', 'temperature stress', 'heat stress protein',
        'cold stress protein', 'thermal stability', 'thermotolerance'
    ],
    'halogen_related': [
        'halogen', 'chloride', 'bromide', 'iodide', 'fluoride',
        'halide', 'dehalogenation', 'haloperoxidase', 'haloacid',
        'chlorination', 'bromination', 'halomethane', 'haloalkane',
        'organohalide', 'halotolerance', 'salt tolerance',
        'halophilic', 'chloride transport', 'halide channel'
    ],
    'methanogenesis': [
        'methanogenesis', 'methanobacterium', 'archaea', 'methane production',  'methyl-coenzyme M reductase', 'methanogenic', 'coenzyme F420',
        'methyl-H4MPT', 'CO2 reduction', 'acetoclastic methanogenesis',  'hydrogenotrophic methanogenesis', 'methylotrophic methanogenesis', 'methanotrophy' 
    ],
} 
# 4. Define organic matter categories - Types of processes organisms perform on organic matter
organic_categories = {
    'degradation': ['degradation', 'breakdown', 'catabolism', 'hydrolysis', 'digestion',
                    'decomposition', 'mineralization', 'dissolution'],
    'synthesis': ['biosynthesis', 'anabolism', 'synthesis', 'polymerization', 'assembly',
                    'construction', 'formation', 'production'],
    'transport': ['transport', 'uptake', 'export', 'secretion', 'efflux', 'influx',
                    'extrusion', 'import', 'transmembrane transport', 'translocation'],
    'modification': ['modification', 'conversion', 'transformation', 'transmutation',
                    'alteration', 'rearrangement', 'isomerization', 'conjugation'],
    'respiration': ['respiration', 'electron transport chain', 'oxidative phosphorylation',
                    'terminal electron acceptor', 'cytochrome'],
    'fermentation': ['fermentation', 'anaerobic metabolism', 'mixed acid fermentation',
                    'alcohol fermentation', 'lactic acid fermentation'],
    'oxidation': ['oxidation', 'oxidase', 'oxidoreductase', 'dehydrogenase',
                    'hydroxylation', 'oxygenase'],
    'reduction': ['reduction', 'reductase', 'hydrogenase', 'dehydrogenase',
                    'nitrate reduction', 'sulfate reduction']
}
# 5. Corrosion-specific metal combinations
corrosion_synergies = {
    'Fe-S': ['iron sulfur', 'Fe-S', 'iron sulfide', 'FeS', 'Fe-S cluster', 'pyrite', 'pyrrhotite'],  # Added sulfide minerals
    'Fe-Cl': ['iron chloride', 'FeCl', 'iron halide', 'ferric chloride', 'ferrous chloride'],  # Added ferrous chloride
    'Fe-C': ['iron carbon', 'FeC', 'iron carbonate', 'siderite', 'carbonate corrosion'],  # Added carbonate corrosion
    'Cu-Fe': ['copper iron', 'Cu-Fe', 'bimetallic', 'galvanic couple', 'bronze', 'brass'],  # Added copper alloys
    'Mn-Fe': ['manganese iron', 'Mn-Fe', 'iron manganese oxide', 'steel manganese'],  # Added steel manganese
    'Ni-Fe': ['nickel iron', 'Ni-Fe', 'stainless steel', 'alloy corrosion'],  # Added MISSING synergy
    'Cr-Fe': ['chromium iron', 'Cr-Fe', 'stainless steel', 'chromium passivation'],  # Added MISSING synergy
    'Al-Cu': ['aluminum copper', 'Al-Cu', 'aluminum brass', 'galvanic corrosion'],  # Added MISSING synergy
    'Zn-Fe': ['zinc iron', 'Zn-Fe', 'galvanized steel', 'sacrificial anode']  # Added MISSING synergy
}
# 6. FUNCTIONAL CATEGORIES - 
functional_categories = {
    'iron/sulfur_redox': {'terms': ['iron_metabolism', 'sulfur_metabolism', 'iron_oxidation', 'iron_reduction',  'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 
                                    'thiosulfate', 'sulfur oxidation', 'SRB'],'score': 1.5},
    'ocre': {'terms': ['ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation'],'score': 1.5},   
    'acid_production': {'terms': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 
                                    'lactate metabolism', 'formate production'],'score': 1.5},
    'electron transfer & redox': {'terms': ['direct_eet', 'redox', 'electron_transfer', 'omc', 'deet',  'cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 
                                            'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase'],'score': 1.1},
    'biofilm_formation': {'terms': ['biofilm_formation', 'metal_chelation', 'quorum_sensing',  'extracellular_matrix', 'EPS', 'surface_disruption', 'polysaccharide', 'adhesin', 'biofilm', 
                                    'EPS', 'extracellular polymeric substance', 'curli',  'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 'attachment'], 'score': 1.2},
    'sulfide_production': {'terms': ['sulfide', 'sulfur_reduction', 'desulfovibrio'], 'score': 1.3},
    'metal binding / chelation': {'terms': ['metal_chelation', 'metal_binding', 'siderophore', 'complexation',  'chelator', 'metallophore', 'iron complex', 'metal transport'], 'score': 1.0},
    'nitrogen_reduction': {'terms': ['denitrification', 'nitrate_reduction', 'nitrite_reduction'], 'score': 1.0},
    'manganese_reduction': {'terms': ['manganese_reduction', 'mn_redox', 'manganese oxidation'], 'score': 1.0},
    'methanogenesis': {'terms': ['methanogenesis', 'methanobacterium', 'archaea'], 'score': 0.6},
    'fumarate_formation': {'terms': ['fumarate', 'propionibacterium'], 'score': 0.5},
    'h2_consumption': {'terms': ['h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase'], 'score': 0.5},
    'o2_consumption': {'terms': ['o2_consumption', 'aerobic_respiration', 'oxygen reduction', 'oxygen consumption'],'score': 0.6}
}

# 7. CORROSION KEYWORD GROUPS - 
corrosion_keyword_groups = {
    'iron_sulfur_redox': ['iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 
                            'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 
                            'thiosulfate', 'sulfur oxidation', 'SRB'],
    'ocre': ['ocre', 'iron oxide', 'iron deposit', 'metal oxide', 'ochre formation', 
            'iron oxide deposits', 'iron precipitation', 'rust formation'],
    'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 
                        'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 
                        'lactate metabolism', 'formate production'],
    'electron_transfer': ['direct eet', 'redox', 'electron transfer', 'omc', 'deet', 'cytochrome', 
                            'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 
                            'oxidoreductase', 'redox', 'reductase', 'oxidase'],
    'biofilm': ['biofilm formation', 'quorum sensing', 'extracellular matrix', 'EPS', 
                'surface disruption', 'polysaccharide', 'adhesin', 'biofilm', 
                'extracellular polymeric substance', 'curli', 'exopolymer', 'attachment', 'colonization'],
    'sulfide': ['sulfide', 'sulfur reduction', 'desulfovibrio', 'h2s', 'hydrogen sulfide'],
    'metal_binding': ['metal chelation', 'metal binding', 'siderophore', 'complexation', 
                        'chelator', 'metallophore', 'iron complex', 'metal transport'],
    'nitrogen': ['denitrification', 'nitrate reduction', 'nitrite reduction', 'nitrogen metabolism', 'nitrate', 'nitrite', 'nitrification', 'ammonification'],
    'manganese': ['manganese reduction', 'mn redox', 'manganese oxidation', 'mn reduction'],
    'methanogenesis': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production'],
    'fumarate': ['fumarate', 'propionibacterium', 'fumarate reduction'],
    'hydrogen': ['h2 consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 
                'h2', 'h2 oxidation', 'h2ase'],
    'oxygen': ['o2 consumption', 'aerobic respiration', 'oxygen reduction', 'oxygen consumption'],
    'corrosion_general': ['corrosion', 'deterioration', 'pitting', 
                            'microbially influenced corrosion','MIC', 'metal deterioration']
}

# Standard mapping: lower-case keys for matching, with standard symbols as values
metal_mapping = {
    'iron': 'Fe',
    'fe': 'Fe',
    'ferrous': 'Fe',
    'ferric': 'Fe',
    'heme': 'Fe',
    'iron-sulfur': 'Fe',
    'fe2+': 'Fe',
    'fe3+': 'Fe',
    'manganese': 'Mn',
    'mn': 'Mn',
    'manganous': 'Mn',
    'manganic': 'Mn',
    'manganese oxidation': 'Mn',
    'metal oxide': 'Mn',
    'copper': 'Cu',
    'cu+': 'Cu',
    'cu2+': 'Cu',
    'nickel': 'Ni',
    'ni2+': 'Ni',
    'cobalt': 'Co',
    'co2+': 'Co',
    'zinc': 'Zn',
    'zn2+': 'Zn',
    'calcium': 'Ca',
    'ca2+': 'Ca',
    'molybdenum': 'Mo',
    'mo': 'Mo',
    'vanadium': 'V5+',
    'v5+': 'V5+',
    'aluminum': 'Al3+',
    'al3+': 'Al3+',
    'chromium': 'Cr3+',
    'cr3+': 'Cr3+',
    'sodium': 'Na',
    'na+': 'Na',
    'nacl': 'Na',
    'potassium': 'K',
    'k+': 'K',
    'kcl': 'K',
    'selenium': 'Se',
    'se': 'Se',
    'barium': 'Ba2+',
    'ba2+': 'Ba2+',
    'sulfate': 'S',
    'sulfide': 'S',
    'thiosulfate': 'S',
    's-s': 'S',
    'sulfur': 'S',
    'sulfur oxidation': 'S',
    'srb': 'S',
    'hydrogen': 'H',
    'h2': 'H',
    'h2o': 'H',
    'h2s': 'H',
    'phosfate': 'po4-3',
    'nitrate': 'NO3-',
    'nitrite': 'NO2',
    'chloride': 'Cl-'
 }
