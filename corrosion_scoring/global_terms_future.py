#===============================================
# 1. Common metal ions to look for Types of metal ions and compounds
metal_terms = {
    'iron': ['Fe2+', 'Fe3+', 'iron', 'ferrous', 'ferric', 'heme', 'iron-sulfur', 'rust', 'ochre', 'iron oxide', 'iron precipitation', 'siderophore', 'ferritin'],
    'manganese': ['Mn2+', 'Mn3+', 'Mn4+', 'manganese', 'mn', 'manganous', 'manganic', 'manganese oxidation', 'manganese oxide', 'MnO2', 'birnessite', 'pyrolusite'],  # Added Mn3+, Mn4+ for different oxidation states
    'copper': ['Cu+', 'Cu2+', 'copper', 'cupric', 'cuprous', 'copper oxide', 'copper corrosion', 'copper carbonate', 'cuprite', 'malachite'],  # Added copper minerals
    'nickel': ['Ni2+', 'Ni3+', 'nickel', 'nickelous', 'nickel oxidation', 'nickel reduction', 'nickel hydroxide'],  # Added Ni3+ and nickel hydroxide
    'cobalt': ['Co2+', 'Co3+', 'cobalt', 'cobaltous', 'cobalamin', 'vitamin B12', 'cobalt chloride'],  # Added Co3+ and cobalt chloride
    'magnesium': ['Mg2+', 'magnesium', 'magnesium oxide', 'magnesium carbonate', 'magnesium hydroxide'],  # Added Mg compounds
    'calcium': ['Ca2+', 'calcium', 'calcium carbonate', 'calcite', 'calcium precipitation', 'calcium sulfate', 'gypsum'],  # Added calcium sulfate/gypsum
    'Mo': ['Mo', 'Mo4+', 'Mo5+', 'Mo6+', 'molybdenum', 'molybdopterin', 'molybdenum cofactor', 'molybdate', 'molybdic acid'],  # Added Mo oxidation states and molybdate
    'V5+': ['V5+', 'V4+', 'V3+', 'vanadium', 'vanadate', 'vanadyl', 'vanadium oxide'],  # Added multiple vanadium states
    'Al3+': ['Al3+', 'aluminum', 'aluminate', 'aluminum oxide', 'alumina', 'aluminum hydroxide'],  # Added aluminum compounds
    'Cr3+': ['Cr3+', 'Cr6+', 'chromium', 'chromate', 'dichromate', 'chromium oxide', 'hexavalent chromium'],  # Added Cr6+ for hexavalent chromium
    'zinc': ['Zn2+', 'zinc', 'zinc finger', 'zinc oxide', 'zinc sulfide', 'zinc carbonate'],  # Added zinc compounds
    'sodium': ['Na+', 'sodium', 'NaCl', 'sodium transport', 'sodium gradient'],
    'potassium': ['K+', 'potassium', 'KCl', 'potassium transport', 'potassium channel'],
    'selenium': ['selenium', 'Se', 'selenocysteine', 'selenoprotein', 'selenite', 'selenate'],  # Added selenate
    'barium': ['Ba2+', 'barium', 'barium sulfate', 'barite'],
    'strontium': ['Sr2+', 'strontium', 'strontium carbonate', 'strontium sulfate'],  # Added strontium - missing element
    'lead': ['Pb2+', 'Pb4+', 'lead', 'plumbous', 'plumbic', 'lead oxide'],  # Added lead - missing element
    'arsenic': ['As3+', 'As5+', 'arsenic', 'arsenite', 'arsenate', 'arsenic oxidation'],  # Added arsenic - missing element
    'mercury': ['Hg2+', 'Hg+', 'mercury', 'mercuric', 'mercurous', 'mercury methylation'],  # Added mercury - missing element
    'phosphate': ['HPO4-2', 'PO4-3', 'phosphate', 'phosphates', 'orthophosphate'],  # Added orthophosphate
    'nitrate': ['NO3-', 'nitrate', 'nitrates'],
    'nitrite': ['NO2-', 'nitrite', 'nitrites'],
    'chloride': ['Cl-', 'chloride', 'chlorine', 'hypochlorite', 'chlorate'],  # Added hypochlorite and chlorate
    'sulfate': ['SO4-2', 'sulfate', 'sulfates'],
    'sulfide': ['S2-', 'sulfide', 'sulfides', 'H2S', 'hydrogen sulfide', 'pyrite', 'pyrrhotite'],  # Added sulfide minerals
    'thiosulfate': ['S2O3-2', 'thiosulfate'],
    'oxygen': ['O2', 'oxygen', 'oxidase', 'superoxide', 'peroxide', 'hydroxyl radical'],  # Added reactive oxygen species
    'hydrogen': ['H2', 'hydrogen', 'hydrogenase', 'hydrogen uptake', 'hydrogen evolution'],
    'organics': ['methane', 'CH4', 'methane', 'methanogenic', 'methanogenesis', 'formate', 'formate', 'formic acid', 'HCOO-', 'acetate', 'acetate', 'acetic acid', 'CH3COO-', 
                    'propionate','propionate', 'propionic acid', 'butyrate', 'butyrate', 'butyric acid', 'lactate', 'lactic acid', 'mercaptans', 'mercaptan', 'thiol', 'methanethiol',
                    'ethanethiol', 'H2S', 'alcohol', 'ethanol', 'methanol', 'propanol', 'alcohol']
}


# 2. Comprehensive corrosion mechanisms classification Specific biochemical processes causing corrosion
corrosion_mechanisms = {
    'direct_eet': ['cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase', 'electron conduit', 
                    'direct electron transfer', 'extracellular electron transfer'],  # Added extracellular electron transfer
    'indirect_eet': ['shuttle', 'mediator', 'redox mediator', 'electron shuttle', 'flavin', 'quinone', 'humic substance', 'pyocyanin', 'phenazine', 'riboflavin'],  # Added specific mediators
    'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 'lactate metabolism', 
                        'formate production', 'proton generation', 'low pH', 'carbonic acid', 'citric acid', 'gluconic acid'],  # Added more organic acids
    'h2_consumption': ['hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism'],
    'o2_consumption': ['oxidase', 'oxygen reduction', 'aerobic respiration', 'oxygen reduc', 'oxygen consum', 'oxygen scavenging', 'oxygen stress', 'oxidative phosphorylation', 
                        'respiratory burst', 'oxygen depletion'],  # Added respiratory burst and oxygen depletion
    'biofilm_formation': ['polysaccharide', 'adhesin', 'biofilm', 'EPS', 'extracellular polymeric substance', 'curli', 'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 
                            'attachment', 'surface', 'adherence', 'biofilm maturation', 'quorum sensing', 'alginate production', 'cellulose production'],  # Added specific biofilm components
    'sulfur_metabolism': ['sulfate reduc', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB', 'sulfur disproportionation', 'sulfate-reducing bacteria', 
                            'sulfur respiration', 'sulfite reduction', 'elemental sulfur', 'polysulfide'],  # Added sulfite reduction and elemental sulfur
    'metal_transformation': ['iron reduction', 'manganese oxidation', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'metal deposition',    
                                'metal solubilization', 'mineral dissolution', 'mineral precipitation', 'metal mobilization', 'co-precipitation', 'biosorption'],  # Added metal mobilization and biosorption
    'iron_metabolism': ['iron reduc', 'ferric reduc', 'iron oxid', 'ferrous oxid', 'iron uptake', 'iron transport', 'iron storage', 'iron homeostasis', 'siderophore production', 
                        'iron-sulfur cluster', 'ferritin', 'bacterioferritin', 'heme biosynthesis'],  # Added iron-specific proteins
    'metal_chelation': ['siderophore', 'metal binding', 'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelation', 'metal complexation', 'metal sequestration', 
                        'chelate formation', 'metal ligand'],  # Added chelate formation and metal ligand
    'carbon_metabolism': ['carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'carbon flux', 'carbon assimilation', 
                            'carbon catabolite repression', 'carbonation', 'carbonate precipitation'],  # Added carbon catabolite repression and carbonation
    'ph_modulation': ['acid tolerance', 'alkaline tolerance', 'proton pump', 'pH homeostasis', 'pH stress', 'pH regulation', 'acid resistance', 
                        'pH buffering', 'pH gradient', 'acidophile', 'alkaliphile'],  # Added pH buffering and pH gradient
    'nitrogen_metabolism': ['denitrification', 'nitrification', 'nitrate reduction', 'nitrite reduction', 'nitrate respiration', 'nitrite respiration', 
                            'nitrous oxide reduction', 'ammonia oxidation', 'anammox', 'nitrogen fixation', 'ammonification'],  # Added MISSING mechanism
    'manganese_metabolism': ['manganese oxidation', 'manganese reduction', 'Mn-oxide formation', 'Mn-oxide dissolution', 'Mn cycling', 
                            'birnessite formation', 'pyrolusite formation', 'Mn precipitation'],  # Added MISSING mechanism
    'exoelectrogenesis': ['exoelectrogen', 'electrochemically active bacteria', 'EAB', 'extracellular respiration', 'electrode respiration', 
                            'bioelectrosynthesis', 'microbial fuel cell'],  # Added MISSING mechanism - important for direct corrosion
    'enzymatic_corrosion': ['enzyme-mediated corrosion', 'metalloenzyme', 'enzyme-catalyzed oxidation', 'enzyme-catalyzed reduction', 
                            'peroxidase', 'laccase', 'oxidoreductase activity'],  # Added MISSING mechanism
    'dealloying': ['selective corrosion', 'dezincification', 'dealuminification', 'parting', 'selective leaching', 
                    'alloy corrosion', 'preferential dissolution'],  # Added MISSING mechanism
    'galvanic_corrosion': ['bimetallic corrosion', 'galvanic couple', 'differential aeration', 'local cell formation', 
                            'microelectrochemical cells', 'electrolytic corrosion'],  # Added MISSING mechanism
    'chloride_attack': ['chloride-induced corrosion', 'pitting initiation', 'localized corrosion', 'crevice corrosion', 
                        'stress corrosion cracking', 'chloride penetration'],  # Added MISSING mechanism
    'microbe_metal_synergy': ['metal-microbe interaction', 'biomineralization', 'biosorption', 'bioaccumulation', 
                                'metal-binding proteins', 'metallothionein', 'metal-enzyme complex']  # Added MISSING mechanism
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
        'oxidative stress', 'oxygen sensor', 'oxygen tolerance',
        'reactive oxygen species', 'peroxidase', 'oxygen-limited growth'  # Added ROS and peroxidase
    ],
    'nitrogen_metabolism': [
        'denitrification', 'nitrate_reduction', 'nitrite_reduction',
        'nitrogen metabolism', 'nitrate', 'nitrite',
        'ammonia oxidation', 'nitrification', 'nitrogen fixation',
        'ammonification', 'nitrous oxide reduction', 'nitric oxide reduction',
        'dissimilatory nitrate reduction', 'nitrite reductase', 'nitrate reductase', 
        'anammox', 'nitrate respiration', 'nitrite respiration', 'nitric oxide metabolism'  # Added anammox and NO metabolism
    ],
    'manganese_processes': [
        'manganese_reduction', 'mn_redox', 'manganese oxidation',
        'manganese oxide', 'pyrolusite', 'birnessite',
        'manganese cycling', 'manganese mineral', 'manganese transport',
        'Mn-oxide formation', 'Mn-oxide reduction', 'Mn precipitation', 'Mn dissolution'  # Added Mn-oxide processes
    ], 
    'iron_sulfur_redox': [
        'iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 
        'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 
        'thiosulfate', 'sulfur oxidation', 'SRB',
        'Fe-S cluster', 'iron-sulfur cluster', 'ferredoxin',
        'rubredoxin', 'ferritin', 'bacterioferritin', 'PWY-7221', 'HEME-BIOSYNTHESIS-II', 'P125-PWY',
        'iron mobilization', 'iron immobilization', 'iron precipitation'  # Added iron mobilization processes
    ],
    'ocre_formation': [
        'ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 
        'iron oxide deposits', 'iron precipitation', 'rust formation',
        'ferrihydrite', 'goethite', 'magnetite', 'hematite',
        'mineral precipitation', 'iron mineral', 'PWY-7219',
        'biogenic iron oxides', 'stalactite formation', 'ochre mats'  # Added biogenic iron oxides
    ],
    'sulfur_metabolism': [
        'sulfur', 'sulfate', 'sulfide', 'thiosulfate', 'sulfite', 'sulfonate',
        'sulfate reduction', 'sulfur oxidation', 'sulfur respiration',
        'SRB', 'dsrAB', 'APS reductase', 'sulfide:quinone oxidoreductase',
        'sulfur disproportionation', 'dissimilatory sulfate reduction',
        'PWY-6932', 'SO4ASSIM-PWY', 'SULFATE-CYS-PWY',
        'sulfur globules', 'elemental sulfur', 'polysulfide metabolism', 'sulfur granules'  # Added sulfur storage
    ],
    'electron_transfer': [
        'cytochrome', 'electron transport', 'oxidoreductase', 'redox',
        'electron transfer', 'direct EET', 'indirect EET', 'nanowire',
        'conductive pili', 'electron shuttle', 'mtrABC', 'omcS',
        'c-type cytochrome', 'multi-heme cytochrome', 'flavin',
        'electron bridge', 'solid-state electron transfer', 'interfacial electron transfer'  # Added interfacial ET
    ],
    'organic_acid_metabolism': [
        'acetate', 'acetic acid', 'acetyl', 'acetate metabolism', 'acetate production',
        'oxalate', 'oxalic acid', 'oxalate metabolism', 'oxalate production',
        'organic acid', 'fatty acid', 'butyric acid', 'butyrate', 'propionate', 'propionic acid',
        'carboxylic acid', 'lactate', 'lactic acid', 'formate', 'formic acid',
        'citrate', 'citric acid', 'succinate', 'succinic acid', 'fumarate', 'fumaric acid',
        'malate', 'malic acid', 'pyruvate', 'pyruvic acid', 'CENTFERM-PWY', 'FERMENTATION-PWY',
        'GLYCOLYSIS', 'PWY-5100', 'GALACTUROCAT-PWY',
        'gluconic acid', 'tartaric acid', 'caproic acid', 'valeric acid'  # Added more organic acids
    ],
    'metal_organic_interaction': [
        'siderophore', 'metal binding', 'metal chelation', 'iron chelation',
        'iron complex', 'metal transport', 'metallophore', 'metal organic',
        'organometallic', 'iron uptake', 'metal uptake', 'metalloprotein',
        'iron-sulfur cluster', 'metal coordination', 'metal sequestration',
        'ferric reductase', 'ferrous oxidase', 'metal homeostasis',
        'metallophore production', 'metal-ligand complex', 'chelator secretion', 'metal bioavailability'  # Added metal bioavailability
    ],
    'biofilm_formation': [
        'biofilm', 'exopolysaccharide', 'EPS production', 'extracellular matrix',
        'adhesion', 'attachment', 'colonization', 'surface attachment',
        'polysaccharide biosynthesis', 'extracellular polymeric substance',
        'cell aggregation', 'quorum sensing', 'biofilm maturation',
        'biofilm regulation', 'biofilm dispersion', 'cell-cell adhesion',
        'matrix production', 'pellicle', 'floc formation', 'COLANSYN-PWY',
        'EXOPOLYSACC-PWY', 'GLUCOSE1PMETAB-PWY',
        'alginate synthesis', 'cellulose biosynthesis', 'amyloid fiber', 'biofilm architecture'  # Added biofilm components
    ],
    'carbon_metabolism': [
        'carbon fixation', 'carbon utilization', 'carbohydrate metabolism',
        'glycolysis', 'TCA cycle', 'pentose phosphate pathway',
        'gluconeogenesis', 'carbon assimilation', 'carbon flux',
        'Calvin cycle', 'reductive acetyl-CoA pathway', 'carbon monoxide dehydrogenase',
        'carbon catabolite repression', 'carbon storage', 'glycogen metabolism'  # Added carbon storage
    ],
    'ph_modulation': [
        'acid', 'alkaline', 'proton pump', 'pH homeostasis',
        'pH stress', 'acid tolerance', 'alkaline tolerance',
        'proton motive force', 'pH regulation', 'acidic environment',
        'alkaline environment', 'acid resistance', 'proton antiporter',
        'pH gradient maintenance', 'pH sensing', 'pH buffering system'  # Added pH sensing and buffering
    ],
    'temp_response': [
        'heat shock', 'cold shock', 'temperature response',
        'thermophilic', 'psychrophilic', 'mesophilic',
        'thermal adaptation', 'temperature stress', 'heat stress protein',
        'cold stress protein', 'thermal stability', 'thermotolerance'
    ],
    'halogen_related': [
        'halogen', 'chloride', 'bromide', 'iodide', 'fluoride',
        'halide', 'dehalogenation', 'haloperoxidase', 'haloacid',
        'chlorination', 'bromination', 'halomethane', 'haloalkane',
        'organohalide', 'halotolerance', 'salt tolerance',
        'halophilic', 'chloride transport', 'halide channel',
        'hypochlorite', 'chlorate reduction', 'perchlorate reduction'  # Added chlorine oxyanions
    ],
    'methanogenesis': [
        'methanogenesis', 'methanobacterium', 'archaea', 'methane production',
        'methyl-coenzyme M reductase', 'methanogenic', 'coenzyme F420',
        'methyl-H4MPT', 'CO2 reduction', 'acetoclastic methanogenesis',
        'hydrogenotrophic methanogenesis', 'methylotrophic methanogenesis', 'methanotrophy'  # Added methanogenesis types and methanotrophy
    ],
    'silicon_metabolism': [  # Added MISSING pathway
        'silicon uptake', 'silicate transport', 'silicon cycling', 'diatom frustule',
        'silicification', 'biogenic silica', 'silicon precipitation'
    ],
    'phosphorus_metabolism': [  # Added MISSING pathway
        'phosphate uptake', 'phosphate transport', 'polyphosphate', 'phosphorus cycling',
        'phosphate precipitation', 'struvite formation', 'phosphate solubilization'
    ]
} 

# 4. Define organic matter categories - Types of processes organisms perform on organic matter
organic_categories = {
    'degradation': ['degradation', 'breakdown', 'catabolism', 'hydrolysis', 'digestion',
                    'decomposition', 'mineralization', 'dissolution', 'depolymerization', 'lysis'],  # Added depolymerization and lysis
    'synthesis': ['biosynthesis', 'anabolism', 'synthesis', 'polymerization', 'assembly',
                    'construction', 'formation', 'production', 'biogenesis', 'neosynthesis'],  # Added biogenesis
    'transport': ['transport', 'uptake', 'export', 'secretion', 'efflux', 'influx',
                    'extrusion', 'import', 'transmembrane transport', 'translocation', 'diffusion', 'facilitated transport'],  # Added transport types
    'modification': ['modification', 'conversion', 'transformation', 'transmutation',
                    'alteration', 'rearrangement', 'isomerization', 'conjugation', 'methylation', 'hydroxylation'],  # Added specific modifications
    'respiration': ['respiration', 'electron transport chain', 'oxidative phosphorylation',
                    'terminal electron acceptor', 'cytochrome', 'anaerobic respiration', 'substrate-level phosphorylation'],  # Added anaerobic respiration
    'fermentation': ['fermentation', 'anaerobic metabolism', 'mixed acid fermentation',
                    'alcohol fermentation', 'lactic acid fermentation', 'homolactic fermentation', 'heterolactic fermentation'],  # Added fermentation types
    'oxidation': ['oxidation', 'oxidase', 'oxidoreductase', 'dehydrogenase',
                    'hydroxylation', 'oxygenase', 'monooxygenase', 'dioxygenase'],  # Added oxygenase types
    'reduction': ['reduction', 'reductase', 'hydrogenase', 'dehydrogenase',
                    'nitrate reduction', 'sulfate reduction', 'electron donation', 'reduction potential']  # Added electron donation
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
                                    'thiosulfate', 'sulfur oxidation', 'SRB', 'iron-sulfur cluster', 'rubredoxin', 'ferredoxin'], 'score': 1.5},  # Added Fe-S proteins
    'ocre': {'terms': ['ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'ferrihydrite', 'goethite', 'magnetite'], 'score': 1.2},  # Added iron minerals, changed score
    'acid_production': {'terms': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 
                                    'lactate metabolism', 'formate production', 'carbonic acid', 'citric acid'], 'score': 1.3},  # Added more acids, changed score
    'electron transfer & redox': {'terms': ['direct_eet', 'redox', 'electron_transfer', 'omc', 'deet',  'cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 
                                            'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase', 'extracellular electron transfer'], 'score': 1.4},  # Added extracellular ET, changed score
    'biofilm_formation': {'terms': ['biofilm_formation', 'metal_chelation', 'quorum_sensing',  'extracellular_matrix', 'EPS', 'surface_disruption', 'polysaccharide', 'adhesin', 'biofilm', 
                                    'EPS', 'extracellular polymeric substance', 'curli',  'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 'attachment', 'alginate', 'cellulose'], 'score': 1.0},  # Added alginate and cellulose, changed score
    'sulfide_production': {'terms': ['sulfide', 'sulfur_reduction', 'desulfovibrio', 'hydrogen sulfide', 'H2S', 'sulfide precipitation'], 'score': 1.0},  # Added sulfide precipitation, changed score
    'metal binding / chelation': {'terms': ['metal_chelation', 'metal_binding', 'siderophore', 'complexation',  'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelate formation', 'metal ligand'], 'score': 1.2},  # Added chelate terms, changed score
    'nitrogen_reduction': {'terms': ['denitrification', 'nitrate_reduction', 'nitrite_reduction', 'nitrate respiration', 'nitrous oxide reduction', 'nitric oxide metabolism'], 'score': 1.1},  # Added NO metabolism, changed score
    'manganese_reduction': {'terms': ['manganese_reduction', 'mn_redox', 'manganese oxidation', 'birnessite formation', 'pyrolusite formation'], 'score': 1.2},  # Added Mn minerals, changed score
    'methanogenesis': {'terms': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 'methyl-coenzyme M reductase'], 'score': 0.7},  # Changed score
    'fumarate_formation': {'terms': ['fumarate', 'propionibacterium', 'fumarate reduction', 'fumarate respiration'], 'score': 0.5},  # Added fumarate respiration
    'h2_consumption': {'terms': ['h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism'], 'score': 0.6},  # Changed score
    'o2_consumption': {'terms': ['o2_consumption', 'aerobic_respiration', 'oxygen reduction', 'oxygen consumption', 'oxygen depletion', 'oxygen-limited growth'], 'score': 0.5},  # Added O2 depletion, changed score
    'nitrogen_metabolism': {'terms': ['nitrogen_metabolism', 'nitrification', 'ammonification', 'nitrogen fixation', 'anammox'], 'score': 1.1},  # Added MISSING functional category
    'exoelectrogenesis': {'terms': ['exoelectrogen', 'electrochemically active bacteria', 'EAB', 'extracellular respiration', 'electrode respiration'], 'score': 1.3},  # Added MISSING functional category
    'enzymatic_metal_oxid': {'terms': ['metalloenzyme', 'enzyme-catalyzed oxidation', 'peroxidase', 'laccase', 'oxidoreductase activity', 'enzyme-mediated corrosion'], 'score': 1.2},  # Added MISSING functional category
    'metal_precipitation': {'terms': ['metal precipitation', 'biomineralization', 'mineral formation', 'crystal nucleation', 'metal immobilization'], 'score': 1.1},  # Added MISSING functional category
    'ph_buffering': {'terms': ['pH buffering', 'pH homeostasis', 'acid resistance', 'alkaline tolerance', 'pH gradient'], 'score': 0.8},  # Added MISSING functional category
    'chloride_interactions': {'terms': ['chloride attack', 'chloride-induced corrosion', 'pitting initiation', 'chloride penetration', 'halide corrosion'], 'score': 1.2},  # Added MISSING functional category
    'dealloying_mechanisms': {'terms': ['selective corrosion', 'dezincification', 'dealuminification', 'preferential dissolution', 'parting'], 'score': 1.1}  # Added MISSING functional category
}

# 7. CORROSION KEYWORD GROUPS - 
corrosion_keyword_groups = {
    'iron_sulfur_redox': ['iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 
                            'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 
                            'thiosulfate', 'sulfur oxidation', 'SRB', 'pyrite', 'pyrrhotite'],  # Added minerals
    'ocre': ['ocre', 'iron oxide', 'iron deposit', 'metal oxide', 'ochre formation', 
            'iron oxide deposits', 'iron precipitation', 'rust formation', 'ferrihydrite', 'goethite'],  # Added minerals
    'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 
                        'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 
                        'lactate metabolism', 'formate production', 'carbonic acid', 'citric acid'],  # Added more acids
    'electron_transfer': ['direct eet', 'redox', 'electron transfer', 'omc', 'deet', 'cytochrome', 
                            'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 
                            'oxidoreductase', 'redox', 'reductase', 'oxidase', 'extracellular electron transfer'],  # Added extracellular ET
    'biofilm': ['biofilm formation', 'quorum sensing', 'extracellular matrix', 'EPS', 
                'surface disruption', 'polysaccharide', 'adhesin', 'biofilm', 
                'extracellular polymeric substance', 'curli', 'exopolymer', 'attachment', 'colonization', 
                'alginate', 'cellulose', 'biofilm architecture'],  # Added biofilm components
    'sulfide': ['sulfide', 'sulfur reduction', 'desulfovibrio', 'h2s', 'hydrogen sulfide', 'sulfide precipitation'],  # Added sulfide precipitation
    'metal_binding': ['metal chelation', 'metal binding', 'siderophore', 'complexation', 
                        'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelate formation'],  # Added chelate terms
    'nitrogen': ['denitrification', 'nitrate reduction', 'nitrite reduction', 
                'nitrogen metabolism', 'nitrate', 'nitrite', 'nitrification', 'ammonification'],  # Added nitrification and ammonification
    'manganese': ['manganese reduction', 'mn redox', 'manganese oxidation', 'mn reduction', 
                    'birnessite', 'pyrolusite', 'manganese mineral'],  # Added manganese minerals
    'methanogenesis': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 
                        'hydrogenotrophic methanogenesis', 'acetoclastic methanogenesis'],  # Added methanogenesis types
    'fumarate': ['fumarate', 'propionibacterium', 'fumarate reduction', 'fumarate respiration'],  # Added fumarate respiration
    'hydrogen': ['h2 consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 
                'h2', 'h2 oxidation', 'h2ase', 'hydrogen evolution', 'hydrogen metabolism'],  # Added hydrogen evolution
    'oxygen': ['o2 consumption', 'aerobic respiration', 'oxygen reduction', 'oxygen consumption', 
                'oxygen depletion', 'oxygen-limited', 'oxygen stress'],  # Added O2-related terms
    'corrosion_general': ['corrosion', 'deterioration', 'pitting', 
                            'microbially influenced corrosion','MIC', 'metal deterioration', 
                            'galvanic corrosion', 'crevice corrosion', 'stress corrosion'],  # Added more corrosion types
    'exoelectrogen': ['exoelectrogen', 'electrochemically active', 'EAB', 'extracellular respiration', 
                        'electrode respiration', 'bioelectrosynthesis'],  # Added MISSING keyword group
    'enzymatic_corrosion': ['enzyme-mediated corrosion', 'metalloenzyme', 'laccase', 'peroxidase', 
                            'enzyme-catalyzed', 'oxidoreductase'],  # Added MISSING keyword group
    'metal_precipitation': ['metal precipitation', 'biomineralization', 'mineral formation', 
                            'mineral precipitation', 'crystal nucleation'],  # Added MISSING keyword group
    'chloride_attack': ['chloride attack', 'chloride corrosion', 'pitting initiation', 'chloride penetration', 
                        'hypochlorite', 'halide corrosion'],  # Added MISSING keyword group
    'ph_regulation': ['pH buffering', 'pH homeostasis', 'acid resistance', 'pH gradient', 
                        'proton pump', 'pH regulation'],  # Added MISSING keyword group
    'dealloying': ['dealloying', 'selective corrosion', 'dezincification', 'preferential dissolution', 
                    'parting', 'alloy corrosion']  # Added MISSING keyword group
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

## Maybe to update 
global_terms = { # it is possible to change the whole conseption to avoid overlapping and to many dictionaries 


    "biochemical_mechanisms": {
        "sulfate_reduction": ["sulfate reduction", "SRB", "sulfidogenesis"],
        "iron_oxidation": ["iron oxidizing", "Fe(II) oxidation", "ferric iron generation"],
        "nitrate_reduction": ["denitrification", "nitrate reduction"],
        "methanogenesis": ["methane production", "methanogens"],
        "manganese_reduction": ["Mn reduction", "manganese respiration"],
        # Clear biological causatives
    },

    "functional_keywords": {
        "cytochromes": ["cytochrome", "c-type cytochrome"],
        "hydrogenases": ["hydrogenase", "NiFe hydrogenase"],
        "oxidoreductases": ["oxidoreductase", "dehydrogenase"],
        "monooxygenases": ["monooxygenase"],
        "peroxidases": ["peroxidase", "catalase"],
    },

    "corrosion_pathways": {
        "sulfur_cycle": ["sulfur cycle", "sulfide oxidation", "sulfur metabolism"],
        "nitrogen_cycle": ["nitrate metabolism", "denitrification pathway"],
        "iron_cycle": ["iron metabolism", "iron cycling"],
        "methanogenesis": ["methanogenesis pathway"],
    },

    "organismal_processes": {
        "biofilm_formation": ["biofilm", "sessile community"],
        "quorum_sensing": ["quorum sensing", "AHL signaling"],
        "motility": ["flagella", "motility"],
        "EPS_production": ["extracellular polymeric substances", "EPS"]
    },

    "corrosion_modes": {
        "pitting": ["pitting", "localized corrosion"],
        "underdeposit": ["underdeposit corrosion", "crevice corrosion"],
        "dealloying": ["dealloying", "selective leaching"]
    },

    "synergistic_interactions": {
        "Fe_S": ["Fe-S interaction", "iron-sulfur"],
        "Mn_N": ["Mn-N coupling", "manganese-nitrate"],
        "Zn_S": ["zinc-sulfide interaction"]
        # Only keep pairs with biological evidence or presence in dataset
    }
}

# Comments:
# - Split "corrosion_mechanisms" into "biochemical_mechanisms" and "corrosion_modes"
# - Renamed "functional_categories" to "functional_keywords" for clarity
# - Merged "organic_categories" into "organismal_processes"
# - Removed "corrosion_keyword_groups" (too broad, overlaps)
# - Curated only biologically relevant metals based on measured elements
# - Added dealloying under corrosion modes (kept niche but valid if detected)
# - All entries designed to reduce overlap and keep categories orthogonal
