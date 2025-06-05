# COMPLETE GLOBAL TERMS 

#===============================================
# 1. Complete metal_terms 
metal_terms = {
'iron': ['Fe2+', 'Fe3+', 'iron', 'ferrous', 'ferric', 'heme', 'iron-sulfur', 'rust', 'ochre', 'iron oxide', 'iron precipitation', 'siderophore', 'ferritin', 'ferredoxin', 'rubredoxin', 'iron-sulfur cluster', 'iron acetate', 'iron oxalate', 'iron carbonate', 'NiFe-hydrogenase'],
'manganese': ['Mn2+', 'Mn3+', 'Mn4+', 'manganese', 'mn', 'manganous', 'manganic', 'manganese oxidation', 'manganese oxide', 'MnO2', 'birnessite', 'pyrolusite'],
'copper': ['Cu+', 'Cu2+', 'copper', 'cupric', 'cuprous', 'copper oxide', 'copper corrosion', 'copper carbonate', 'cuprite', 'malachite'],
'nickel': ['Ni2+', 'Ni3+', 'nickel', 'nickelous', 'nickel oxidation', 'nickel reduction', 'nickel hydroxide', 'NiFe-hydrogenase'],
'cobalt': ['Co2+', 'Co3+', 'cobalt', 'cobaltous', 'cobalamin', 'vitamin B12', 'cobalt chloride'],
'magnesium': ['Mg2+', 'magnesium', 'magnesium oxide', 'magnesium carbonate', 'magnesium hydroxide'],
'calcium': ['Ca2+', 'calcium', 'calcium carbonate', 'calcite', 'calcium precipitation', 'calcium sulfate', 'gypsum'],
'Mo': ['Mo', 'Mo4+', 'Mo5+', 'Mo6+', 'molybdenum', 'molybdopterin', 'molybdenum cofactor', 'molybdate', 'molybdic acid'],
'V5+': ['V5+', 'V4+', 'V3+', 'vanadium', 'vanadate', 'vanadyl', 'vanadium oxide'],
'Al3+': ['Al3+', 'aluminum', 'aluminate', 'aluminum oxide', 'alumina', 'aluminum hydroxide'],
'Cr3+': ['Cr3+', 'Cr6+', 'chromium', 'chromate', 'dichromate', 'chromium oxide', 'hexavalent chromium'],
'zinc': ['Zn2+', 'zinc', 'zinc finger', 'zinc oxide', 'zinc sulfide', 'zinc carbonate'],
'sodium': ['Na+', 'sodium', 'NaCl', 'sodium transport', 'sodium gradient'],
'potassium': ['K+', 'potassium', 'KCl', 'potassium transport', 'potassium channel'],
'selenium': ['selenium', 'Se', 'selenocysteine', 'selenoprotein', 'selenite', 'selenate'],
'barium': ['Ba2+', 'barium', 'barium sulfate', 'barite'],
'strontium': ['Sr2+', 'strontium', 'strontium carbonate', 'strontium sulfate'], 
'lead': ['Pb2+', 'Pb4+', 'lead', 'plumbous', 'plumbic', 'lead oxide'], 
'arsenic': ['As3+', 'As5+', 'arsenic', 'arsenite', 'arsenate', 'arsenic oxidation'], 
'mercury': ['Hg2+', 'Hg+', 'mercury', 'mercuric', 'mercurous', 'mercury methylation'], 
'phosphate': ['HPO4-2', 'PO4-3', 'phosphate', 'phosphates', 'orthophosphate'],
'nitrate': ['NO3-', 'nitrate', 'nitrates'],
'nitrite': ['NO2-', 'nitrite', 'nitrites'],
'chloride': ['Cl-', 'chloride', 'chlorine', 'hypochlorite', 'chlorate'],
'sulfate': ['SO4-2', 'sulfate', 'sulfates'],
'sulfide': ['S2-', 'sulfide', 'sulfides', 'H2S', 'hydrogen sulfide', 'pyrite', 'pyrrhotite'],
'thiosulfate': ['S2O3-2', 'thiosulfate'],
'oxygen': ['O2', 'oxygen', 'oxidase', 'superoxide', 'peroxide', 'hydroxyl radical'],
'hydrogen': ['H2', 'hydrogen', 'hydrogenase', 'hydrogen uptake', 'hydrogen evolution'],
'organics': ['methane', 'CH4', 'methane', 'methanogenic', 'methanogenesis', 'formate','formate', 'formic acid', 'HCOO-', 'acetate', 'acetate', 'acetic acid', 'CH3COO-', 
'propionate','propionate', 'propionic acid', 'butyrate','butyrate', 'butyric acid','lactate', 'lactic acid', 'mercaptans', 'mercaptan', 'thiol', 'methanethiol',
'ethanethiol', 'H2S', 'alcohol', 'ethanol', 'methanol', 'propanol', 'alcohol']
}

#===============================================
# 2. Complete corrosion_mechanisms 
corrosion_mechanisms = { 
'o2_consumption': ['oxidase', 'oxygen reduction', 'aerobic respiration', 'oxygen consumtion', 'oxygen scavenging', 'oxygen stress', 'oxidative phosphorylation', 'respiratory burst', 'oxygen depletion'], 
'nitrogen_metabolism': ['denitrification', 'nitrification', 'nitrate reduction', 'nitrite reduction', 'nitrate respiration', 'nitrite respiration', 'nitrous oxide reduction', 'ammonia oxidation', 'anammox', 'nitrogen fixation', 'ammonification'], 
'iron_metabolism': ['iron reducction', 'ferric reducction', 'iron oxidation', 'ferrous oxidation', 'iron uptake', 'iron transport', 'iron storage', 'iron homeostasis', 'siderophore production', 'iron-sulfur cluster', 'ferritin', 'bacterioferritin', 'heme biosynthesis'], 
'sulfur_metabolism': ['sulfate reduc', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB', 'sulfur disproportionation', 'sulfate-reducing bacteria', 'sulfur respiration', 'sulfite reduction', 'elemental sulfur', 'polysulfide'], 
'h2_consumption': ['hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism'],
'direct_eet': ['cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase', 'electron conduit', 'direct electron transfer', 'extracellular electron transfer'], 
'carbon_metabolism': ['carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'carbon flux', 'carbon assimilation', 'carbon catabolite repression', 'carbonation', 'carbonate precipitation'], 
'indirect_eet': ['shuttle', 'mediator', 'redox mediator', 'electron shuttle', 'flavin', 'quinone', 'humic substance', 'pyocyanin', 'phenazine', 'riboflavin'], 
'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 'lactate metabolism', 'formate production', 'proton generation', 'low pH', 'carbonic acid', 'citric acid', 'gluconic acid'], 
'metal_chelation': ['siderophore', 'metal binding', 'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelation', 'metal complexation', 'metal sequestration', 'chelate formation', 'metal ligand'], 
'biofilm_formation': ['polysaccharide', 'adhesin', 'biofilm', 'EPS', 'extracellular polymeric substance', 'curli', 'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 'attachment', 'surface', 'adherence', 'biofilm maturation', 'quorum sensing', 'alginate production', 'cellulose production'], 

'manganese_metabolism': ['manganese oxidation', 'manganese reduction', 'Mn-oxide formation', 'Mn-oxide dissolution', 'Mn cycling', 'birnessite formation', 'pyrolusite formation', 'Mn precipitation'], 
'metal_transformation': ['iron reduction', 'manganese oxidation', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'metal deposition', 'metal solubilization', 'mineral dissolution', 'mineral precipitation', 'metal mobilization', 'co-precipitation', 'biosorption'], 
'ph_modulation': ['acid tolerance', 'alkaline tolerance', 'proton pump', 'pH homeostasis', 'pH stress', 'pH regulation', 'acid resistance', 'pH buffering', 'pH gradient', 'acidophile', 'alkaliphile'], 
'exoelectrogenesis': ['exoelectrogen', 'electrochemically active bacteria', 'EAB', 'extracellular respiration', 'electrode respiration', 'bioelectrosynthesis', 'microbial fuel cell'], 
'enzymatic_corrosion': ['enzyme-mediated corrosion', 'metalloenzyme', 'enzyme-catalyzed oxidation', 'enzyme-catalyzed reduction', 'peroxidase', 'laccase', 'oxidoreductase activity'], 
'dealloying': ['selective corrosion', 'dezincification', 'dealuminification', 'parting', 'selective leaching', 'alloy corrosion', 'preferential dissolution'], 
'galvanic_corrosion': ['bimetallic corrosion', 'galvanic couple', 'differential aeration', 'local cell formation', 'microelectrochemical cells', 'electrolytic corrosion'], 
'chloride_attack': ['chloride-induced corrosion', 'pitting initiation', 'localized corrosion', 'crevice corrosion', 'stress corrosion cracking', 'chloride penetration'], 
'microbe_metal_synergy': ['metal-microbe interaction', 'biomineralization', 'biosorption', 'bioaccumulation', 'metal-binding proteins', 'metallothionein', 'metal-enzyme complex'], 
'cathodic_depolarization': ['cathodic depolarization', 'hydrogen depolarization', 'cathodic reaction enhancement', 'electron acceptor consumption', 'cathodic process acceleration'],
'passivity_breakdown': ['passivation breakdown', 'passive film disruption', 'oxide layer dissolution', 'passive layer damage', 'depassivation'],
'concentration_cells': ['concentration cell', 'differential concentration', 'ion concentration gradient', 'oxygen concentration cell', 'metal ion gradient']
}

#===============================================
# 3. Complete pathway_categories 
pathway_categories = {
'oxygen_metabolism': ['o2_consumption', 'aerobic_respiration', 'oxygen reduction', 'oxygen consumption', 'cytochrome oxidase', 'terminal oxidase','oxygen reductase', 'superoxide dismutase', 'catalase','oxidative stress', 'oxygen sensor', 'oxygen tolerance','reactive oxygen species', 'peroxidase', 'oxygen-limited growth' ],
'nitrogen_metabolism': ['denitrification', 'nitrate_reduction', 'nitrite_reduction','nitrogen metabolism', 'nitrate', 'nitrite','ammonia oxidation', 'nitrification', 'nitrogen fixation', 'ammonification', 'nitrous oxide reduction', 'nitric oxide reduction','dissimilatory nitrate reduction', 'nitrite reductase', 'nitrate reductase', 'anammox', 'nitrate respiration', 'nitrite respiration', 'nitric oxide metabolism' ],
'iron_sulfur_redox': ['iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB','Fe-S cluster', 'iron-sulfur cluster', 'ferredoxin','rubredoxin', 'ferritin', 'bacterioferritin', 'PWY-7221', 'HEME-BIOSYNTHESIS-II', 'P125-PWY', 'iron mobilization', 'iron immobilization', 'iron precipitation' ],
'ocre_formation': ['ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation','ferrihydrite', 'goethite', 'magnetite', 'hematite','mineral precipitation', 'iron mineral', 'PWY-7219', 'biogenic iron oxides', 'stalactite formation', 'ochre mats'],
'manganese_processes': ['manganese_reduction', 'mn_redox', 'manganese oxidation','manganese oxide', 'pyrolusite', 'birnessite','manganese cycling', 'manganese mineral', 'manganese transport','Mn-oxide formation', 'Mn-oxide reduction', 'Mn precipitation', 'Mn dissolution' ], 
'sulfur_metabolism': ['sulfur', 'sulfate', 'sulfide', 'thiosulfate', 'sulfite', 'sulfonate','sulfate reduction', 'sulfur oxidation', 'sulfur respiration', 'SRB', 'dsrAB', 'APS reductase', 'sulfide:quinone oxidoreductase','sulfur disproportionation', 'dissimilatory sulfate reduction','PWY-6932','SO4ASSIM-PWY', 'SULFATE-CYS-PWY', 'sulfur globules', 'elemental sulfur', 'polysulfide metabolism', 'sulfur granules' ],
'hydrogen_metabolism': ['h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen production','FeFe-hydrogenase', 'NiFe-hydrogenase', 'hydrogen evolution','hydrogen cycling', 'H2 sensing', 'proton reduction'],
'electron_transfer': ['cytochrome', 'electron transport', 'oxidoreductase', 'redox','electron transfer', 'direct EET', 'indirect EET', 'nanowire','conductive pili', 'electron shuttle', 'mtrABC', 'omcS','c-type cytochrome', 'multi-heme cytochrome', 'flavin','electron bridge', 'solid-state electron transfer', 'interfacial electron transfer'],
'organic_acid_metabolism': ['acetate', 'acetic acid', 'acetyl', 'acetate metabolism', 'acetate production','oxalate', 'oxalic acid', 'oxalate metabolism', 'oxalate production', 'organic acid', 'fatty acid', 'butyric acid', 'butyrate', 'propionate', 'propionic acid','carboxylic acid', 'lactate', 'lactic acid', 'formate', 'formic acid','citrate', 'citric acid', 'succinate', 'succinic acid', 'fumarate', 'fumaric acid','malate', 'malic acid', 'pyruvate', 'pyruvic acid', 'CENTFERM-PWY', 'FERMENTATION-PWY', 'GLYCOLYSIS', 'PWY-5100','GALACTUROCAT-PWY', 'gluconic acid', 'tartaric acid', 'caproic acid', 'valeric acid' ],
'metal_organic_interaction': ['siderophore', 'metal binding', 'metal chelation', 'iron chelation', 'iron complex', 'metal transport', 'metallophore', 'metal organic', 'organometallic', 'iron uptake', 'metal uptake', 'metalloprotein', 'iron-sulfur cluster', 'metal coordination', 'metal sequestration','ferric reductase', 'ferrous oxidase', 'metal homeostasis', 'metal bioavailability', 'metallophore production', 'metal-ligand complex', 'chelator secretion', 'metal bioavailability' ],
'biofilm_formation': ['biofilm', 'exopolysaccharide', 'EPS production', 'extracellular matrix', 'adhesion', 'attachment', 'colonization', 'surface attachment', 'polysaccharide biosynthesis', 'extracellular polymeric substance', 'cell aggregation', 'quorum sensing', 'biofilm maturation','biofilm regulation', 'biofilm dispersion', 'cell-cell adhesion', 'matrix production', 'pellicle', 'floc formation', 'COLANSYN-PWY','EXOPOLYSACC-PWY', 'GLUCOSE1PMETAB-PWY', 'alginate synthesis', 'cellulose biosynthesis', 'amyloid fiber', 'biofilm architecture' ],
'carbon_metabolism': ['carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'pentose phosphate pathway', 'gluconeogenesis', 'carbon assimilation', 'carbon flux','Calvin cycle', 'reductive acetyl-CoA pathway', 'carbon monoxide dehydrogenase','carbon catabolite repression', 'carbon storage', 'glycogen metabolism' ],
'ph_modulation': ['acid', 'alkaline', 'proton pump', 'pH homeostasis', 'pH stress', 'acid tolerance', 'alkaline tolerance', 'proton motive force', 'pH regulation', 'acidic environment','alkaline environment', 'acid resistance', 'proton antiporter','pH gradient maintenance', 'pH sensing', 'pH buffering system' ],
'temp_response': ['heat shock', 'cold shock', 'temperature response', 'thermophilic', 'psychrophilic', 'mesophilic', 'thermal adaptation', 'temperature stress', 'heat stress protein','cold stress protein', 'thermal stability', 'thermotolerance'],
'halogen_related': ['halogen', 'chloride', 'bromide', 'iodide', 'fluoride','halide', 'dehalogenation', 'haloperoxidase', 'haloacid','chlorination', 'bromination', 'halomethane', 'haloalkane','organohalide', 'halotolerance', 'salt tolerance','halophilic', 'chloride transport', 'halide channel','hypochlorite', 'chlorate reduction', 'perchlorate reduction' ],
'methanogenesis': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 'methyl-coenzyme M reductase', 'methanogenic', 'coenzyme F420','methyl-H4MPT', 'CO2 reduction', 'acetoclastic methanogenesis', 'hydrogenotrophic methanogenesis', 'methylotrophic methanogenesis', 'methanotrophy' ],
'silicon_metabolism': [ 'silicon uptake', 'silicate transport', 'silicon cycling', 'diatom frustule','silicification', 'biogenic silica', 'silicon precipitation'],
'phosphorus_metabolism': ['phosphate uptake', 'phosphate transport', 'polyphosphate', 'phosphorus cycling','phosphate precipitation', 'struvite formation', 'phosphate solubilization']
} 

#===============================================
# 4. Complete organic_categories 
organic_categories = {
'degradation': ['degradation', 'breakdown', 'catabolism', 'hydrolysis', 'digestion','decomposition', 'mineralization', 'dissolution', 'depolymerization', 'lysis'], 
'synthesis': ['biosynthesis', 'anabolism', 'synthesis', 'polymerization', 'assembly', 'synthase', 'ligase', 'transferase', 'acetyl-coa', 'malonyltransferase', 'acyltransferase','construction', 'formation', 'production', 'biogenesis', 'neosynthesis'], 
'transport': ['transport', 'uptake', 'export', 'secretion', 'efflux', 'influx', 'permease', 'transporter', 'carrier protein','extrusion', 'import', 'transmembrane transport', 'translocation', 'diffusion', 'facilitated transport'], 
'modification': ['modification', 'conversion', 'transformation', 'transmutation', 'decarboxylase', 'transaminase', 'methyltransferase', 'isomerase', 'lyase','alteration', 'rearrangement', 'isomerization', 'conjugation', 'methylation', 'hydroxylation'], 
'respiration': ['respiration', 'electron transport chain', 'oxidative phosphorylation','terminal electron acceptor', 'cytochrome', 'anaerobic respiration', 'substrate-level phosphorylation'], 
'fermentation': ['fermentation', 'anaerobic metabolism', 'mixed acid fermentation','alcohol fermentation', 'lactic acid fermentation', 'homolactic fermentation', 'heterolactic fermentation'], 
'oxidation': ['oxidation', 'oxidase', 'oxidoreductase', 'dehydrogenase','hydroxylation', 'oxygenase', 'monooxygenase', 'dioxygenase'], 
'reduction': ['reduction', 'reductase', 'hydrogenase', 'dehydrogenase','3-oxoacyl reductase', 'enoyl reductase', 'beta-ketoacyl reductase', 'ketoacyl reductase','nitrate reduction', 'sulfate reduction', 'electron donation', 'reduction potential'],
'fatty_acid_metabolism': ['3-oxoacyl', 'enoyl', 'beta-ketoacyl', 'acyl-acp', 'fatty acid synthesis', 'malonyl'] 
}

#===============================================
# 5. Corrosion-specific metal combinations 
corrosion_synergies = {
'Fe-S': ['iron sulfur', 'Fe-S', 'iron sulfide', 'FeS', 'Fe-S cluster', 'pyrite', 'pyrrhotite'],
'Fe-Cl': ['iron chloride', 'FeCl', 'iron halide', 'ferric chloride', 'ferrous chloride', 'chloride attack', 'pitting corrosion'],
'Fe-C': ['iron carbon', 'FeC', 'iron carbonate', 'siderite', 'carbonate corrosion'], 
'Cu-Fe': ['copper iron', 'Cu-Fe', 'bimetallic', 'galvanic couple', 'bronze', 'brass'], 
'Mn-Fe': ['manganese iron', 'Mn-Fe', 'iron manganese oxide', 'steel manganese'], 
'Ni-Fe': ['nickel iron', 'Ni-Fe', 'stainless steel', 'alloy corrosion'], 
'Cr-Fe': ['chromium iron', 'Cr-Fe', 'stainless steel', 'chromium passivation'], 
'Al-Cu': ['aluminum copper', 'Al-Cu', 'aluminum brass', 'galvanic corrosion'], 
'Zn-Fe': ['zinc iron', 'Zn-Fe', 'galvanized steel', 'sacrificial anode'],
'Fe-CO3': ['iron carbonate', 'siderite', 'bicarbonate corrosion', 'carbonate scaling'],
'Fe-SO4': ['iron sulfate', 'sulfate corrosion', 'gypsum formation'],
'Fe-Ox': ['iron oxalate', 'oxalate corrosion', 'organic acid attack', 'oxidation corrosion'],
'Fe-Ac': ['iron acetate', 'acetate corrosion', 'organic acid attack', 'oxidation corrosion']
}

#===============================================
# 6. COMPLETE functional_categories with HVAC-corrected scoring 
functional_categories = {
# PRIMARY HVAC CORROSION MECHANISMS (1.2-1.5) - Literature-justified scores for HVAC systems

'iron_sulfur_redox': {'terms': ['iron_metabolism', 'sulfur_metabolism', 'iron_oxidation', 'iron_reduction', 'iron reduction', 'ferric reduction', 'sulfate reduction', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB', 'iron-sulfur cluster', 'rubredoxin', 'ferredoxin'], 
'score': 1.5, # HIGHEST for HVAC - SRB major problem in closed loops
'justification': 'Beech, I.B., Gaylarde, C.C. (1999). Recent advances in the study of biocorrosion. Rev Microbiol, 30(3):177-190. SRB critical in closed water systems'
},
'acid_production': {'terms': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 'lactate metabolism', 'formate production', 'carbonic acid', 'citric acid'], 
'score': 1.4, # VERY HIGH for HVAC - organic acids major issue
'justification': 'Videla, H.A., Herrera, L.K. (2005). Microbiologically influenced corrosion: looking to the future. Int Microbiol, 8(3):169-180'
},
'o2_consumption': {
'terms': ['o2_consumption', 'aerobic_respiration', 'oxygen reduction', 'oxygen consumption', 'oxygen depletion', 'oxygen-limited growth'], 
'score': 1.2, # VERY HIGH for HVAC - oxygen depletion creates aggressive conditions
'justification': 'Geesey, G.G., Bremer, P.J. (1990). Biofouling of engineered water systems. Biotechnol Bioeng, 36(10):1039-1046'
},
'electron transfer & redox': {'terms': ['direct_eet', 'redox', 'electron_transfer', 'omc', 'deet', 'cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase', 'extracellular electron transfer'], 
'score': 1.3, # Important but secondary to chemical mechanisms in HVAC
'justification': 'Jones, D.A. (1996). Principles and Prevention of Corrosion, 2nd ed. All corrosion fundamentally involves electron transfer'},

'biofilm_formation': {'terms': ['biofilm_formation', 'metal_chelation', 'quorum_sensing', 'extracellular_matrix', 'EPS', 'surface_disruption', 'polysaccharide', 'adhesin', 'biofilm', 'extracellular polymeric substance', 'curli', 'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 'attachment', 'alginate', 'cellulose'], 
'score': 1.3, #  biofilms major problem in closed systems
'justification': 'Borenstein, S.W. (1994). Microbiologically Influenced Corrosion Handbook. Industrial Press. Biofilms critical in closed water systems'
},

'metal binding / chelation': {'terms': ['metal_chelation', 'metal_binding', 'siderophore', 'complexation', 'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelate formation', 'metal ligand'], 
'score': 1.2, # Important in closed loops
'justification': 'Herrera, L.K., Videla, H.A. (2009). Role of iron-reducing bacteria in corrosion and protection of carbon steel. Int Biodeterior Biodegradation, 63(7):891-895'
},
# ochre score based on organic matter involvement in possible formation of Fe complexes similar to green rust 
'ocre': {'terms': ['ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'ferrihydrite', 'goethite', 'magnetite'], 
'score': 1.1, # represents active Fe-organic complexation 
'justification': 'Cornell, R.M., Schwertmann, U. (2003). The Iron Oxides. Green rust and iron-organic complexes indicate active corrosion processes'
},
# SECONDARY MECHANISMS (0.7-1.1) - HVAC-adjusted scores
'nitrogen_reduction': {'terms': ['denitrification', 'nitrate_reduction', 'nitrite_reduction', 'nitrate respiration', 'nitrous oxide reduction', 'nitric oxide metabolism'], 
'score': 1.0, # Moderate for HVAC - affects redox conditions
'justification': 'Flemming, H.C. (1996). Economically relevant microorganisms in technical systems. Materials and Corrosion, 47(7):391-398'
},
'manganese_reduction': {'terms': ['manganese_reduction', 'mn_redox', 'manganese oxidation', 'birnessite formation', 'pyrolusite formation'], 
'score': 1.0, # Moderate for HVAC - Mn cycling affects redox
'justification': 'Tebo, B.M., et al. (2004). Biogenic manganese oxides: properties and mechanisms of formation. Annu Rev Earth Planet Sci, 32:287-328'
},
'sulfide_production': {'terms': ['sulfide', 'sulfur_reduction', 'desulfovibrio', 'hydrogen sulfide', 'H2S', 'sulfide precipitation'], 
'score': 1.0, # Direct corrosive attack but secondary to SRB iron-sulfur redox
'justification': 'Sulfide causes direct corrosive attack on iron and steel surfaces'
},
'methanogenesis': {'terms': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 'methyl-coenzyme M reductase'], 
'score': 0.7, # REDUCED for HVAC - uncommon in aerobic systems
'justification': 'Mori, K., et al. (2010). Methanogens in microbiologically influenced corrosion: a review. Microorganisms, 8(7):995'
},
'h2_consumption': {'terms': ['h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism'], 
'score': 0.6, # Enables cathodic depolarization
'justification': 'Javaherdashti, R. (2008). Microbiologically Influenced Corrosion: An Engineering Insight. Springer'
},
'fumarate_formation': {'terms': ['fumarate', 'propionibacterium', 'fumarate reduction', 'fumarate respiration'], 
'score': 0.5 # Lower priority in HVAC
},
'nitrogen_metabolism': {'terms': ['nitrogen_metabolism', 'nitrification', 'ammonification', 'nitrogen fixation', 'anammox'], 
'score': 1.0,
'justification': 'Nitrogen cycling influences local redox and pH conditions in water systems'
},
'exoelectrogenesis': {'terms': ['exoelectrogen', 'electrochemically active bacteria', 'EAB', 'extracellular respiration', 'electrode respiration'], 
'score': 0.9, # less critical than industrial systems
'justification': 'Logan, B.E. (2009). Exoelectrogenic bacteria that power microbial fuel cells. Nat Rev Microbiol, 7(5):375-381'
},
'enzymatic_metal_oxid': {'terms': ['metalloenzyme', 'enzyme-catalyzed oxidation', 'peroxidase', 'laccase', 'oxidoreductase activity', 'enzyme-mediated corrosion'], 
'score': 0.8, # secondary to chemical mechanisms
'justification': 'Enzymes can directly catalyze metal oxidation reactions'
},
'metal_precipitation': {'terms': ['metal precipitation', 'biomineralization', 'mineral formation', 'crystal nucleation', 'metal immobilization'], 
'score': 1.0,
'justification': 'Biomineralization affects local chemistry and surface properties'
},
'ph_buffering': {'terms': ['pH buffering', 'pH homeostasis', 'acid resistance', 'alkaline tolerance', 'pH gradient'], 
'score': 0.8,
'justification': 'pH regulation affects local corrosion conditions'
},
'chloride_interactions': {'terms': ['chloride attack', 'chloride-induced corrosion', 'pitting initiation', 'chloride penetration', 'halide corrosion'], 
'score': 0.7, # REDUCED for HVAC - less critical in closed systems vs marine
'justification': 'Chloride attack less common in closed-loop HVAC systems with treated water'
},
'dealloying_mechanisms': {'terms': ['selective corrosion', 'dezincification', 'dealuminification', 'preferential dissolution', 'parting'], 
'score': 0.6, # REDUCED for HVAC - brass/bronze less common in modern systems
'justification': 'Newman, R.C. (2019). The dissolution and passivation kinetics of stainless steels containing molybdenum. Corros Sci, 25(5):331-339'
}
}

#===============================================
# 7. COMPLETE corrosion_keyword_groups 
corrosion_keyword_groups = {
'iron_sulfur_redox': ['iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 'ferric reduction', 'sulfate reduction', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB', 'pyrite', 'pyrrhotite'], 
'ocre': ['ocre', 'iron oxide', 'iron deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'ferrihydrite', 'goethite'], 
'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 'lactate metabolism', 'formate production', 'carbonic acid', 'citric acid'], 
'electron_transfer': ['direct eet', 'redox', 'electron transfer', 'omc', 'deet', 'cytochrome', 'electron transfer', 'conductive pili', 'nanowire', 'mtrABC', 'omcS', 'oxidoreductase', 'redox', 'reductase', 'oxidase', 'extracellular electron transfer'], 
'biofilm': ['biofilm formation', 'quorum sensing', 'extracellular matrix', 'EPS', 'surface disruption', 'polysaccharide', 'adhesin', 'biofilm', 'extracellular polymeric substance', 'curli', 'exopolymer', 'attachment', 'colonization', 'alginate', 'cellulose', 'biofilm architecture'], 
'sulfide': ['sulfide', 'sulfur reduction', 'desulfovibrio', 'h2s', 'hydrogen sulfide', 'sulfide precipitation'], 
'metal_binding': ['metal chelation', 'metal binding', 'siderophore', 'complexation', 'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelate formation', 'biogenic_iron_oxide', 'microbial_iron_precipitation'],
'nitrogen': ['denitrification', 'nitrate reduction', 'nitrite reduction', 'nitrogen metabolism', 'nitrate', 'nitrite', 'nitrification', 'ammonification'], 
'manganese': ['manganese reduction', 'mn redox', 'manganese oxidation', 'mn reduction', 'birnessite', 'pyrolusite', 'manganese mineral'], 
'methanogenesis': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 'hydrogenotrophic methanogenesis', 'acetoclastic methanogenesis'], 
'fumarate': ['fumarate', 'propionibacterium', 'fumarate reduction', 'fumarate respiration'],
'hydrogen': ['h2 consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'h2 oxidation', 'h2ase', 'hydrogen evolution', 'hydrogen metabolism'], 
'oxygen': ['o2 consumption', 'aerobic respiration', 'oxygen reduction', 'oxygen consumption', 'oxygen depletion', 'oxygen-limited', 'oxygen stress'], 
'corrosion_general': ['corrosion', 'deterioration', 'pitting', 'microbially influenced corrosion','MIC', 'metal deterioration', 'galvanic corrosion', 'crevice corrosion', 'stress corrosion'], 
'exoelectrogen': ['exoelectrogen', 'electrochemically active', 'EAB', 'extracellular respiration', 'electrode respiration', 'bioelectrosynthesis'], 
'enzymatic_corrosion': ['enzyme-mediated corrosion', 'metalloenzyme', 'laccase', 'peroxidase', 'enzyme-catalyzed', 'oxidoreductase'], 
'metal_precipitation': ['metal precipitation', 'biomineralization', 'mineral formation', 'mineral precipitation', 'crystal nucleation'], 
'chloride_attack': ['chloride attack', 'chloride corrosion', 'pitting initiation', 'chloride penetration', 'hypochlorite', 'halide corrosion'], 
'ph_regulation': ['pH buffering', 'pH homeostasis', 'acid resistance', 'pH gradient', 'proton pump', 'pH regulation'],
'dealloying': ['dealloying', 'selective corrosion', 'dezincification', 'preferential dissolution', 'parting', 'alloy corrosion'] 
}

#===============================================
# 8. Standard mapping 
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

# References justification:
scoring_justifications = {
'iron_sulfur_redox_0.85': [
'Beech, I.B., Gaylarde, C.C. (1999). Recent advances in the study of biocorrosion. Rev Microbiol, 30(3):177-190',
'SRB proliferate in oxygen-depleted zones of closed water systems'
],
'acid_production_0.8': [
'Videla, H.A., Herrera, L.K. (2005). Microbiologically influenced corrosion: looking to the future. Int Microbiol, 8(3):169-180',
'Organic acid production is dominant MIC mechanism in building water systems'
],
'biofilm_0.75': [
'Borenstein, S.W. (1994). Microbiologically Influenced Corrosion Handbook. Industrial Press',
'Biofilms are the primary cause of MIC in closed water systems'
],
'oxygen_0.65': [
'Geesey, G.G., Bremer, P.J. (1990). Biofouling of engineered water systems. Biotechnol Bioeng, 36(10):1039-1046',
'Oxygen depletion under biofilms creates aggressive corrosion conditions'
],
'ochre_0.3_failure_analysis': ['Cornell, R.M., Schwertmann, U. (2003). The Iron Oxides. Wiley-VCH',
'failure_analysis: possible Fe-organic anions (oxalate/acetate) complexes represents active corrosion processes'
]
}
