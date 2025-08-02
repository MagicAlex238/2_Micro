# 1. Complete metal_terms 
metal_terms = {	
'iron': ['Fe2+', 'Fe3+', 'iron', 'ferrous', 'ferric', 'heme', 'iron-sulfur', 'rust', 'ochre', 'iron oxide', 'siderophore', 'ferritin'],	
'manganese': ['Mn2+',  'manganese', 'mn', 'manganous', 'manganic', 'manganese oxidation', 'manganese oxide', 'MnO2'],	
'copper': ['Cu+', 'Cu2+', 'copper', 'cupric', 'cuprous', 'copper oxide', 'copper corrosion'],	
'nickel': ['Ni2+', 'nickel', 'nickelous', 'nickel oxidation', 'nickel reduction'],	
'cobalt': ['Co2+',  'cobalt', 'cobaltous', 'cobalamin', 'vitamin B12'],	
'magnesium': ['Mg2+', 'magnesium', 'magnesium oxide'],	
'calcium': ['Ca2+', 'calcium', 'calcium carbonate', 'calcite', 'calcium precipitation'],	
'Mo': ['Mo',  'molybdenum', 'molybdopterin', 'molybdenum cofactor'],	
'V5+': ['V5+', 'vanadium', 'vanadate', 'vanadyl'],	
'Al3+': ['Al3+', 'aluminum', 'aluminate', 'aluminum oxide'],	
'Cr3+': ['Cr3+', 'Cr6+', 'chromium', 'chromate', 'dichromate', 'chromium oxide'],	
'zinc': ['Zn2+', 'zinc', 'zinc finger', 'zinc oxide'],	
'sodium': ['Na+', 'sodium', 'NaCl', 'sodium transport', 'sodium gradient'],	
'potassium': ['K+', 'potassium', 'KCl', 'potassium transport', 'potassium channel'],	
'selenium': ['selenium', 'Se', 'selenocysteine', 'selenoprotein', 'selenite'],	
'barium': ['Ba2+', 'barium', 'barium sulfate', 'barite'],	
'strontium': ['Sr2+', 'strontium', 'strontium carbonate', 'strontium sulfate'], 	
'lead': ['Pb2+', 'Pb4+', 'lead', 'plumbous', 'plumbic', 'lead oxide'], 	
'arsenic': ['As3+', 'As5+', 'arsenic', 'arsenite', 'arsenate', 'arsenic oxidation'], 	
'mercury': ['Hg2+', 'Hg+', 'mercury', 'mercuric', 'mercurous', 'mercury methylation'], 	
}
# Corrosion-specific metal combinations 
corrosion_synergies= {
'Fe-S': ['iron_sulfur', 'Fe-S','iron sulfide','FeS', 'Fe-S cluster'],	
'Fe-Cl': ['iron chloride', 'FeCl', 'iron halide', 'ferric chloride'],	
'Fe-C': ['iron carbon', 'FeC', 'iron carbonate', 'siderite'],	
'Cu-Fe': ['copper iron', 'Cu-Fe', 'bimetallic', 'galvanic couple'],	
'Mn-Fe': ['manganese iron', 'Mn-Fe', 'iron manganese oxide'],
'Ni-Fe': ['Ni-Fe'],  	
'Cr-Fe': ['chromium iron', 'Cr-Fe', 'stainless steel', 'chromium passivation'], 	
'Al-Cu': ['aluminum copper', 'Al-Cu', 'aluminum brass', 'galvanic corrosion'], 	
'Zn-Fe': ['zinc iron', 'Zn-Fe', 'galvanized steel', 'sacrificial anode'],	
'Fe-CO3': ['iron carbonate', 'siderite', 'bicarbonate corrosion', 'carbonate scaling'],	
'Fe-SO4': ['iron sulfate', 'sulfate corrosion', 'gypsum formation'],	
'Fe-Ox': ['iron oxalate', 'oxalate corrosion', 'organic acid', 'oxidation corrosion'],	
'Fe-Ac': ['iron acetate', 'acetate corrosion', 'organic acid', 'oxidation corrosion']}	

# 3. COMPLETE functional_categories with HVAC-corrected scoring and justification
functional_categories =	{
'o2_consumption': {'terms': ['o2_consumption', 'aerobic_respiration', 'oxygen reduction', 'oxygen consumption', 'cytochrome oxidase', 'oxidase', 'terminal oxidase',  'oxygen reductase',  'superoxide dismutase',  'catalase',  'oxidative stress', 'oxygen sensor',  'oxygen tolerance', 'oxygen consum', 'oxygen scavenging', 'oxygen stress', 'peroxidase', 'oxidative phosphorylation', 'NADH dehydrogenase', 'CYTOCHROME-C-OXIDASE','ATPSYN-RXN', "TCA cycle IV (2-oxoglutarate decarboxylase)","TCA cycle V (2-oxoglutarate:ferredoxin oxidoreductase)",  "TCA cycle VI (obligate autotrophs)", "TCA cycle VIII (helicobacter)", "superpathway of glyoxylate bypass and TCA"], 'score': 1.2, 'justification': 'Geesey, G.G., Bremer, P.J. (1990). Biofouling of engineered water systems. Biotechnol Bioeng, 36(10):1039-1046'},# VERY HIGH for HVAC - oxygen depletion creates aggressive conditions
'nitrogen_metabolism':  {'terms': ['nitrate_reduction', 'nitrite_reduction', 'denitrification', 'nitrification', 'nitrate respiration', 'nitrite respiration', 'nitrous oxide reduction', 'ammonia oxidation', 'anammox', 'nitrogen fixation', 'ammonification', 'nitrogen metabolism',  'nitrate',  'nitrite','dissimilatory nitrate reduction',  'nitrite reductase', 'nitrate reductase', "L-arginine biosynthesis II (acetyl cycle)", "L-ornithine biosynthesis", "L-histidine biosynthesis", "L-methionine biosynthesis III", "L-isoleucine biosynthesis I (from threonine)", "L-isoleucine biosynthesis II", "L-isoleucine biosynthesis III", "L-isoleucine biosynthesis IV", "L-lysine biosynthesis I", "L-lysine biosynthesis III", "L-lysine biosynthesis VI", "L-valine biosynthesis", "L-tryptophan biosynthesis", "superpathway of L-isoleucine biosynthesis I", "superpathway of L-phenylalanine biosynthesis", "superpathway of L-tyrosine biosynthesis", "superpathway of L-threonine biosynthesis" ], 'score': 1.0,'justification': 'Flemming, H.C. (1996). Economically relevant microorganisms in technical systems. Materials and Corrosion, 47(7):391-398 and  Washizu, N., et al. (2004). MIC by nitrate-reducing bacteria in water injection systems. Corrosion, 60(4):336-342'}, # Moderate for HVAC - affects redox conditions	
'iron_metabolism': {'terms': ['corrosion', 'MIC', 'ferric reduc', 'SRB', 'ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'iron oxid', 'ferrous oxid', 'ferric', 'iron uptake', 'iron transport', 'iron storage',  'iron homeostasis', 'siderophore production', 'iron_sulfur_redox', 'ferredoxin',  'rubredoxin', 'ferritin', 'bacterioferritin',  'PWY-7221',  'PWY-7219', 'HEME-BIOSYNTHESIS-II',  'P125-PWY', 'iron mobilization', 'iron immobilization', 'ferrihydrite', 'goethite', 'magnetite', 'hematite', 'iron mineral', 'biogenic iron oxides', 'stalactite formation', 'ochre mats', 'superpathway of tetrahydrofolate biosynthesis and salvage',  "heme biosynthesis II (anaerobic)", "tetrapyrrole biosynthesis I (from glutamate)", "tetrapyrrole biosynthesis II (from glycine)","flavin biosynthesis I (bacteria and plants)","chorismate biosynthesis I", "chorismate biosynthesis from 3-dehydroquinate", 'enzyme-mediated corrosion'], 'score': 1.5, 'justification': 'Beech, I.B., Gaylarde, C.C. (1999). Recent advances in the study of biocorrosion. Rev Microbiol, 30(3):177-190. SRB critical in closed water systems AND Cornell, R.M., Schwertmann, U. (2003). The Iron Oxides. Green rust and iron-organic complexes indicate active corrosion processes'},# represents active Fe-organic complexation '},	# HIGHEST for HVAC - SRB major problem in closed loops
'sulfur_metabolism': {'terms': ['sulphur_metabolism', 'sulfur_metabolism', 'sulfate reduc', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB', 'dsrAB', 'APS reductase', 'sulfide','quinone oxidoreductase', 'dissimilatory sulfate reduction', 'sulfur globules', 'elemental sulfur', 'polysulfide metabolism', 'sulfur granules', 'PWY-6932', 'SO4ASSIM-PWY', 'SULFATE-CYS-PWY', 'sulfide_production', 'sulfonate',  'sulfur_reduction', 'desulfovibrio', 'sulfur disproportionation', 'sulfate-reducing bacteria', 'sulfur respiration', "sulfate reduction I (assimilatory)", "superpathway of sulfate assimilation and cysteine biosynthesis"], 'score': 1.5, 'justification': 'Sulfide causes direct corrosive attack on iron and steel surfaces'},# Direct corrosive attack but secondary to SRB iron-sulfur redox
'h2_consumption': {'terms': ['h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism', 'hydrogen production',  'FeFe-hydrogenase',  'NiFe-hydrogenase',  'hydrogen evolution',  'hydrogen cycling',  'H2 sensing',  'proton reduction', "methylerythritol phosphate pathway I", "methylerythritol phosphate pathway II"], 'score': 0.7, 'justification': 'Javaherdashti, R. (2008). Microbiologically Influenced Corrosion: An Engineering Insight. Springer and Enning, D., Garrelfs, J. (2014). Corrosion of iron by sulfate-reducing bacteria: new views of an old problem. Appl Environ Microbiol, 80(4):1226-1236'},# H2 consumption by SRB accelerates corrosion through cathodic depolarization	
'carbon_metabolism': {'terms': ['carbon_metabolism', 'carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'carbon flux', 'carbon assimilation', 'pentose phosphate pathway', 'gluconeogenesis', 'Calvin cycle', 'reductive acetyl-CoA pathway', 'carbon monoxide dehydrogenase', 'hydrocarbon degradation', 'aromatic degradation', 'alcohol metabolism', 'organic matter degradation', 'VFA production', 'propionate', 'butyrate', 'valerate', 'caproate', 'ACETYL-COA-ACETYLTRANSFER-RXN','METHYLACETOACETYLCOYTHIOL-RXN','ACETOLACTSYN-RXN','ACETOOHBUTSYN-RXN','ACETYL-COA-CARBOXYLTRANSFER-RXN','ACYLCOASYN-RXN','N10-formyl-tetrahydrofolate biosynthesis', 'glycolysis III','Calvin-Benson-Bassham cycle', 'gluconeogenesis', "glycolysis I (from glucose 6-phosphate)", "glycolysis II (from fructose 6-phosphate)","glycolysis III (from glucose)","gluconeogenesis I","pentose phosphate pathway (non-oxidative branch)","Calvin-Benson-Bassham cycle","pyruvate fermentation to isobutanol (engineered)","superpathway of branched amino acid biosynthesis","superpathway of aromatic amino acid biosynthesis","superpathway of adenosine nucleotides de novo biosynthesis I","superpathway of adenosine nucleotides de novo biosynthesis II","superpathway of guanosine nucleotides de novo biosynthesis I","superpathway of guanosine nucleotides de novo biosynthesis II"  ], 'score': 0.5,'justification': 'Pope, D.H. (1986). A study of microbiologically influenced corrosion in nuclear power plants. Electric Power Research Institute'},
'organic_acid_metabolism': {'terms':  ['acetate', 'acetic acid', 'acetyl', 'acetate metabolism', 'acetate production', 'oxalate', 'oxalic acid', 'oxalate metabolism', 'oxalate production', 'organic acid', 'fatty acid', 'butyric acid', 'butyrate', 'propionate', 'propionic acid', 'carboxylic acid', 'lactate', 'lactic acid', 'formate', 'formic acid', 'citrate', 'citric acid', 'succinate', 'succinic acid', 'fumarate', 'fumaric acid', 'malate', 'malic acid', 'pyruvate', 'pyruvic acid', 'acidification', 'fermentation', 'CENTFERM-PWY', 'FERMENTATION-PWY', 'GLYCOLYSIS', 'PWY-5100', 'GALACTUROCAT-PWY', 'fatty acid β-oxidation I', 'fatty acid elongation',"fatty acid salvage","stearate biosynthesis II (bacteria and plants)","palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)","cis-vaccenate biosynthesis","oleate biosynthesis IV (anaerobic)","gondoate biosynthesis (anaerobic)","mycolate biosynthesis"], 'score': 1.4, 'justification': 'Videla, H.A., Herrera, L.K. (2005). Microbiologically influenced corrosion: looking to the future. Int Microbiol, 8(3):169-180'}, # VERY HIGH for HVAC - organic acids major issue	
'metal_binding_chelation': {'terms': ['metal_chelation', 'metal_binding', 'siderophore', 'complexation', 'iron chelation', 'enzymatic_metal_oxid',   'chelator', 'metallophore', 'iron complex', 'metal transport', 'metal oxide', 'iron oxide deposits', 'metal deposition', 'metal solubilization', 'mineral dissolution', 'mineral precipitation', 'chelation', 'metal complexation', 'metal sequestration' , 'metal_organic_interaction', 'metal organic', 'metal homeostasis', 'organometallic',  'iron uptake', 'metal uptake', 'metalloprotein',  'iron-sulfur cluster', 'metal coordination', 'ferric reductase', 'ferrous oxidase', 'metal homeostasis', 'mineral dissolution', 'laccase', 'mineral precipitation', 'copper reduction', 'nickel oxidation', 'chromium reduction', 'crystal nucleation', 'metal immobilization', "coenzyme A biosynthesis I","pantothenate and coenzyme A biosynthesis I","phosphopantothenate biosynthesis I","NAD biosynthesis I (from aspartate)","thiamin salvage II", 'metalloenzyme', 'enzyme-catalyzed oxidation'], 'score': 1.2, 'justification': 'Herrera, L.K., Videla, H.A. (2009). Role of iron-reducing bacteria in corrosion and protection of carbon steel. Int Biodeterior Biodegradation, 63(7):891-895'},# Important in closed loops
'biofilm_formation': {'terms': ['biofilm_formation', 'metal_chelation', 'quorum_sensing', 'extracellular_matrix', 'exopolysaccharide', 'EPS production', 'EPS', 'surface_disruption', 'polysaccharide', 'adhesin', 'biofilm', 'EPS', 'extracellular polymeric substance', 'curli', 'exopolymer','extracellular matrix', 'adhesion', 'colonization', 'attachment', 'surface', 'adherence', 'biofilm maturation', 'biofilm regulation', 'biofilm dispersion', 'cell-cell adhesion', 'surface attachment', 'polysaccharide biosynthesis', 'cell aggregation', 'matrix production', 'pellicle', 'floc formation', 'COLANSYN-PWY', 'EXOPOLYSACC-PWY', 'GLUCOSE1PMETAB-PWY', 'alginate', 'cellulose', 'lipid metabolism', 'fatty acid synthesis', 'fatty acid degradation', 'biosurfactant', 'VFA', 'volatile fatty acid', 'propionate', 'butyrate', 'oleaginous', 'lipid accumulation', '3-oxoacyl', '3-oxoacyl-(acyl-carrier-protein)','3-oxocerotoyl-[acp] reductase','3-oxo-cis-Δ7-tetradecenoyl-[acp] reductase','3-oxo-cis-Δ9-hexadecenoyl-[acp] reductase','3-oxo-glutaryl-[acp] methyl ester reductase','3-oxo-pimeloyl-[acp] methyl ester reductase','3-oxo-docosapentaenoyl [acp][c]','(5Z)-3-oxo-tetradec-5-enoyl-[acyl-carrier-protein] reductase','(7Z)-3-oxo-hexadec-7-enoyl-[acp] reductase','(9Z)-3-oxo-octadec-9-enoyl-[acp] reductase','(11Z)-3-oxo-icos-11-enoyl-[acp] reductase','acetoacetyl-[acyl-carrier protein] reductase','3-hydroxyhexanoyl-[acyl-carrier protein] reductase','3-oxo-octanoyl-[acyl-carrier protein] reductase','3-oxo-decanoyl-[acyl-carrier protein] reductase','LINOLENOYL-RXN', 'coenzyme A biosynthesis I', "peptidoglycan biosynthesis I (meso-diaminopimelate containing)","peptidoglycan biosynthesis III (mycobacteria)","UDP-N-acetylmuramoyl-pentapeptide biosynthesis I (meso-diaminopimelate containing)","UDP-N-acetylmuramoyl-pentapeptide biosynthesis II (lysine-containing)","dTDP-L-rhamnose biosynthesis I","O-antigen building blocks biosynthesis (E. coli)","phosphatidylglycerol biosynthesis I (plastidic)","phosphatidylglycerol biosynthesis II (non-plastidic)","CDP-diacylglycerol biosynthesis I","CDP-diacylglycerol biosynthesis II"], 'score': 1.3, 'justification': 'Borenstein, S.W. (1994). Microbiologically Influenced Corrosion Handbook. Industrial Press. Biofilms critical in closed water systems'},#  biofilms major problem in closed systems	
'manganese_processes': {'terms': ['manganese_reduction', 'mn_redox', 'manganese oxidation', 'manganese oxide',  'pyrolusite',  'birnessite',  'manganese cycling',  'manganese mineral',  'manganese transport', 'Mn-oxide formation', 'Mn-oxide reduction', 'Mn precipitation', 'Mn dissolution'],'score': 1.0,'justification': 'Tebo, B.M., et al. (2004). Biogenic manganese oxides: properties and mechanisms of formation. Annu Rev Earth Planet Sci, 32:287-328'}, 
'methanogenesis': {'terms': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 'methyl-coenzyme M reductase', 'methanogenic', 'coenzyme F420',  'methyl-H4MPT', 'CO2 reduction', 'acetoclastic methanogenesis'],'score': 0.6,'justification': 'Mori, K., et al. (2010). Methanogens in microbiologically influenced corrosion: a review. Microorganisms, 8(7):995'},# REDUCED for HVAC - uncommon in aerobic systems	
'fumarate_formation': {'terms': ['fumarate', 'propionibacterium'],  'score': 0.5},# Lower priority in HVAC. 
'phosphorus_metabolism': {'terms': ['phosphate transport',  'polyphosphate', 'phosphite oxidation', 'organophosphonate metabolism',  "UMP biosynthesis", "pyrimidine deoxyribonucleotides de novo biosynthesis I","pyrimidine deoxyribonucleotide phosphorylation","superpathway of pyrimidine nucleobases salvage","superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis","adenosine ribonucleotides de novo biosynthesis","adenosine deoxyribonucleotides de novo biosynthesis II","guanosine ribonucleotides de novo biosynthesis","guanosine deoxyribonucleotides de novo biosynthesis II","superpathway of pyrimidine ribonucleotides de novo biosynthesis","superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis (E. coli)"], 'score': 0.5,'justification': 'Beech, I.B., Sunner, J. (2004). Biocorrosion: towards understanding interactions between biofilms and metals. Curr Opin Biotechnol, 15(3):181-186'},	
}

# 4.The following is a mapping no a dictionary 
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
'phosphate': 'po4-3',
'nitrate': 'NO3-',
'nitrite': 'NO2',
'chloride': 'Cl-'
}
# 5. Complete corrosion_mechanisms, no to be use on scoring system
corrosion_mechanisms = { 
'o2_consumption': ['oxidase', 'oxygen reduction', 'aerobic respiration', 'oxygen consumtion', 'oxygen scavenging', 'oxygen stress', 'oxidative phosphorylation', 'respiratory burst', 'oxygen depletion'], 
'nitrogen_metabolism': ['denitrification', 'nitrification', 'nitrate reduction', 'nitrite reduction', 'nitrate respiration', 'nitrite respiration', 'nitrous oxide reduction', 'ammonia oxidation', 'anammox', 'nitrogen fixation', 'ammonification'], 
'iron_metabolism': ['iron reducction', 'ferric reducction', 'iron oxidation', 'ferrous oxidation', 'iron uptake', 'iron transport', 'iron storage', 'iron homeostasis', 'siderophore production', 'iron-sulfur cluster', 'ferritin', 'bacterioferritin', 'heme biosynthesis'], 
'sulfur_metabolism': ['sulfate reduc', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB', 'sulfur disproportionation', 'sulfate-reducing bacteria', 'sulfur respiration', 'sulfite reduction', 'elemental sulfur', 'polysulfide'], 
'h2_consumption': ['hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen metabolism'],
'carbon_metabolism': ['carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'carbon flux', 'carbon assimilation', 'carbon catabolite repression', 'carbonation', 'carbonate precipitation'], 
'acid_production': ['acid', 'acidification', 'fermentation', 'lactic acid', 'formic acid', 'acetic acid', 'oxalic acid', 'organic acid', 'acetate production', 'lactate metabolism', 'formate production', 'proton generation', 'low pH', 'carbonic acid', 'citric acid', 'gluconic acid'], 
'metal_chelation': ['siderophore', 'metal binding', 'chelator', 'metallophore', 'iron complex', 'metal transport', 'chelation', 'metal complexation', 'metal sequestration', 'chelate formation', 'metal ligand'], 
'biofilm_formation': ['polysaccharide', 'adhesin', 'biofilm', 'EPS', 'extracellular polymeric substance', 'curli', 'exopolymer', 'extracellular matrix', 'adhesion', 'colonization', 'attachment', 'surface', 'adherence', 'biofilm maturation', 'quorum sensing', 'alginate production', 'cellulose production',
                          'quorum sensing', 'autoinducer', 'AI-2', 'AHL', 'acyl homoserine lactone','c-di-GMP', 'biofilm matrix', 'eDNA', 'extracellular DNA', 'biofilm dispersal', 'matrix metalloprotease', 'biofilm regulation', 'stringent response', '(p)ppGpp', 'biofilm architecture'], 
'manganese_metabolism': ['manganese oxidation', 'manganese reduction', 'Mn-oxide formation', 'Mn-oxide dissolution', 'Mn cycling', 'birnessite formation', 'pyrolusite formation', 'Mn precipitation'], 
'metal_transformation': ['iron reduction', 'manganese oxidation', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation', 'metal deposition', 'metal solubilization', 'mineral dissolution', 'mineral precipitation', 'metal mobilization', 'co-precipitation', 'biosorption'], 

}

#===============================================
# 6. Complete pathway_categories, , no to be use on scoring system
pathway_categories = {
'oxygen_metabolism': ['o2_consumption', 'aerobic_respiration', 'oxygen reduction', 'oxygen consumption', 'cytochrome oxidase', 'terminal oxidase','oxygen reductase', 'superoxide dismutase', 'catalase','oxidative stress', 'oxygen sensor', 'oxygen tolerance','reactive oxygen species', 'peroxidase', 'oxygen-limited growth', 'solid-state respiration', 'electrode respiration'],
'nitrogen_metabolism': ['denitrification', 'nitrate_reduction', 'nitrite_reduction','nitrogen metabolism', 'nitrate', 'nitrite','ammonia oxidation', 'nitrification', 'nitrogen fixation', 'ammonification', 'nitrous oxide reduction', 'nitric oxide reduction','dissimilatory nitrate reduction', 'nitrite reductase', 'nitrate reductase', 'anammox', 'nitrate respiration', 'nitrite respiration', 'nitric oxide metabolism',   'L-arginine biosynthesis II (acetyl cycle)','L-ornithine biosynthesis',  'L-histidine biosynthesis','L-methionine biosynthesis III','L-isoleucine biosynthesis I (from threonine)', 'L-isoleucine biosynthesis II','L-isoleucine biosynthesis III','L-isoleucine biosynthesis IV','L-lysine biosynthesis I','L-lysine biosynthesis III', 'L-lysine biosynthesis VI','L-valine biosynthesis', 'L-tryptophan biosynthesis', 'superpathway of L-isoleucine biosynthesis I', 'superpathway of L-phenylalanine biosynthesis','superpathway of L-tyrosine biosynthesis',  'superpathway of L-threonine biosynthesis'],
'iron_sulfur_redox': ['iron metabolism', 'sulfur metabolism', 'iron oxidation', 'iron reduction', 'iron reduc', 'ferric reduc', 'sulfate reduc', 'sulfide', 'sulfite', 'thiosulfate', 'sulfur oxidation', 'SRB','Fe-S cluster', 'iron-sulfur cluster', 'ferredoxin','rubredoxin', 'ferritin', 'bacterioferritin', 'PWY-7221', 'HEME-BIOSYNTHESIS-II', 'P125-PWY', 'iron mobilization', 'iron immobilization', 'iron precipitation',   'PWY-7221', 'HEME-BIOSYNTHESIS-II', 'P125-PWY', 'heme biosynthesis II (anaerobic)', 'tetrapyrrole biosynthesis I (from glutamate)', 'tetrapyrrole biosynthesis II (from glycine)' ],
'ocre_formation': ['ocre', 'iron_oxide', 'iron_deposit', 'metal oxide', 'ochre formation', 'iron oxide deposits', 'iron precipitation', 'rust formation','ferrihydrite', 'goethite', 'magnetite', 'hematite','mineral precipitation', 'iron mineral', 'PWY-7219', 'biogenic iron oxides', 'stalactite formation', 'ochre mats'],
'manganese_processes': ['manganese_reduction', 'mn_redox', 'manganese oxidation','manganese oxide', 'pyrolusite', 'birnessite','manganese cycling', 'manganese mineral', 'manganese transport','Mn-oxide formation', 'Mn-oxide reduction', 'Mn precipitation', 'Mn dissolution' ], 
'sulfur_metabolism': ['sulfur', 'sulfate', 'sulfide', 'thiosulfate', 'sulfite', 'sulfonate','sulfate reduction', 'sulfur oxidation', 'sulfur respiration', 'SRB', 'dsrAB', 'APS reductase', 'sulfide:quinone oxidoreductase','sulfur disproportionation', 'dissimilatory sulfate reduction','PWY-6932','SO4ASSIM-PWY', 'SULFATE-CYS-PWY', 'sulfur globules', 'elemental sulfur', 'polysulfide metabolism', 'sulfur granules',   'sulfate reduction I (assimilatory)', 'superpathway of sulfate assimilation and cysteine biosynthesis' ],
'hydrogen_metabolism': ['h2_consumption', 'hydrogenase', 'hydrogen uptake', 'hydrogen consumption', 'h2', 'H2 oxidation', 'H2ase', 'hydrogen production','FeFe-hydrogenase', 'NiFe-hydrogenase', 'hydrogen evolution','hydrogen cycling', 'H2 sensing', 'proton reduction'],
'organic_acid_metabolism': ['acetate', 'acetic acid', 'acetyl', 'acetate metabolism', 'acetate production','oxalate', 'oxalic acid', 'oxalate metabolism', 'oxalate production', 'organic acid', 'fatty acid', 'butyric acid', 'butyrate', 'propionate', 'propionic acid','carboxylic acid', 'lactate', 'lactic acid', 'formate', 'formic acid','citrate', 'citric acid', 'succinate', 'succinic acid', 'fumarate', 'fumaric acid','malate', 'malic acid', 'pyruvate', 'pyruvic acid', 'CENTFERM-PWY', 'FERMENTATION-PWY', 'GLYCOLYSIS', 'PWY-5100','GALACTUROCAT-PWY', 'gluconic acid', 'tartaric acid', 'caproic acid', 'valeric acid', 'fatty acid β-oxidation I', 'fatty acid elongation', 'fatty acid salvage','stearate biosynthesis II (bacteria and plants)', 'palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)', 'cis-vaccenate biosynthesis', 'oleate biosynthesis IV (anaerobic)', 'gondoate biosynthesis (anaerobic)', 'mycolate biosynthesis' ],
'metal_organic_interaction': ['siderophore', 'metal binding', 'metal chelation', 'iron chelation', 'iron complex', 'metal transport', 'metallophore', 'metal organic', 'organometallic', 'iron uptake', 'metal uptake', 'metalloprotein', 'iron-sulfur cluster', 'metal coordination', 'metal sequestration','ferric reductase', 'ferrous oxidase', 'metal homeostasis', 'metal bioavailability', 'metallophore production', 'metal-ligand complex', 'chelator secretion', 'metal bioavailability' ],
'biofilm_formation': ['biofilm', 'exopolysaccharide', 'EPS production', 'extracellular matrix', 'adhesion', 'attachment', 'colonization', 'surface attachment', 'polysaccharide biosynthesis', 'extracellular polymeric substance', 'cell aggregation', 'quorum sensing', 'biofilm maturation','biofilm regulation', 'biofilm dispersion', 'cell-cell adhesion', 'matrix production', 'pellicle', 'floc formation', 'COLANSYN-PWY','EXOPOLYSACC-PWY', 'GLUCOSE1PMETAB-PWY', 'alginate synthesis', 'cellulose biosynthesis', 'amyloid fiber', 'biofilm architecture', 'bioelectricity',     'peptidoglycan biosynthesis I (meso-diaminopimelate containing)', 'peptidoglycan biosynthesis III (mycobacteria)','UDP-N-acetylmuramoyl-pentapeptide biosynthesis I (meso-diaminodiaminopimelate containing)', 'UDP-N-acetylmuramoyl-pentapeptide biosynthesis II (lysine-containing)', 'dTDP-L-rhamnose biosynthesis I', 'O-antigen building blocks biosynthesis (E. coli)', 'phosphatidylglycerol biosynthesis I (plastidic)', 'phosphatidylglycerol biosynthesis II (non-plastidic)', 'CDP-diacylglycerol biosynthesis I', 'CDP-diacylglycerol biosynthesis II'],
'carbon_metabolism': ['carbon fixation', 'carbon utilization', 'carbohydrate metabolism', 'glycolysis', 'TCA cycle', 'pentose phosphate pathway', 'gluconeogenesis', 'carbon assimilation', 'carbon flux','Calvin cycle', 'reductive acetyl-CoA pathway', 'carbon monoxide dehydrogenase','carbon catabolite repression', 'carbon storage', 'glycogen metabolism',     'TCA cycle IV (2-oxoglutarate decarboxylase)',  'TCA cycle V (2-oxoglutarate:ferredoxin oxidoreductase)', 'TCA cycle VI (obligate autotrophs)', 'TCA cycle VIII (helicobacter)', 'superpathway of glyoxylate bypass and TCA', 'methylerythritol phosphate pathway I', 'methylerythritol phosphate pathway II', 'glycolysis I (from glucose 6-phosphate)', 'glycolysis II (from fructose 6-phosphate)', 'glycolysis III (from glucose)',    'gluconeogenesis I', 'pentose phosphate pathway (non-oxidative branch)', 'Calvin-Benson-Bassham cycle', 'pyruvate fermentation to isobutanol (engineered)', 'superpathway of branched amino acid biosynthesis', 'superpathway of aromatic amino acid biosynthesis',  'superpathway of adenosine nucleotides de novo biosynthesis I', 'superpathway of adenosine nucleotides de novo biosynthesis II', 'superpathway of guanosine nucleotides de novo biosynthesis I',    'superpathway of guanosine nucleotides de novo biosynthesis II', 'coenzyme A biosynthesis I',  'pantothenate and coenzyme A biosynthesis I', 'phosphopantothenate biosynthesis I','NAD biosynthesis I (from aspartate)', 'thiamin salvage II',  'superpathway of tetrahydrofolate biosynthesis and salvage',  'flavin biosynthesis I (bacteria and plants)', 'chorismate biosynthesis I',    'chorismate biosynthesis from 3-dehydroquinate' ],
'methanogenesis': ['methanogenesis', 'methanobacterium', 'archaea', 'methane production', 'methyl-coenzyme M reductase', 'methanogenic', 'coenzyme F420','methyl-H4MPT', 'CO2 reduction', 'acetoclastic methanogenesis', 'hydrogenotrophic methanogenesis', 'methylotrophic methanogenesis', 'methanotrophy' ],
#'silicon_metabolism': [ 'silicon uptake', 'silicate transport', 'silicon cycling', 'diatom frustule','silicification', 'biogenic silica', 'silicon precipitation'],
'phosphorus_metabolism': ['phosphate uptake', 'phosphate transport', 'polyphosphate', 'phosphorus cycling','phosphate precipitation', 'struvite formation', 'phosphate solubilization',     'UMP biosynthesis', 'pyrimidine deoxyribonucleotides de novo biosynthesis I','pyrimidine deoxyribonucleotide phosphorylation','superpathway of pyrimidine nucleobases salvage','superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis','adenosine ribonucleotides de novo biosynthesis', 'adenosine deoxyribonucleotides de novo biosynthesis II', 'guanosine ribonucleotides de novo biosynthesis','guanosine deoxyribonucleotides de novo biosynthesis II','superpathway of pyrimidine ribonucleotides de novo biosynthesis','superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis (E. coli)']
} 

#7 Environmental factors/ Operational Factors
'halogen_related': {'terms': ['halogen', 'chloride', 'bromide', 'iodide', 'fluoride', 'halide', 'dehalogenation', 'haloperoxidase', 'haloacid', 'chlorination', 'bromination', 'organohalide', 'halomethane', 'haloalkane', 'organohalide', 'halotolerance', 'salt tolerance', 'halophilic', 'chloride transport', 'halide channel', 'chloride attack', 'chloride-induced corrosion', 'pitting initiation', 'chloride penetration', 'halide corrosion', 'perchlorate reduction', 'halorespiration', 'organohalide'],
'microaerobic_conditions': ['microaerophilic', 'oxygen-limited', 'hypoxic', 'microaerobic'],
'ph_modulation': ['acid tolerance', 'alkaline tolerance', 'proton pump', 'pH homeostasis', 'pH stress', 'pH regulation', 'acid resistance', 'pH buffering', 'pH gradient', 'acidophile', 'alkaliphile'], 
'direct_eet': ['electron transfer', 'extracellular electron transfer', 'direct electron transfer',  'solid-state respiration', 'interfacial electron transfer', '',  'electrode respiration', 'bioelectricity', 'electron shuttle', 'electron bridge','electron conduit', 'oxidoreductase', 'reductase', 'oxidase', 'redox',  'cytochrome', 'c-type cytochrome', 'multi-heme cytochrome', 'cytochrome c oxidase',  'quinol oxidase', 'NADH:quinone oxidoreductase', 'succinate dehydrogenase',  'fumarate reductase', 'mtrABC', 'omc', 'omcS', 'nanowire', 'microbial nanowires',  'conductive pili', 'type IV pili', 'conductive biofilms', 'flavin',  'riboflavin shuttle', 'flavin secretion', 'histidine kinase', 'PROTEIN-KINASE-RXN',  'L-arginine biosynthesis II', 'chorismate biosynthesis I',  'superpathway of branched amino acid biosynthesis',  'superpathway of aromatic amino acid biosynthesis',  'L-lysine biosynthesis I', 'L-ornithine biosynthesis', 'oxidoreductase activity']
'indirect_eet': ['shuttle', 'mediator', 'redox mediator', 'electron shuttle', 'flavin', 'quinone', 'humic substance', 'pyocyanin', 'phenazine', 'riboflavin'], 
'exoelectrogenesis': ['exoelectrogen', 'electrochemically active bacteria', 'EAB', 'extracellular respiration', 'electrode respiration', 'bioelectrosynthesis', 'microbial fuel cell'], 
'enzymatic_corrosion': ['enzyme-mediated corrosion', 'metalloenzyme', 'enzyme-catalyzed oxidation', 'enzyme-catalyzed reduction', 'peroxidase', 'laccase', 'oxidoreductase activity'], 
'dealloying': ['selective corrosion', 'dezincification', 'dealuminification', 'parting', 'selective leaching', 'alloy corrosion', 'preferential dissolution'], 
'galvanic_corrosion': ['bimetallic corrosion', 'galvanic couple', 'differential aeration', 'local cell formation', 'microelectrochemical cells', 'electrolytic corrosion'], 
'chloride_attack': ['chloride-induced corrosion', 'pitting initiation', 'localized corrosion', 'crevice corrosion', 'stress corrosion cracking', 'chloride penetration'], 
'microbe_metal_synergy': ['metal-microbe interaction', 'biomineralization', 'biosorption', 'bioaccumulation', 'metal-binding proteins', 'metallothionein', 'metal-enzyme complex'], 
'cathodic_depolarization': ['cathodic depolarization', 'hydrogen depolarization', 'cathodic reaction enhancement', 'electron acceptor consumption', 'cathodic process acceleration'],
'passivity_breakdown': ['passivation breakdown', 'passive film disruption', 'oxide layer dissolution', 'passive layer damage', 'depassivation'],
'concentration_cells': ['concentration cell', 'differential concentration', 'ion concentration gradient', 'oxygen concentration cell', 'metal ion gradient'],
'syntrophic_interactions': ['syntrophic', 'cross-feeding', 'metabolic coupling', 'interspecies', 'consortia', 'metabolic cooperation']
}
  