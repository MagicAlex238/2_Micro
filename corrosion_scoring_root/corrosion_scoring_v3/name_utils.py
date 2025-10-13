import re
import unicodedata

from .utils_ec import normalize_ec_id, strip_all_ec_tokens

GARBAGE_PREFIXES = (
    'transferred to', 'transferred to and', 'deleted', 'obsolete',
    'reclassified as', 'renamed to', 'see '
)

def enhanced_clean_protein_name(name: str) -> str:
    """Canonicalize DE/protein names for grouping/aggregation."""
    if name is None:
        return "uncharacterized protein"
    if not isinstance(name, str):
        name = str(name)

    # Unicode normalize (e.g., fancy dashes)
    name = unicodedata.normalize("NFKD", name)

    # Remove EC tokens (keep one in case we need a fallback)
    ec_token = normalize_ec_id(name)
    name = strip_all_ec_tokens(name)

    # Remove bracketed qualifiers
    name = re.sub(r'\([^)]*\)', '', name)
    name = re.sub(r'\[[^\]]*\]', '', name)

    # Normalize separators and punctuation
    name = re.sub(r'\s*[-/]\s*', '-', name)
    name = re.sub(r'[,:;|]+', ' ', name)
    name = re.sub(r'\s+', ' ', name).strip()

    # Lowercase
    name = name.lower()

    # Standardize common patterns
    replacements = [
        (r'\bacyl carrier protein\b', 'acp'),
        (r'\balcohol\s+dehydrogenase\b', 'alcohol-dehydrogenase'),
        (r'\bglutathione\s+dehydrogenase\b', 'glutathione-dehydrogenase'),
        (r'\b(?:l-)?threonine\s+dehydrogenase\b', 'threonine-dehydrogenase'),
        (r'\b3\-oxoacyl\-acp\s+reductase\b', '3-oxoacyl-acp-reductase'),
        (r'\b(\w+)\s+dehydrogenase\b', r'\1-dehydrogenase'),
        (r'\b(\w+)\s+reductase\b', r'\1-reductase'),
        (r'\b(\w+)\s+synthase\b', r'\1-synthase'),
        (r'\b(\w+)\s+synthetase\b', r'\1-synthetase'),
    ]
    for pat, repl in replacements:
        name = re.sub(pat, repl, name)

    # Light suffix cleanup
    name = re.sub(r'\s+(domain|fragment|precursor)$', '', name)

    # Remove immediate duplicate tokens only
    tokens = name.split()
    dedup = []
    prev = None
    for t in tokens:
        if t != prev:
            dedup.append(t)
        prev = t
    name = ' '.join(dedup)

    # Final cleanup
    name = re.sub(r'-{2,}', '-', name).strip(' -')

    if not name and ec_token:
        name = f"ec {ec_token}"
    if not name:
        name = "uncharacterized protein"
    return name
#===================================

custom_dict = {
    "bifunctional glutamine-synthetase adenylyltransferase-adenylyl-removing enzyme": "bifunctional glutamine-synthetase adenylyltransferase",
    "glutamine--fructose-6-phosphate aminotransferase": "glutamine--fructose-6-phosphate aminotransferase",
    "udp-3-o-acyl-n-acetylglucosamine deacetylase": "udp-3-o-acyl-n-acetylglucosamine deacetylase",
    "multifunctional oxoglutarate decarboxylase": "decarboxylase-oxoglutarate-dehydrogenase thiamine pyrophosphate",
    "aldo-keto-reductase 2 family": "aldo-keto-reductase", 
    "coenzyme a biosynthesis bifunctional protein coabc (dna-pantothenate metabolism flavoprotein) (phosphopantothenoylcysteine-synthetase-decarboxy": "coenzyme a biosynthesis protein ppcs-ppcdc",
    "bacterial ferritin": "bacterioferritin",
    "bacterial non-heme ferritin": "bacterioferritin",
    "bacterioferritin comigratory protein": "bacterioferritin",
    "ferredoxin-nad+ reductase": "ferredoxin-nad-reductase",
    "formate acetyltransferase": "formate c-acetyltransferase"
}

def clean_protein_name(name: str) -> str:
    """Simplify the names of the protein without losing the biological meaning as literature recommends."""
    if not isinstance(name, str):
        name = '' if name is None else str(name)
    # normalise name
    name = name.strip()
    name_lower = name.lower()
    #drop garbage prefixes early
    if any(name_lower.startswith(p) for p in GARBAGE_PREFIXES):
        return ""
    # Step 2: Remove uncertainty terms at the beginning
    uncertainty_terms = list(GARBAGE_PREFIXES) + ['probable','putative','possible','uncharacterized','hypothetical']
    pattern_uncertainty = r'^(?:' + '|'.join(uncertainty_terms) + r')\s+'
    name = re.sub(pattern_uncertainty, '', name, flags=re.IGNORECASE)

    # Step 3: Apply custom dictionary based on prefix match
    for original_start, short_name in custom_dict.items():
        if name_lower.startswith(original_start.lower()):
            name = short_name
            name_lower = name.lower()
            break

    # step 5 cut at ~50 chars but finish the current word
    if len(name) > 50:
        tail = name[50:]
        m = re.search(r'[\s\.\-\[\]\(\)]', tail)
        if m:
            end = 50 + (m.start()) 
        else:               
            end = min(50, len(name))  # hard cut at 60 if no separator found
        name = name[:end].strip()

    return name

