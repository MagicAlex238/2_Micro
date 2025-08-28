import re
import unicodedata

from .utils_ec import normalize_ec_id, strip_all_ec_tokens

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

custom_dict = {
    "bifunctional glutamine-synthetase adenylyltransferase-adenylyl-removing enzyme": "bifunctional glutamine-synthetase adenylyltransferase",
    "glutamine--fructose-6-phosphate aminotransferase": "glutamine--fructose-6-phosphate aminotransferase",
    "udp-3-o-acyl-n-acetylglucosamine deacetylase": "udp-3-o-acyl-n-acetylglucosamine deacetylase",
    "multifunctional oxoglutarate decarboxylase": "decarboxylase-oxoglutarate-dehydrogenase thiamine pyrophosphate",
    "aldo-keto-reductase 2 family": "aldo-keto-reductase", 
    "coenzyme a biosynthesis bifunctional protein coabc (dna-pantothenate metabolism flavoprotein) (phosphopantothenoylcysteine-synthetase-decarboxy": "coenzyme a biosynthesis protein ppcs-ppcdc"
}

def clean_protein_name(name: str) -> str:
    """Simplify the names of the protein without losing the biological meaning as literature recommends."""
    if not isinstance(name, str):
        name = '' if name is None else str(name)
    # normalise name
    name = name.strip()
    name_lower = name.lower()

    # Step 2: Remove uncertainty terms at the beginning
    uncertainty_terms = ['probable', 'putative', 'possible', 'uncharacterized', 'hypothetical']
    pattern_uncertainty = r'^(?:' + '|'.join(uncertainty_terms) + r')\s+'
    name = re.sub(pattern_uncertainty, '', name, flags=re.IGNORECASE)

    # Step 3: Apply custom dictionary based on prefix match
    for original_start, short_name in custom_dict.items():
        if name_lower.startswith(original_start.lower()):
            return short_name.lower()

    # step 5 cut at ~50 chars but finish the current word
    if len(name) > 50:
        tail = name[50:]
        m = re.search(r'[\s\.\-\[\]\(\)]', tail)
        end = 50 + (m.start() if m else len(name))
        name = name[:end].strip()

    return name
