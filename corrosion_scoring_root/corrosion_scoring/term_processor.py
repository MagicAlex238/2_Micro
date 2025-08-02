import re
from collections import OrderedDict

class TermProcessor:
    """Normalizes and matches corrosion terms with priority handling"""
    
    def __init__(self, taxonomy: dict):
        self.priority_order = ['iron_metabolism', 'sulfur_metabolism', 'organic_acid_metabolism',
        'biofilm_formation', 'o2_consumption', 'metal_binding_chelation',
        'nitrogen_metabolism', 'manganese_processes', 'h2_consumption',
        'methanogenesis', 'carbon_metabolism',  'fumarate_formation',
        'phosphorus_metabolism']
        
        self.normalized_taxonomy = self._create_priority_dict(taxonomy)

    def _normalize_term(self, term: str) -> str:
        """Enhanced normalization with corrosion-specific substitutions"""
        substitutions = {
            'reduc': 'reduction',
            'oxid': 'oxidation', 
            'ferri': 'iron',
            'ferro': 'iron',
            'sulph': 'sulf',
            'metallo': 'metal',
            'corrosi': 'corrosion'
        }
        
        base_term = re.sub(r'[\W_]+', '', term).lower()
        for pattern, replacement in substitutions.items():
            base_term = re.sub(f'^{pattern}', replacement, base_term)
        return base_term

    def matches_normalized(self, term: str, text: str) -> bool:
        """Check if normalized term appears in text"""
        norm_term = self._normalize_term(term)
        norm_text = self._normalize_term(text)
        return norm_term in norm_text
    
    def _create_priority_dict(self, taxonomy: dict) -> OrderedDict:
        """Create search priority structure"""
        priority_dict = OrderedDict()
        
        for cat in self.priority_order:
            if cat in taxonomy:
                if isinstance(taxonomy[cat], dict) and 'terms' in taxonomy[cat]:
                    # Handle functional_categories format
                    priority_dict[cat] = {self._normalize_term(t) for t in taxonomy[cat]['terms']}
                elif isinstance(taxonomy[cat], list):
                    # Handle simple list format
                    priority_dict[cat] = {self._normalize_term(t) for t in taxonomy[cat]}
        
        # Add any remaining categories not in priority order
        for cat, terms in taxonomy.items():
            if cat not in priority_dict:
                if isinstance(terms, dict) and 'terms' in terms:
                    priority_dict[cat] = {self._normalize_term(t) for t in terms['terms']}
                elif isinstance(terms, list):
                    priority_dict[cat] = {self._normalize_term(t) for t in terms}
        
        return priority_dict
    
    def find_first_category(self, term: str) -> str:
        """Returns first matching category based on priority"""
        norm_term = self._normalize_term(term)
        for category, terms in self.normalized_taxonomy.items():
            if norm_term in terms:
                return category
        return None
#=================================================================================
def enhance_pathway_extraction(record, ec_pathway_mapping, ipath_mapping, ko_ec):
    """
    Enhanced pathway extraction that better integrates all sources and normalizes terms.
    """
    ec_number = record['ec_number']
    
    # Start with existing pathways
    all_pathways = set(record.get('pathways', []))
    
    # 1. Add pathways from EC-pathway mapping with better processing
    if ec_number in ec_pathway_mapping:
        for pathway_id in ec_pathway_mapping[ec_number]:
            # Standardize pathway ID
            std_id = pathway_id
            if pathway_id.startswith('ec'):
                std_id = 'map' + pathway_id[2:]
            
            # Look up pathway name and normalize
            if std_id in pathway_data:
                pathway_name = pathway_data[std_id]
                # Use pathway processor to normalize
                normalized_pathway = _pathway_processor._normalize_term(pathway_name)
                all_pathways.add(pathway_name)  # Keep original for display
    
    # 2. Add ipath data (ground truth) - this have priority
    if ec_number in ipath_mapping:
        for ipath_pathway in ipath_mapping[ec_number]:
            # Normalize ipath terms using pathway processor
            if isinstance(ipath_pathway, str):
                normalized_ipath = _pathway_processor._normalize_term(ipath_pathway)
                all_pathways.add(ipath_pathway)
    
    # 3. Add pathways from KO data with normalization
    if ec_number in ko_ec:
        ko_data = ko_ec[ec_number]
        if isinstance(ko_data, list):
            for path in ko_data:
                if isinstance(path, str):
                    all_pathways.add(path)
        elif isinstance(ko_data, dict) and 'pathway' in ko_data:
            path = ko_data['pathway']
            if isinstance(path, str):
                all_pathways.add(path)
    
    # 4. Extract mechanisms from all collected pathways
    pathway_text = ' '.join(all_pathways)
    extracted_mechanisms = assign_mechanism_from_pathway(pathway_text)
    
    # Update record
    record['pathways'] = list(all_pathways)
    existing_mechanisms = record.get('corrosion_mechanisms', [])
    record['corrosion_mechanisms'] = list(set(existing_mechanisms + extracted_mechanisms))
    
    return record

def validate_against_ipath(record, ipath_data):
    """
    Validates detected pathways, mechanisms, and functional categories against ipath ground truth.
    Returns validation metrics and suggestions.
    """
    ec_number = record['ec_number']
    validation_results = {
        'pathway_validation': {},
        'mechanism_validation': {},
        'functional_category_validation': {},
        'overall_confidence': 0.0
    }
    
    if ec_number not in ipath_data:
        validation_results['overall_confidence'] = 0.5  # No ground truth available
        return validation_results
    
    ipath_pathways = ipath_data[ec_number]
    detected_pathways = record.get('pathways', [])
    detected_mechanisms = record.get('corrosion_mechanisms', [])
    detected_fc = [fc['category'] for fc in record.get('functional_categories', [])]
    
    # Pathway validation
    if ipath_pathways and detected_pathways:
        # Normalize both for comparison
        norm_ipath = {_pathway_processor._normalize_term(p) for p in ipath_pathways}
        norm_detected = {_pathway_processor._normalize_term(p) for p in detected_pathways}
        
        overlap = norm_ipath.intersection(norm_detected)
        pathway_precision = len(overlap) / len(norm_detected) if norm_detected else 0
        pathway_recall = len(overlap) / len(norm_ipath) if norm_ipath else 0
        
        validation_results['pathway_validation'] = {
            'precision': pathway_precision,
            'recall': pathway_recall,
            'f1_score': 2 * (pathway_precision * pathway_recall) / (pathway_precision + pathway_recall) if (pathway_precision + pathway_recall) > 0 else 0,
            'overlap_terms': list(overlap)
        }
    
    # Mechanism validation (infer expected mechanisms from ipath)
    if ipath_pathways:
        expected_mechanisms = set()
        for pathway in ipath_pathways:
            expected_mechanisms.update(assign_mechanism_from_pathway(pathway))
        
        if expected_mechanisms and detected_mechanisms:
            norm_expected = {_mechanism_term_processor._normalize_term(m) for m in expected_mechanisms}
            norm_detected = {_mechanism_term_processor._normalize_term(m) for m in detected_mechanisms}
            
            overlap = norm_expected.intersection(norm_detected)
            mech_precision = len(overlap) / len(norm_detected) if norm_detected else 0
            mech_recall = len(overlap) / len(norm_expected) if norm_expected else 0
            
            validation_results['mechanism_validation'] = {
                'precision': mech_precision,
                'recall': mech_recall,
                'f1_score': 2 * (mech_precision * mech_recall) / (mech_precision + mech_recall) if (mech_precision + mech_recall) > 0 else 0,
                'expected_mechanisms': list(expected_mechanisms),
                'overlap_terms': list(overlap)
            }
    
    # Calculate overall confidence
    pathway_f1 = validation_results['pathway_validation'].get('f1_score', 0)
    mechanism_f1 = validation_results['mechanism_validation'].get('f1_score', 0)
    validation_results['overall_confidence'] = (pathway_f1 + mechanism_f1) / 2
    
    return validation_results