import re
from corrosion_scoring.utils_ec import normalize_ec_id
from corrosion_scoring.name_utils import enhanced_clean_protein_name

def read_enzyme_names(unique_ecs_to_filter=None):
    paths = setup_paths()
    enzyme_path = paths['enzyme']

    # Normalize filter set
    ec_filter_norm = None
    if unique_ecs_to_filter is not None:
        ec_filter_norm = {normalize_ec_id(ec) for ec in unique_ecs_to_filter}
        ec_filter_norm.discard(None)

    ec_to_names = {}
    with open(enzyme_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t', 1)
            if len(parts) < 2:
                continue
            raw_ec, raw_names = parts[0], parts[1]
            ec_norm = normalize_ec_id(raw_ec)
            if not ec_norm:
                continue
            if ec_filter_norm is not None and ec_norm not in ec_filter_norm:
                continue
            names = [n.strip() for n in re.split(r';\s*', raw_names) if n.strip()]
            cleaned = [enhanced_clean_protein_name(n) for n in names]
            ec_to_names[ec_norm] = cleaned
    return ec_to_names