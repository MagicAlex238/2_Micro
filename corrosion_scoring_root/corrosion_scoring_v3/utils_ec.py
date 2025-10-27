import re
from typing import Optional
#====================== Normalizing EC =====================

# Accept "EC 1.1.1.1", "1.1.1.1", and hyphens "1.1.1.-"
_EC_RE = re.compile(r'\b(?:EC\s*)?((?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-))\b', re.IGNORECASE)

def normalize_ec_id(s: str) -> Optional[str]:
    """Return normalized EC id 'x.x.x.x' (digits or '-') or None if not found."""
    if not isinstance(s, str):
        return None
    m = _EC_RE.search(s.strip())
    return m.group(1) if m else None

def strip_all_ec_tokens(text: str) -> str:
    """Remove all EC tokens from text."""
    if not isinstance(text, str):
        return ""
    return _EC_RE.sub("", text).strip()

#====================== Normalizing normalize_listlike =====================
import re, ast, numpy as np, pandas as pd

def normalize_listlike(val):
    """
    Return a clean list[str] from mixed inputs:
    - list/tuple/set
    - numpy/pandas arrays
    - semicolon/comma-separated strings
    - stringified lists/dicts: "['Fe','Ni']" or numpy repr "['Fe' 'Ni']"
    - None/NaN -> []
    Deduplicates while preserving order.
    """
    def unique_preserve(seq):
        seen, out = set(), []
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    # None/NaN
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return []

    # list-like
    if isinstance(val, (list, tuple, set)):
        return unique_preserve([str(x).strip() for x in val if str(x).strip()])

    # numpy/pandas arrays
    if isinstance(val, (np.ndarray, pd.api.extensions.ExtensionArray, pd.Series)):
        return unique_preserve([str(x).strip() for x in list(np.array(val).tolist()) if str(x).strip()])

    # strings
    if isinstance(val, str):
        s = val.strip()
        if not s or s.lower() in {"nan", "none", "[]"}:
            return []
        # try literal eval for normal Python reprs
        if (s.startswith('[') and s.endswith(']')) or (s.startswith('{') and s.endswith('}')):
            try:
                parsed = ast.literal_eval(s)
                return normalize_listlike(parsed)
            except Exception:
                # numpy-like repr with quoted tokens and no commas
                quoted = re.findall(r"'([^']+)'|\"([^\"]+)\"", s)
                if quoted:
                    return unique_preserve([(a or b).strip() for a, b in quoted if (a or b).strip()])
        # semicolon/comma separated fallback
        if ';' in s or ',' in s:
            return unique_preserve([p.strip() for p in re.split(r'[;,]', s) if p.strip()])
        # single token string
        return [s]

    # anything else
    s = str(val).strip()
    return [s] if s else []