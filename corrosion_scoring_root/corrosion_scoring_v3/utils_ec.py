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
