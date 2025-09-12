"""
Input validation utilities for the scoring system.
"""

import re
from typing import Any, Dict, List, Optional, Union
from .exceptions import ValidationError

def validate_text(text: Any, max_length: int = 10000, allow_empty: bool = False) -> str:
    """
    Validate and sanitize input text.
    
    Args:
        text: Input text to validate
        max_length: Maximum allowed text length
        allow_empty: Whether to allow empty text
        
    Returns:
        Sanitized text string
        
    Raises:
        ValidationError: If text is invalid
    """
    if text is None:
        if allow_empty:
            return ""
        raise ValidationError("Text cannot be None")
    
    if not isinstance(text, str):
        try:
            text = str(text)
        except Exception as e:
            raise ValidationError(f"Cannot convert text to string: {e}") from e
    
    if not allow_empty and not text.strip():
        raise ValidationError("Text cannot be empty")
    
    if len(text) > max_length:
        text = text[:max_length]
    
    # Sanitize text - remove potentially problematic characters for regex
    text = re.sub(r'[\\^$.*+?{}[\]|()\x00-\x1f]', ' ', text)
    text = re.sub(r'\s+', ' ', text).strip()
    
    return text

def validate_processor(processor: Any, processor_name: str) -> None:
    """
    Validate that a processor has required methods.
    
    Args:
        processor: Processor object to validate
        processor_name: Name of the processor for error messages
        
    Raises:
        ValidationError: If processor is invalid
    """
    if processor is None:
        raise ValidationError(f"{processor_name} processor cannot be None")
    
    required_methods = ['find_all_matches', 'matches_normalized']
    for method in required_methods:
        if not hasattr(processor, method):
            raise ValidationError(f"{processor_name} processor missing required method: {method}")

def validate_processors(processors: Dict[str, Any]) -> None:
    """
    Validate all required processors are present and valid.
    
    Args:
        processors: Dictionary of processors
        
    Raises:
        ValidationError: If processors are invalid
    """
    required_processors = ['fc_processor', 'metal_processor', 'synergy_processor']
    
    if not isinstance(processors, dict):
        raise ValidationError("Processors must be provided as a dictionary")
    
    for proc_name in required_processors:
        if proc_name not in processors:
            raise ValidationError(f"Missing required processor: {proc_name}")
        validate_processor(processors[proc_name], proc_name)

def validate_metals_list(metals: Any) -> List[str]:
    """
    Validate and clean metals list.
    
    Args:
        metals: Input metals list
        
    Returns:
        Cleaned list of metal strings
    """
    if metals is None:
        return []
    
    if not isinstance(metals, (list, tuple)):
        metals = [metals]
    
    cleaned_metals = []
    for metal in metals:
        if metal is not None:
            metal_str = str(metal).strip()
            if metal_str and metal_str.lower() not in ['not detected', 'none', '']:
                cleaned_metals.append(metal_str)
    
    return cleaned_metals