"""
Custom exceptions for the corrosion scoring system.
"""

class ScoringError(Exception):
    """Base exception for scoring system errors."""
    pass

class TextMiningError(ScoringError):
    """Exception raised during text mining operations."""
    pass

class SynergyDetectionError(ScoringError):
    """Exception raised during synergy detection."""
    pass

class ProcessorError(ScoringError):
    """Exception raised when processors are invalid or missing."""
    pass

class ValidationError(ScoringError):
    """Exception raised when input validation fails."""
    pass