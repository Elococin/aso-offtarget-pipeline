"""
Edit Distance Computation Module

Implements substitution-only distance (Hamming distance) for ASO off-target screening.
This measures mismatches between equal-length sequences only - insertions and deletions
are ignored because they strongly disrupt ASO hybridization.
"""


def edit_distance(seq1: str, seq2: str) -> int:
    """
    Compute substitution-only distance (Hamming distance) between two sequences.
    
    This function counts only mismatches (substitutions) between equal-length sequences.
    Insertions and deletions are not considered, as they strongly disrupt ASO hybridization.
    
    Args:
        seq1: First sequence (typically the ASO)
        seq2: Second sequence (typically a transcript window, must be same length as seq1)
    
    Returns:
        Integer mismatch count (0 = identical, higher = more mismatches)
        Returns -1 if sequences have different lengths (invalid for substitution-only distance)
    
    Example:
        >>> edit_distance("ATCG", "ATCG")
        0
        >>> edit_distance("ATCG", "ATCC")
        1
        >>> edit_distance("ATCG", "ATC")
        -1  # Different lengths - not valid for substitution-only distance
    """
    # Convert to uppercase for consistency
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Substitution-only distance requires equal-length sequences
    if len(seq1) != len(seq2):
        return -1  # Invalid - sequences must be equal length
    
    # Count mismatches (substitutions only)
    mismatches = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatches += 1
    
    return mismatches


def is_valid_hit(edit_dist: int, threshold: int = 2) -> bool:
    """
    Determine if a substitution-only distance represents a valid off-target hit.
    
    Args:
        edit_dist: Computed mismatch count (substitution-only distance)
        threshold: Maximum mismatches to consider a hit (default: 2)
    
    Returns:
        True if edit_dist <= threshold and edit_dist >= 0 (valid), False otherwise
    """
    # Reject invalid distances (different-length sequences return -1)
    if edit_dist < 0:
        return False
    return edit_dist <= threshold

