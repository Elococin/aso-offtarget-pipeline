"""
Edit Distance Computation Module

Implements substitution-only distance (Hamming distance) for ASO off-target screening.
This measures mismatches between equal-length sequences only - insertions and deletions
are ignored because they strongly disrupt ASO hybridization.
"""


def edit_distance(seq1: str, seq2: str, threshold: int = None) -> int:
    """
    Compute substitution-only distance (Hamming distance) between two sequences.

    This function counts only mismatches (substitutions) between equal-length sequences.
    Insertions and deletions are not considered, as they strongly disrupt ASO hybridization.

    Args:
        seq1: First sequence (typically the ASO, should be pre-uppercased)
        seq2: Second sequence (typically a transcript window, should be pre-uppercased)
        threshold: Optional early termination threshold. If provided, stops counting
                   once mismatches exceed this value (returns early for efficiency).

    Returns:
        Integer mismatch count (0 = identical, higher = more mismatches)
        Returns -1 if sequences have different lengths (invalid for substitution-only distance)

    Note:
        For performance, caller should pre-uppercase sequences before calling this
        function in a loop. This avoids redundant .upper() calls.

    Example:
        >>> edit_distance("ATCG", "ATCG")
        0
        >>> edit_distance("ATCG", "ATCC")
        1
        >>> edit_distance("ATCG", "ATC")
        -1  # Different lengths - not valid for substitution-only distance
    """
    # Substitution-only distance requires equal-length sequences
    if len(seq1) != len(seq2):
        return -1  # Invalid - sequences must be equal length

    # Count mismatches using zip (faster than indexing)
    mismatches = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != c2:
            mismatches += 1
            # Early termination: stop if we've exceeded the threshold
            if threshold is not None and mismatches > threshold:
                return mismatches

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

