"""
Scanning Module

Implements sliding window approach to scan transcripts for off-target hits.
Uses substitution-only distance (no indels) - only compares windows of exactly ASO length.
"""

from .edit_distance import edit_distance, is_valid_hit
from .annotate import parse_fasta_header, is_functional_region


def scan_transcript(aso_id: str, aso_sequence: str, transcript_header: str,
                    transcript_sequence: str, max_edit_distance: int = 2) -> list:
    """
    Scan a single transcript for potential off-target hits of an ASO.
    
    Uses a sliding window approach: slides the ASO length across the transcript
    and computes substitution-only distance (mismatch count) at each position.
    Only windows of exactly ASO length are compared - insertions and deletions
    are ignored as they strongly disrupt ASO hybridization.
    
    Args:
        aso_id: Identifier for the ASO
        aso_sequence: ASO nucleotide sequence
        transcript_header: FASTA header of the transcript
        transcript_sequence: Transcript nucleotide sequence
        max_edit_distance: Maximum mismatches (substitution-only) to consider a hit (default: 2)
    
    Returns:
        List of hit dictionaries, each containing:
            - aso_id
            - aso_sequence
            - transcript_id
            - gene_symbol
            - transcript_type
            - match_start (0-indexed)
            - match_end (0-indexed, exclusive)
            - matched_sequence
            - edit_distance (mismatch count, substitution-only)
    
    Performance Note:
        For very long transcripts (>100kb), this sliding window approach may be slow.
        The substitution-only distance computation is O(n) where n=ASO length.
        With ASO length ~20bp, this is very fast even for long transcripts.
    """
    hits = []
    
    # Parse transcript annotation from header
    annotation = parse_fasta_header(transcript_header)
    gene_symbol = annotation["gene_symbol"]
    transcript_type = annotation["transcript_type"]
    transcript_id = annotation["transcript_id"]
    
    # Only scan functional regions
    if not is_functional_region(transcript_type):
        return hits
    
    aso_len = len(aso_sequence)
    transcript_len = len(transcript_sequence)
    
    # Safety check: skip transcripts shorter than ASO length
    # Substitution-only distance requires equal-length sequences
    if transcript_len < aso_len:
        return hits
    
    # Slide window across transcript
    # For each position, extract a window of exactly ASO length and compute
    # substitution-only distance (mismatch count). Only equal-length windows
    # are compared - indels are ignored as they disrupt ASO hybridization.
    for pos in range(transcript_len - aso_len + 1):
        window = transcript_sequence[pos:pos + aso_len]
        
        # Compute substitution-only distance (mismatch count)
        # Window is guaranteed to be same length as ASO (aso_len)
        dist = edit_distance(aso_sequence, window)
        
        # Check if it's a valid hit (dist >= 0 ensures equal-length sequences)
        if is_valid_hit(dist, max_edit_distance):
            match_end = pos + aso_len  # Exclusive end position
            hits.append({
                "aso_id": aso_id,
                "aso_sequence": aso_sequence,
                "transcript_id": transcript_id,
                "gene_symbol": gene_symbol,
                "transcript_type": transcript_type,
                "match_start": pos,
                "match_end": match_end,
                "matched_sequence": window,
                "edit_distance": dist
            })
    
    return hits


def scan_all_sequences(aso_list: list, transcript_list: list,
                       max_edit_distance: int = 2) -> list:
    """
    Scan all ASOs against all transcripts using substitution-only distance.
    
    Args:
        aso_list: List of (aso_id, sequence) tuples
        transcript_list: List of (header, sequence) tuples
        max_edit_distance: Maximum mismatches (substitution-only) to consider a hit (default: 2)
    
    Returns:
        List of all hit dictionaries from all ASO-transcript pairs
    
    Performance Notes:
        - Memory: Each transcript sequence is held in memory. For ~200k transcripts
          with average length ~2kb, this is ~400MB, which is acceptable.
        - Time: O(ASOs × transcripts × transcript_length × ASO_length)
          With ~40k ASOs, ~200k transcripts, ~2kb avg length, ~20bp ASOs:
          ~40k × 200k × 2k × 20 = 3.2×10^14 operations (theoretical worst case)
          In practice, substitution-only distance computation is O(ASO_length) per window,
          which is very fast (~20 operations per comparison).
        - Very long transcripts (>100kb) may be slow but are handled.
    """
    all_hits = []
    
    # Performance safety: Check for very long transcripts and report
    max_transcript_len = 0
    long_transcript_count = 0
    for _, transcript_sequence in transcript_list:
        tlen = len(transcript_sequence)
        if tlen > max_transcript_len:
            max_transcript_len = tlen
        if tlen > 100000:  # >100kb
            long_transcript_count += 1
    
    print(f"Scanning {len(aso_list)} ASOs against {len(transcript_list)} transcripts...")
    if max_transcript_len > 0:
        print(f"  Transcript length range: max={max_transcript_len:,} bp")
    if long_transcript_count > 0:
        print(f"  Warning: {long_transcript_count} transcripts >100kb (may be slow)")
    
    total_comparisons = len(aso_list) * len(transcript_list)
    if total_comparisons > 1000000:
        print(f"  Note: {total_comparisons:,} total ASO-transcript comparisons")
    
    for aso_idx, (aso_id, aso_sequence) in enumerate(aso_list, 1):
        print(f"  Scanning ASO {aso_idx}/{len(aso_list)}: {aso_id}")
        
        for transcript_header, transcript_sequence in transcript_list:
            hits = scan_transcript(
                aso_id, aso_sequence, transcript_header,
                transcript_sequence, max_edit_distance
            )
            all_hits.extend(hits)
    
    print(f"Found {len(all_hits)} total hits.")
    
    return all_hits

