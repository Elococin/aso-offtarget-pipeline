"""
Input/Output Module

Handles reading ASO sequences and FASTA transcript files.
"""


def read_aso_sequences(filepath: str) -> list:
    """
    Read ASO sequences from a text file.
    
    Expected format: ASO_ID SEQUENCE (space-separated, one per line)
    Example:
        ASO_001 ATCGATCGATCGATCGATCG
        ASO_002 GCTAGCTAGCTAGCTAGCTA
    
    Args:
        filepath: Path to the ASO sequences file
    
    Returns:
        List of tuples: [(aso_id, sequence), ...]
    
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    aso_list = []
    
    try:
        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                # Split by whitespace (allows multiple spaces or tabs)
                parts = line.split()
                
                if len(parts) < 2:
                    raise ValueError(
                        f"Invalid format at line {line_num}: expected 'ASO_ID SEQUENCE', "
                        f"got: {line}"
                    )
                
                aso_id = parts[0]
                sequence = ''.join(parts[1:])  # Join in case sequence has spaces
                
                # Validate sequence contains only valid nucleotides
                valid_bases = set('ATCGN')  # N for ambiguous
                if not all(base.upper() in valid_bases for base in sequence):
                    raise ValueError(
                        f"Invalid nucleotide at line {line_num}: "
                        f"sequence contains non-ATCGN characters"
                    )
                
                aso_list.append((aso_id, sequence.upper()))
    
    except FileNotFoundError:
        raise FileNotFoundError(f"ASO sequences file not found: {filepath}")
    
    if not aso_list:
        raise ValueError(f"No valid ASO sequences found in {filepath}")
    
    return aso_list


def read_fasta(filepath: str) -> list:
    """
    Read sequences from a FASTA file.
    
    Expected format:
        >HEADER1
        SEQUENCE1
        >HEADER2
        SEQUENCE2
    
    Args:
        filepath: Path to the FASTA file
    
    Returns:
        List of tuples: [(header, sequence), ...]
        Headers are stored without the '>' prefix
    
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    sequences = []
    current_header = None
    current_sequence = []
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_header is not None:
                        sequences.append((
                            current_header,
                            ''.join(current_sequence).upper()
                        ))
                    
                    # Start new sequence
                    current_header = line[1:].strip()  # Remove '>'
                    current_sequence = []
                else:
                    # Append to current sequence
                    current_sequence.append(line)
            
            # Don't forget the last sequence
            if current_header is not None:
                sequences.append((
                    current_header,
                    ''.join(current_sequence).upper()
                ))
    
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    
    if not sequences:
        raise ValueError(f"No sequences found in {filepath}")
    
    return sequences


def write_results_csv(results: list, filepath: str):
    """
    Write screening results to a CSV file.
    
    Args:
        results: List of dictionaries, each containing:
            - aso_id
            - aso_sequence
            - transcript_id
            - gene_symbol
            - transcript_type
            - match_start
            - match_end
            - matched_sequence
            - edit_distance
        filepath: Path to output CSV file
    
    Output Schema:
        The CSV contains the following columns (in order):
        - aso_id: ASO identifier
        - aso_sequence: ASO nucleotide sequence
        - transcript_id: RefSeq transcript ID (e.g., NM_000014.6)
        - gene_symbol: Gene symbol (e.g., A2M) or "NA" if not available
        - transcript_type: Transcript type (mRNA, non-coding RNA, lncRNA) or "NA"
        - match_start: Start position of match in transcript (0-indexed)
        - match_end: End position of match in transcript (0-indexed, exclusive)
        - matched_sequence: Actual sequence from transcript that matched
        - edit_distance: Computed edit distance (0-2)
    """
    import csv
    
    # Define fieldnames in desired order for human readability
    fieldnames = [
        'aso_id',
        'aso_sequence',
        'transcript_id',
        'gene_symbol',
        'transcript_type',
        'match_start',
        'match_end',
        'matched_sequence',
        'edit_distance'
    ]
    
    if not results:
        # Create empty file with headers
        with open(filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
        return
    
    # Write results
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

