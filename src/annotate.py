"""
Annotation Module

Parses FASTA headers to extract transcript IDs, gene symbols, and transcript types.
Supports both mock format (GENE_NAME|region_type) and RefSeq format.
"""

import re


def parse_fasta_header(header: str) -> dict:
    """
    Parse a FASTA header to extract transcript ID, gene symbol, and transcript type.
    
    Supports two formats:
    1. Mock format: GENE_NAME|region_type
       Example: "GENE1|exon" -> transcript_id="GENE1|exon", gene_symbol="GENE1", transcript_type="exon"
    
    2. RefSeq format: NM_XXXXX.X Homo sapiens gene name (GENE_SYMBOL), transcript variant X, mRNA
       Example: "NM_000014.6 Homo sapiens alpha-2-macroglobulin (A2M), transcript variant 1, mRNA"
       -> transcript_id="NM_000014.6", gene_symbol="A2M", transcript_type="mRNA"
    
    Args:
        header: FASTA header line (with or without '>' prefix)
    
    Returns:
        Dictionary with keys:
            - transcript_id: RefSeq transcript ID (e.g., NM_000014.6) or full header for mock format
            - gene_symbol: Gene symbol (e.g., A2M) or gene name from mock format
            - transcript_type: Type of transcript (mRNA, non-coding RNA, lncRNA, or mock region_type)
            - full_header: Full header (for reference)
    
    Note:
        If header doesn't match expected formats, attempts to extract what it can,
        defaulting to "NA" for missing parts (not "UNKNOWN" to match CSV convention).
    """
    # Remove '>' prefix if present
    header = header.lstrip('>').strip()
    full_header = header
    
    # Try RefSeq format first (most common for real data)
    # Pattern: NM_XXXXX.X or NR_XXXXX.X or XM_XXXXX.X or XR_XXXXX.X at start
    refseq_pattern = r'^([NX][MR]_\d+\.\d+)\s+(.+)$'
    match = re.match(refseq_pattern, header)
    
    if match:
        # RefSeq format detected
        transcript_id = match.group(1)  # e.g., "NM_000014.6"
        rest_of_header = match.group(2)  # Rest of the header
        
        # Extract gene symbol from parentheses: (GENE_SYMBOL)
        gene_symbol_match = re.search(r'\(([A-Z0-9_-]+)\)', rest_of_header)
        if gene_symbol_match:
            gene_symbol = gene_symbol_match.group(1)
        else:
            gene_symbol = "NA"  # Not available
        
        # Extract transcript type (mRNA, non-coding RNA, lncRNA, etc.)
        transcript_type = "NA"
        rest_lower = rest_of_header.lower()
        if 'non-coding rna' in rest_lower or 'noncoding rna' in rest_lower:
            transcript_type = "non-coding RNA"
        elif 'lncrna' in rest_lower or 'long non-coding' in rest_lower:
            transcript_type = "lncRNA"
        elif 'mrna' in rest_lower:
            transcript_type = "mRNA"
        elif 'rna' in rest_lower:
            transcript_type = "RNA"  # Generic RNA if type unclear
        
        return {
            "transcript_id": transcript_id,
            "gene_symbol": gene_symbol,
            "transcript_type": transcript_type,
            "full_header": full_header
        }
    
    # Try mock format: GENE_NAME|region_type
    if '|' in header:
        parts = header.split('|')
        if len(parts) >= 2:
            gene_symbol = parts[0].strip()
            transcript_type = parts[1].strip()
            return {
                "transcript_id": header,  # Use full header as transcript_id for mock format
                "gene_symbol": gene_symbol,
                "transcript_type": transcript_type,
                "full_header": full_header
            }
    
    # Fallback: couldn't parse, return what we can
    # Try to extract first token as potential transcript ID
    first_token = header.split()[0] if header.split() else "NA"
    return {
        "transcript_id": first_token if first_token else "NA",
        "gene_symbol": "NA",
        "transcript_type": "NA",
        "full_header": full_header
    }


def is_functional_region(transcript_type: str) -> bool:
    """
    Check if a transcript type is considered "functional" for screening.
    
    For RefSeq transcripts: All annotated transcripts (mRNA, non-coding RNA, lncRNA)
    are considered functional since they represent transcribed regions.
    
    For mock format: Accepts exon, intron, lncRNA as before.
    
    Args:
        transcript_type: Transcript type string (e.g., "mRNA", "non-coding RNA", "lncRNA", "exon", "intron")
    
    Returns:
        True if transcript is functional, False otherwise
    
    Note:
        With RefSeq data, we screen all transcripts since they represent functional
        transcribed regions. The distinction between coding and non-coding is
        captured in the transcript_type field for downstream analysis.
    """
    if transcript_type == "NA":
        # If type is unknown, be conservative and include it
        return True
    
    transcript_type_lower = transcript_type.lower()
    
    # RefSeq transcript types (all are functional)
    refseq_functional = {"mrna", "rna", "non-coding rna", "noncoding rna", "lncrna", "long non-coding"}
    
    # Mock format types (for backward compatibility)
    mock_functional = {"exon", "intron", "lncRNA"}
    
    # Check if it matches any functional type
    if transcript_type_lower in refseq_functional:
        return True
    
    if transcript_type_lower in mock_functional:
        return True
    
    # Check for partial matches (e.g., "non-coding RNA" contains "rna")
    if any(func_type in transcript_type_lower for func_type in refseq_functional):
        return True
    
    # Default: be conservative and include unknown types
    return True

