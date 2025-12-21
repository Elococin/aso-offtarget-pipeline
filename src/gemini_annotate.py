"""
Gemini-based Biological Consequence Annotation Module

Adds conservative biological consequence annotations to ASO off-target hits
using Google's Gemini API. Assumes monoallelic disruption and uses cautious language.
"""

import os
from google import genai
from typing import Optional, Dict, Any


def create_gemini_client() -> Optional[genai.Client]:
    """
    Create and return a Gemini client, or None if API key is missing.
    
    Returns:
        genai.Client instance if API key is available, None otherwise
    """
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        return None
    
    try:
        return genai.Client(api_key=api_key)
    except Exception:
        return None


def generate_consequence_annotation(
    client: genai.Client,
    gene_symbol: str,
    transcript_id: str,
    transcript_type: str,
    edit_distance: int,
    model_name: str = "models/gemini-2.5-flash"
) -> str:
    """
    Generate a conservative biological consequence annotation using Gemini.
    
    Args:
        client: Initialized Gemini client
        gene_symbol: Gene symbol (e.g., "A2M")
        transcript_id: Transcript ID (e.g., "NM_000014.6")
        transcript_type: Transcript type (e.g., "mRNA", "non-coding RNA")
        edit_distance: Edit distance of the hit (0-2)
        model_name: Gemini model to use (default: "gemini-1.5-flash")
    
    Returns:
        Conservative sentence describing potential biological consequences, or "NA" on failure
    """
    # Build prompt with expert role and structured guidelines
    allelic_impact = "monoallelic (assumed)"
    
    prompt = (
        f"You are a cautious computational biology assistant supporting an early-stage "
        f"in silico safety screen for antisense oligonucleotide (ASO) design.\n\n"
        f"Given the following information:\n"
        f"- Gene: {gene_symbol}\n"
        f"- Transcript type: {transcript_type}\n"
        f"- Assumed allelic impact: {allelic_impact}\n\n"
        f"Task:\n"
        f"Write ONE conservative sentence (~20–25 words) describing the potential "
        f"biological or functional consequences of disrupting this gene due to an "
        f"off-target ASO hit.\n\n"
        f"Guidelines:\n"
        f"- Use cautious language (e.g., \"may\", \"could\", \"potentially\").\n"
        f"- Focus on general biological role or known mechanisms.\n"
        f"- If the gene is widely known to be essential or highly redundant, you may "
        f"  mention this qualitatively.\n"
        f"- If the gene is poorly characterized or evidence is unclear, explicitly "
        f"  state uncertainty.\n"
        f"- Do NOT imply certainty, diagnosis, lethality, or clinical recommendations.\n"
        f"- Return ONLY the sentence. No headers or explanations."
    )
    
    try:
        response = client.models.generate_content(
            model=model_name,
            contents=prompt
        )
        
        # Extract text from response
        if hasattr(response, 'text'):
            annotation = response.text.strip()
            # Ensure it's roughly the right length (allow some flexibility)
            if len(annotation) > 300:  # Too long, truncate
                annotation = annotation[:297] + "..."
            return annotation
        else:
            # Log that response doesn't have text attribute for debugging
            print(f"  WARNING: Response for {gene_symbol} has no 'text' attribute")
            return "NA"
    
    except Exception as e:
        # Fail gracefully - return "NA" on any error
        # Log first error for debugging (helps identify issues)
        if not hasattr(generate_consequence_annotation, '_error_logged'):
            print(f"  WARNING: Gemini API call failed for {gene_symbol}: {type(e).__name__}: {str(e)[:100]}")
            generate_consequence_annotation._error_logged = True
        return "NA"


def annotate_hit(
    hit: Dict[str, Any],
    client: Optional[genai.Client],
    model_name: str = "models/gemini-2.5-flash"
) -> Dict[str, Any]:
    """
    Annotate a single hit with allelic status and Gemini-generated consequence annotation.
    
    Args:
        hit: Dictionary containing hit information (must have gene_symbol, transcript_id, etc.)
        client: Gemini client (None if API unavailable)
        model_name: Gemini model to use
    
    Returns:
        Updated hit dictionary with 'allelic_status' and 'gemini_annotation' fields
    """
    # Add allelic status (always monoallelic assumed)
    hit['allelic_status'] = "monoallelic (assumed)"
    
    # Generate Gemini annotation if client is available
    if client is None:
        hit['gemini_annotation'] = "NA"
        return hit
    
    # Extract required fields (handle both old and new CSV formats)
    # Old format: gene_name, region_type, hit_position
    # New format: gene_symbol, transcript_type, match_start, match_end
    gene_symbol = hit.get('gene_symbol') or hit.get('gene_name', 'NA')
    transcript_id = hit.get('transcript_id', 'NA')
    transcript_type = hit.get('transcript_type') or hit.get('region_type', 'NA')
    
    # Handle edit_distance (should always be present)
    try:
        edit_distance = int(hit.get('edit_distance', 'NA'))
    except (ValueError, TypeError):
        edit_distance = 'NA'
    
    # Skip if essential information is missing
    if gene_symbol == 'NA' or transcript_id == 'NA':
        hit['gemini_annotation'] = "NA"
        return hit
    
    # Generate annotation
    annotation = generate_consequence_annotation(
        client=client,
        gene_symbol=gene_symbol,
        transcript_id=transcript_id,
        transcript_type=transcript_type,
        edit_distance=edit_distance,
        model_name=model_name
    )
    
    hit['gemini_annotation'] = annotation
    return hit

