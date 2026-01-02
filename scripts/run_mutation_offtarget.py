#!/usr/bin/env python3
"""
Mutation-Based Allele-Specific Off-Target Analysis

Extracts a 51bp window around the pathogenic SYT1 mutation (chr12:79448958),
constructs mutant and WT sequences, and runs the mutant sequence through the
off-target scanning pipeline.

Usage:
    python3 scripts/run_mutation_offtarget.py \
        --mutation-check mutation_check.csv \
        --genome /path/to/hg38.fa \
        --transcripts /path/to/grch38_refseq_transcripts.fa \
        --chrom 12 \
        --pos 79448958 \
        --window 25

Outputs:
    mutation_sequences.csv        - Mutant and WT sequences
    mutation_offtarget_hits.csv  - Off-target hits for mutant sequence
"""

import argparse
import csv
import sys
from pathlib import Path

# Add parent directory to path to import pipeline modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.io import read_fasta, write_results_csv
from src.scan import scan_all_sequences


def read_mutation_check(filepath: str) -> dict:
    """
    Read mutation_check.csv and extract mutation information.
    
    Returns:
        Dictionary with keys: chrom, pos, ref, alt, gt, phased
    """
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        if not rows:
            raise ValueError(f"No rows found in {filepath}")
        return rows[0]  # Return first (and should be only) row


def extract_reference_sequence(genome_fasta: str, chrom: str, start: int, end: int) -> str:
    """
    Extract a reference sequence from a genome FASTA file.
    
    Args:
        genome_fasta: Path to genome FASTA (can be single chromosome or full genome)
        chrom: Chromosome name (e.g., "chr12" or "12")
        start: 0-based start position
        end: 0-based end position (exclusive)
    
    Returns:
        Reference sequence (uppercase)
    
    Note:
        Handles both single-chromosome FASTA files and full genome FASTA files.
        For full genome, expects headers like ">chr12" or ">12".
    """
    # Try to find matching chromosome
    target_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    alt_chrom = chrom if not chrom.startswith("chr") else chrom[3:]  # Remove "chr" prefix
    
    sequences = read_fasta(genome_fasta)
    
    # Look for matching chromosome
    chrom_sequence = None
    for header, seq in sequences:
        header_clean = header.split()[0].lower()  # Take first token, lowercase
        if header_clean == target_chrom.lower() or header_clean == alt_chrom.lower():
            chrom_sequence = seq
            break
    
    if chrom_sequence is None:
        available = [h.split()[0] for h, _ in sequences[:5]]
        raise ValueError(
            f"Chromosome {chrom} not found in {genome_fasta}. "
            f"Available headers (first 5): {available}"
        )
    
    # Extract window (0-based coordinates)
    if start < 0 or end > len(chrom_sequence):
        raise ValueError(
            f"Window [{start}:{end}] out of bounds for chromosome {chrom} "
            f"(length: {len(chrom_sequence)})"
        )
    
    return chrom_sequence[start:end].upper()


def construct_mutation_sequences(ref_seq: str, center_idx: int, ref_allele: str, alt_allele: str) -> tuple:
    """
    Construct mutant and WT sequences with ALT/REF at center position.
    
    Args:
        ref_seq: Reference sequence (51bp window)
        center_idx: Index of center position (25 for 51bp window)
        ref_allele: Reference allele (should match ref_seq[center_idx])
        alt_allele: Alternate allele
    
    Returns:
        (mutant_sequence, wt_sequence) tuple
    """
    # Verify center position matches REF
    if ref_seq[center_idx].upper() != ref_allele.upper():
        raise ValueError(
            f"Reference sequence mismatch at center: "
            f"expected {ref_allele}, got {ref_seq[center_idx]}"
        )
    
    # Construct sequences
    mutant_seq = ref_seq[:center_idx] + alt_allele.upper() + ref_seq[center_idx + 1:]
    wt_seq = ref_seq  # WT is just the reference sequence
    
    return mutant_seq, wt_seq


def save_sequences(mutant_seq: str, wt_seq: str, chrom: str, pos: int, output_file: str):
    """
    Save mutant and WT sequences to CSV.
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'target_name', 'chrom', 'pos', 'allele_type', 'sequence'
        ])
        writer.writeheader()
        writer.writerow({
            'target_name': 'SYT1_pathogenic_mutation',
            'chrom': chrom,
            'pos': pos,
            'allele_type': 'mutant',
            'sequence': mutant_seq
        })
        writer.writerow({
            'target_name': 'SYT1_pathogenic_mutation',
            'chrom': chrom,
            'pos': pos,
            'allele_type': 'wt',
            'sequence': wt_seq
        })


def main():
    parser = argparse.ArgumentParser(
        description='Run off-target analysis for mutation-based allele-specific target'
    )
    parser.add_argument(
        '--mutation-check',
        default='mutation_check.csv',
        help='Path to mutation_check.csv (default: mutation_check.csv)'
    )
    parser.add_argument(
        '--genome',
        required=True,
        help='Path to hg38 reference genome FASTA (can be single chromosome or full genome)'
    )
    parser.add_argument(
        '--transcripts',
        default='data/grch38_refseq_transcripts.fa',
        help='Path to RefSeq transcript FASTA (default: data/grch38_refseq_transcripts.fa)'
    )
    parser.add_argument(
        '--chrom',
        default='12',
        help='Chromosome (default: 12)'
    )
    parser.add_argument(
        '--pos',
        type=int,
        default=79448958,
        help='Mutation position (1-based, default: 79448958)'
    )
    parser.add_argument(
        '--window',
        type=int,
        default=25,
        help='Window size on each side of mutation (default: 25, total length = 51bp)'
    )
    parser.add_argument(
        '--output-sequences',
        default='mutation_sequences.csv',
        help='Output file for sequences (default: mutation_sequences.csv)'
    )
    parser.add_argument(
        '--output-hits',
        default='mutation_offtarget_hits.csv',
        help='Output file for off-target hits (default: mutation_offtarget_hits.csv)'
    )
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("Mutation-Based Allele-Specific Off-Target Analysis")
    print("=" * 80)
    print()
    
    # Step 1: Read mutation information
    print(f"Step 1: Reading mutation information from {args.mutation_check}...")
    try:
        mutation = read_mutation_check(args.mutation_check)
        chrom = mutation['chrom'].replace('chr', '')  # Normalize
        pos = int(mutation['pos'])
        ref_allele = mutation['ref']
        alt_allele = mutation['alt']
        print(f"  Mutation: {mutation['chrom']}:{pos} {ref_allele}>{alt_allele} (GT: {mutation['gt']})")
    except Exception as e:
        print(f"  ERROR: Failed to read mutation: {e}")
        sys.exit(1)
    
    # Override with CLI args if provided
    if args.chrom:
        chrom = args.chrom
    if args.pos:
        pos = args.pos
    
    # Step 2: Extract reference sequence window
    print(f"\nStep 2: Extracting reference sequence window...")
    print(f"  Chromosome: {chrom}, Position: {pos}, Window: ±{args.window}bp")
    
    # Convert 1-based position to 0-based for extraction
    center_0based = pos - 1
    start_0based = center_0based - args.window
    end_0based = center_0based + args.window + 1  # +1 for inclusive end
    center_idx = args.window  # Index within the extracted window
    
    try:
        ref_sequence = extract_reference_sequence(
            args.genome,
            chrom,
            start_0based,
            end_0based
        )
        print(f"  Extracted {len(ref_sequence)}bp reference sequence")
    except Exception as e:
        print(f"  ERROR: Failed to extract reference sequence: {e}")
        print(f"  Make sure --genome points to a valid hg38 FASTA file")
        sys.exit(1)
    
    # Step 3: Construct mutant and WT sequences
    print(f"\nStep 3: Constructing mutant and WT sequences...")
    try:
        mutant_seq, wt_seq = construct_mutation_sequences(
            ref_sequence,
            center_idx,
            ref_allele,
            alt_allele
        )
        print(f"  Mutant sequence: {mutant_seq}")
        print(f"  WT sequence:     {wt_seq}")
    except Exception as e:
        print(f"  ERROR: Failed to construct sequences: {e}")
        sys.exit(1)
    
    # Step 4: Save sequences
    print(f"\nStep 4: Saving sequences to {args.output_sequences}...")
    try:
        save_sequences(mutant_seq, wt_seq, chrom, pos, args.output_sequences)
        print(f"  Sequences saved successfully")
    except Exception as e:
        print(f"  ERROR: Failed to save sequences: {e}")
        sys.exit(1)
    
    # Step 5: Run off-target scan on mutant sequence only
    print(f"\nStep 5: Running off-target scan on mutant sequence...")
    print(f"  Transcript file: {args.transcripts}")
    
    try:
        # Load transcripts
        transcript_list = read_fasta(args.transcripts)
        print(f"  Loaded {len(transcript_list)} transcript sequences")
        
        # Create ASO list (only mutant sequence)
        aso_list = [("SYT1_mutant", mutant_seq)]
        
        # Run scan
        hits = scan_all_sequences(aso_list, transcript_list, max_edit_distance=2)
        print(f"  Found {len(hits)} potential off-target hits")
        
    except Exception as e:
        print(f"  ERROR: Failed during off-target scanning: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 6: Write results
    print(f"\nStep 6: Writing results to {args.output_hits}...")
    try:
        write_results_csv(hits, args.output_hits)
        print(f"  Results written successfully")
    except Exception as e:
        print(f"  ERROR: Failed to write results: {e}")
        sys.exit(1)
    
    # Step 7: Print summary
    print()
    print("=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"Total off-target hits: {len(hits)}")
    
    if hits:
        # Count hits in coding exons
        exon_hits = sum(1 for h in hits if h.get('transcript_type', '').lower() == 'mrna')
        print(f"Hits in coding exons (mRNA): {exon_hits}")
        
        # Top 5 genes
        gene_counts = {}
        for hit in hits:
            gene = hit.get('gene_symbol', 'NA')
            if gene != 'NA':
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        if gene_counts:
            print(f"\nTop 5 genes hit:")
            for gene, count in sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {gene}: {count} hit(s)")
    else:
        print("No off-target hits found.")
    
    print()
    print("=" * 80)
    print("Analysis complete!")
    print(f"Sequences: {args.output_sequences}")
    print(f"Hits: {args.output_hits}")
    print("=" * 80)


if __name__ == "__main__":
    main()

