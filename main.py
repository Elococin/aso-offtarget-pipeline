#!/usr/bin/env python3
"""
ASO Off-Target Screening Pipeline - Main Entry Point

This script orchestrates the entire screening pipeline:
1. Reads ASO sequences
2. Reads transcript sequences
3. Scans for off-target hits (substitution-only distance ≤ 2 mismatches)
4. Outputs results to CSV

Usage:
    python main.py
"""

import os
import sys
from pathlib import Path

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent))

from src.io import read_aso_sequences, read_fasta, write_results_csv
from src.scan import scan_all_sequences


def main():
    """
    Main pipeline execution function.
    """
    # Define file paths
    script_dir = Path(__file__).parent
    data_dir = script_dir / "data"
    aso_file = data_dir / "aso_sequences.txt"
    # Use real RefSeq transcript FASTA (fallback to mock for testing)
    transcript_file = data_dir / "grch38_refseq_transcripts.fa"
    if not transcript_file.exists():
        # Fallback to mock data if RefSeq file not found
        transcript_file = data_dir / "mock_transcripts.fa"
        if transcript_file.exists():
            print("  WARNING: Using mock_transcripts.fa (RefSeq file not found)")
        else:
            print(f"  ERROR: Neither RefSeq nor mock transcript file found")
            sys.exit(1)
    output_file = script_dir / "results.csv"
    
    print("=" * 60)
    print("ASO Off-Target Screening Pipeline")
    print("=" * 60)
    print()
    
    # Step 1: Read ASO sequences
    print(f"Step 1: Reading ASO sequences from {aso_file}...")
    try:
        aso_list = read_aso_sequences(str(aso_file))
        print(f"  Loaded {len(aso_list)} ASO sequences")
    except Exception as e:
        print(f"  ERROR: Failed to read ASO sequences: {e}")
        sys.exit(1)
    
    # Step 2: Read transcript sequences
    print(f"Step 2: Reading transcript sequences from {transcript_file}...")
    try:
        transcript_list = read_fasta(str(transcript_file))
        print(f"  Loaded {len(transcript_list)} transcript sequences")
    except Exception as e:
        print(f"  ERROR: Failed to read transcript sequences: {e}")
        sys.exit(1)
    
    # Step 3: Scan for off-target hits
    print("Step 3: Scanning for off-target hits (substitution-only distance ≤ 2 mismatches)...")
    try:
        hits = scan_all_sequences(aso_list, transcript_list, max_edit_distance=2)
        print(f"  Found {len(hits)} potential off-target hits")
    except Exception as e:
        print(f"  ERROR: Failed during scanning: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Step 4: Write results
    print(f"Step 4: Writing results to {output_file}...")
    try:
        write_results_csv(hits, str(output_file))
        print(f"  Results written successfully")
    except Exception as e:
        print(f"  ERROR: Failed to write results: {e}")
        sys.exit(1)
    
    print()
    print("=" * 60)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {output_file}")
    print("=" * 60)
    
    # Print summary statistics
    if hits:
        print()
        print("Summary Statistics:")
        print(f"  Total hits: {len(hits)}")
        
        # Count by transcript type
        transcript_type_counts = {}
        for hit in hits:
            tt = hit.get("transcript_type", "NA")
            transcript_type_counts[tt] = transcript_type_counts.get(tt, 0) + 1
        
        print("  Hits by transcript type:")
        for transcript_type, count in sorted(transcript_type_counts.items()):
            print(f"    {transcript_type}: {count}")
        
        # Count by edit distance
        dist_counts = {}
        for hit in hits:
            dist = hit["edit_distance"]
            dist_counts[dist] = dist_counts.get(dist, 0) + 1
        
        print("  Hits by mismatch count (substitution-only):")
        for dist in sorted(dist_counts.keys()):
            print(f"    {dist} mismatch(es): {dist_counts[dist]}")
    else:
        print()
        print("No off-target hits found (substitution-only distance ≤ 2 mismatches).")


if __name__ == "__main__":
    main()

