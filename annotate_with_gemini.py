#!/usr/bin/env python3
"""
CLI entry point for Gemini-based biological consequence annotation.

Reads an existing results CSV, adds allelic_status and gemini_annotation columns,
and writes an annotated CSV.

Usage:
    python3 annotate_with_gemini.py input_results.csv annotated_results.csv
"""

import sys
import csv
from pathlib import Path
from src.gemini_annotate import create_gemini_client, annotate_hit


def read_results_csv(filepath: str) -> list:
    """
    Read results CSV into a list of dictionaries.
    
    Args:
        filepath: Path to input CSV file
    
    Returns:
        List of dictionaries, one per row
    """
    results = []
    
    try:
        with open(filepath, 'r', newline='') as f:
            reader = csv.DictReader(f)
            results = list(reader)
    except FileNotFoundError:
        print(f"ERROR: Input file not found: {filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read CSV: {e}")
        sys.exit(1)
    
    if not results:
        print("WARNING: Input CSV is empty")
    
    return results


def write_annotated_csv(results: list, filepath: str):
    """
    Write annotated results to CSV.
    
    Args:
        results: List of dictionaries with annotation fields added
        filepath: Path to output CSV file
    """
    if not results:
        print("WARNING: No results to write")
        return
    
    # Get all fieldnames from first row, ensuring new fields are included
    fieldnames = list(results[0].keys())
    
    # Ensure allelic_status and gemini_annotation are at the end
    if 'allelic_status' in fieldnames:
        fieldnames.remove('allelic_status')
    if 'gemini_annotation' in fieldnames:
        fieldnames.remove('gemini_annotation')
    
    # Add new fields at the end
    fieldnames.extend(['allelic_status', 'gemini_annotation'])
    
    try:
        with open(filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
    except Exception as e:
        print(f"ERROR: Failed to write CSV: {e}")
        sys.exit(1)


def main():
    """
    Main CLI function.
    """
    # Parse command line arguments
    if len(sys.argv) != 3:
        print("Usage: python3 annotate_with_gemini.py input_results.csv annotated_results.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print("=" * 80)
    print("Gemini Biological Consequence Annotation")
    print("=" * 80)
    print()
    
    # Step 1: Read input CSV
    print(f"Step 1: Reading input CSV: {input_file}")
    try:
        results = read_results_csv(input_file)
        print(f"  Loaded {len(results)} hits")
    except Exception as e:
        print(f"  ERROR: {e}")
        sys.exit(1)
    
    # Step 2: Initialize Gemini client
    print(f"\nStep 2: Initializing Gemini client...")
    client = create_gemini_client()
    if client is None:
        print("  WARNING: GEMINI_API_KEY not set or invalid - annotations will be 'NA'")
        use_gemini = False
    else:
        print("  ✓ Gemini client initialized")
        use_gemini = True
    
    # Step 3: Annotate each hit
    print(f"\nStep 3: Annotating hits with Gemini...")
    print(f"  Processing {len(results)} hits...")
    
    annotated_count = 0
    na_count = 0
    
    for idx, hit in enumerate(results, 1):
        if idx % 10 == 0:
            print(f"  Progress: {idx}/{len(results)} hits processed...")
        
        # Annotate hit
        annotated_hit = annotate_hit(hit, client, model_name="models/gemini-2.5-flash")
        
        if annotated_hit['gemini_annotation'] != 'NA':
            annotated_count += 1
        else:
            na_count += 1
        
        results[idx - 1] = annotated_hit
    
    print(f"  Completed: {annotated_count} annotated, {na_count} marked as 'NA'")
    
    # Step 4: Write annotated CSV
    print(f"\nStep 4: Writing annotated CSV: {output_file}")
    try:
        write_annotated_csv(results, output_file)
        print(f"  ✓ Annotated results written successfully")
    except Exception as e:
        print(f"  ERROR: {e}")
        sys.exit(1)
    
    print()
    print("=" * 80)
    print("Annotation completed successfully!")
    print(f"Output saved to: {output_file}")
    print("=" * 80)


if __name__ == "__main__":
    main()

