#!/usr/bin/env python3
"""
Join mutation_check.csv with phased heterozygous SNPs and label SNPs that share
the same haplotype as the pathogenic mutation.

Inputs (defaults expect files in repo root):
  - mutation_check.csv   (from scripts/parse_vcf_targets.py)
  - het_phased_snps.csv  (from scripts/parse_vcf_targets.py)

Outputs:
  - het_same_haplotype_snps.csv  (only SNPs whose ALT is on the same haplotype as the mutation ALT)

Logic:
  - Determine which haplotype (left or right) carries the pathogenic ALT based on GT:
      * GT 0|1  -> ALT on right haplotype
      * GT 1|0  -> ALT on left haplotype
      * Otherwise -> haplotype unknown
  - For each phased SNP (GT with '|'), mark same_haplotype_as_mutation = True/False/NA
  - Output only SNPs where same_haplotype_as_mutation == True

Usage:
    python scripts/join_mutation_phase.py \
        --mutation mutation_check.csv \
        --phased het_phased_snps.csv \
        --output het_same_haplotype_snps.csv
"""

import argparse
import csv
import sys
from typing import Tuple, List, Dict, Optional


def load_single_mutation(path: str) -> Optional[Dict]:
    try:
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            if not rows:
                return None
            return rows[0]
    except FileNotFoundError:
        sys.exit(f"ERROR: mutation file not found: {path}")


def load_phased_snps(path: str) -> List[Dict]:
    try:
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            return list(reader)
    except FileNotFoundError:
        sys.exit(f"ERROR: phased SNP file not found: {path}")


def alt_haplotype(gt: str) -> Optional[str]:
    """
    Return 'left' if ALT on left haplotype (1|0),
    'right' if ALT on right haplotype (0|1),
    None otherwise.
    """
    if gt == "1|0":
        return "left"
    if gt == "0|1":
        return "right"
    return None


def main():
    parser = argparse.ArgumentParser(description="Join mutation check with phased SNPs and infer shared haplotype.")
    parser.add_argument("--mutation", default="mutation_check.csv", help="Path to mutation_check.csv")
    parser.add_argument("--phased", default="het_phased_snps.csv", help="Path to het_phased_snps.csv")
    parser.add_argument("--output", default="het_same_haplotype_snps.csv", help="Output CSV (filtered to same haplotype)")
    args = parser.parse_args()

    mutation = load_single_mutation(args.mutation)
    if mutation is None:
        sys.exit("ERROR: mutation_check.csv is empty or missing required rows.")

    mut_gt = mutation.get("gt", "./.")
    mut_alt_hap = alt_haplotype(mut_gt)

    phased_snps = load_phased_snps(args.phased)

    out_rows = []
    for snp in phased_snps:
        gt = snp.get("gt", "./.")
        snp_alt_hap = alt_haplotype(gt)
        if mut_alt_hap is None or snp_alt_hap is None:
            same = "NA"
        else:
            same = str(snp_alt_hap == mut_alt_hap)
        snp["same_haplotype_as_mutation"] = same
        snp["mutation_gt"] = mut_gt
        snp["mutation_alt_haplotype"] = mut_alt_hap if mut_alt_hap else "NA"

        if same == "True":
            out_rows.append(snp)

    fieldnames = list(out_rows[0].keys()) if out_rows else [
        "chrom",
        "pos",
        "ref",
        "alt",
        "gt",
        "phase",
        "variant_id",
        "info",
        "same_haplotype_as_mutation",
        "mutation_gt",
        "mutation_alt_haplotype",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Mutation GT: {mut_gt} (ALT haplotype: {mut_alt_hap or 'unknown'})")
    print(f"Phased SNPs input: {len(phased_snps)}")
    print(f"SNPs on same haplotype: {len(out_rows)} -> {args.output}")


if __name__ == "__main__":
    main()

