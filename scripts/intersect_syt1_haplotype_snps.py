#!/usr/bin/env python3
"""
Intersect haplotype-correct SNPs with curated SYT1 design list.

Inputs:
  - het_same_haplotype_snps.csv (phased SNPs on pathogenic haplotype)
  - syt1.csv (curated candidate SNP list)

Outputs:
  - syt1_category2_final.csv (SNPs present in both lists, labeled Category 2)

Logic:
  1) Load both CSVs
  2) Normalize chromosome formatting (e.g., "12" vs "chr12")
  3) Join by (chrom, pos)
  4) Keep only SNPs present in BOTH files
  5) Exclude indels (len(ref)!=1 or len(alt)!=1) and rows missing coordinates
  6) Label retained rows:
       category = "Category 2: intronic haplotype SNP"
       reason   = "Intronic SNP phased to pathogenic SYT1 allele"

Usage:
    python scripts/intersect_syt1_haplotype_snps.py \
        --hap het_same_haplotype_snps.csv \
        --syt syt1.csv \
        --output syt1_category2_final.csv
"""

import argparse
import csv
import sys
from typing import Dict, List, Optional, Tuple


def read_csv(path: str) -> List[Dict]:
    try:
        with open(path, "r", newline="") as f:
            return list(csv.DictReader(f))
    except FileNotFoundError:
        sys.exit(f"ERROR: file not found: {path}")


def get_field(row: Dict, candidates: List[str]) -> Optional[str]:
    for key in candidates:
        if key in row and row[key] not in (None, ""):
            return row[key]
    # try case-insensitive
    lower_map = {k.lower(): v for k, v in row.items()}
    for key in candidates:
        if key.lower() in lower_map and lower_map[key.lower()] not in (None, ""):
            return lower_map[key.lower()]
    return None


def normalize_chrom(chrom: str) -> str:
    if chrom is None:
        return ""
    c = chrom.strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def to_int(value: str) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def is_biallelic_snp(ref: str, alt: str) -> bool:
    return ref is not None and alt is not None and len(ref) == 1 and len(alt) == 1


def main():
    parser = argparse.ArgumentParser(description="Intersect haplotype SNPs with SYT1 design list.")
    parser.add_argument("--hap", default="het_same_haplotype_snps.csv", help="Haplotype SNPs CSV")
    parser.add_argument("--syt", default="syt1.csv", help="Curated SYT1 design CSV")
    parser.add_argument("--output", default="syt1_category2_final.csv", help="Output CSV")
    args = parser.parse_args()

    hap_rows = read_csv(args.hap)
    syt_rows = read_csv(args.syt)

    print(f"Loaded haplotype SNPs: {len(hap_rows)}")
    print(f"Loaded syt1 candidates: {len(syt_rows)}")

    # Candidate column names
    chrom_keys = ["chrom", "chr", "#chrom", "chromosome"]
    pos_keys = ["pos", "position", "start"]
    ref_keys = ["ref", "ref_allele"]
    alt_keys = ["alt", "alt_allele"]
    gt_keys = ["gt", "genotype"]
    variant_id_keys = ["variant_id", "id", "rsid", "name"]
    mut_gt_keys = ["mutation_gt"]
    mut_hap_keys = ["mutation_alt_haplotype"]

    # Index SYT1 by (chrom, pos)
    syt_index: Dict[Tuple[str, int], Dict] = {}
    for row in syt_rows:
        chrom = normalize_chrom(get_field(row, chrom_keys))
        pos = to_int(get_field(row, pos_keys))
        if not chrom or pos is None:
            continue
        ref = get_field(row, ref_keys)
        alt = get_field(row, alt_keys)
        if not is_biallelic_snp(ref, alt):
            continue
        syt_index[(chrom, pos)] = row

    out_rows: List[Dict] = []
    for row in hap_rows:
        chrom = normalize_chrom(get_field(row, chrom_keys))
        pos = to_int(get_field(row, pos_keys))
        if not chrom or pos is None:
            continue
        ref = get_field(row, ref_keys)
        alt = get_field(row, alt_keys)
        if not is_biallelic_snp(ref, alt):
            continue
        key = (chrom, pos)
        if key not in syt_index:
            continue

        gt = get_field(row, gt_keys) or ""
        variant_id = get_field(row, variant_id_keys) or ""
        mut_gt = get_field(row, mut_gt_keys) or ""
        mut_hap = get_field(row, mut_hap_keys) or ""

        out_rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gt": gt,
                "mutation_gt": mut_gt,
                "mutation_alt_haplotype": mut_hap,
                "variant_id": variant_id,
                "category": "Category 2: intronic haplotype SNP",
                "reason": "Intronic SNP phased to pathogenic SYT1 allele",
            }
        )

    # Write output
    fieldnames = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "gt",
        "mutation_gt",
        "mutation_alt_haplotype",
        "variant_id",
        "category",
        "reason",
    ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Retained SNPs: {len(out_rows)} -> {args.output}")
    print("Summary:")
    print(f"  syt1.csv rows: {len(syt_rows)}")
    print(f"  haplotype SNPs input: {len(hap_rows)}")
    print(f"  intersected SNPs: {len(out_rows)}")


if __name__ == "__main__":
    main()

