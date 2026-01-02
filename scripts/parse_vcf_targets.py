#!/usr/bin/env python3
"""
VCF Target Extraction for Allele-Specific SYT1 Project

This script parses a phased VCF and produces three CSV outputs:
1) mutation_check.csv      - Confirms the pathogenic mutation at chr12:79448958 (GRCh38)
2) het_phased_snps.csv     - Phased heterozygous SNPs (0|1 or 1|0), biallelic SNVs only
3) het_unphased_snps.csv   - Unphased heterozygous SNPs (0/1), biallelic SNVs only

Defaults (can be overridden via CLI):
  VCF input:  /mnt/data/GB-01.joint.GRCh38.small_variants.phased.vcf.gz_chr12_78859774-79457008.vcf
  Outputs:    mutation_check.csv, het_phased_snps.csv, het_unphased_snps.csv (written to repo root)

Usage:
    python scripts/parse_vcf_targets.py \
        --vcf /path/to/file.vcf.gz \
        --mutation-pos 79448958 \
        --chrom 12

Notes:
- Pure Python; no heavy dependencies.
- Handles missing fields gracefully.
- Does not log full VCF contents.
"""

import argparse
import csv
import gzip
import sys
from typing import Dict, List, Tuple, Optional


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(sample_field: str, fmt: str) -> str:
    """
    Extract GT from a sample field given the FORMAT string.
    Returns './.' if not found.
    """
    fmt_keys = fmt.split(":")
    sample_values = sample_field.split(":")
    mapping = dict(zip(fmt_keys, sample_values))
    return mapping.get("GT", "./.")


def normalize_chrom(chrom: str) -> str:
    return chrom.replace("chr", "").replace("CHR", "")


def is_biallelic_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and "," not in alt


def is_het(gt: str) -> bool:
    alleles = gt.replace("|", "/").split("/")
    return set(alleles) == {"0", "1"}


def is_phased(gt: str) -> bool:
    return "|" in gt


def load_vcf_records(vcf_path: str) -> Tuple[List[str], List[Dict]]:
    """
    Load VCF records into a list of dicts with basic fields.
    """
    records = []
    samples: List[str] = []
    try:
        with open_maybe_gzip(vcf_path) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    parts = line.strip().split("\t")
                    samples = parts[9:]
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
                fmt = parts[8] if len(parts) > 8 else ""
                sample_fields = parts[9:] if len(parts) > 9 else []
                sample = sample_fields[0] if sample_fields else ""
                gt = parse_gt(sample, fmt) if fmt and sample else "./."

                records.append(
                    {
                        "chrom": chrom,
                        "pos": int(pos),
                        "id": vid if vid != "." else "",
                        "ref": ref,
                        "alt": alt,
                        "info": info,
                        "format": fmt,
                        "sample_raw": sample,
                        "gt": gt,
                    }
                )
    except FileNotFoundError:
        sys.exit(f"ERROR: VCF not found: {vcf_path}")
    return samples, records


def write_csv(path: str, rows: List[Dict], fieldnames: List[str]):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Parse phased VCF for allele-specific SYT1 targets")
    parser.add_argument(
        "--vcf",
        default="/mnt/data/GB-01.joint.GRCh38.small_variants.phased.vcf.gz_chr12_78859774-79457008.vcf",
        help="Input VCF (can be .vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--chrom",
        default="12",
        help="Target chromosome for pathogenic mutation (default: 12)",
    )
    parser.add_argument(
        "--mutation-pos",
        type=int,
        default=79448958,
        help="Position of pathogenic mutation (GRCh38) (default: 79448958)",
    )
    parser.add_argument(
        "--mutation-out",
        default="mutation_check.csv",
        help="Output CSV for pathogenic mutation check",
    )
    parser.add_argument(
        "--het-phased-out",
        default="het_phased_snps.csv",
        help="Output CSV for phased heterozygous SNPs",
    )
    parser.add_argument(
        "--het-unphased-out",
        default="het_unphased_snps.csv",
        help="Output CSV for unphased heterozygous SNPs",
    )

    args = parser.parse_args()

    samples, records = load_vcf_records(args.vcf)
    sample_name = samples[0] if samples else "SAMPLE"

    target_chrom = args.chrom.replace("chr", "").replace("CHR", "")

    mutation_rows = []
    het_phased_rows = []
    het_unphased_rows = []

    for rec in records:
        chrom_norm = normalize_chrom(rec["chrom"])
        pos = rec["pos"]
        ref = rec["ref"]
        alt = rec["alt"]
        gt = rec["gt"]
        vid = rec["id"]
        info = rec["info"]

        # Goal A: pathogenic mutation check
        if chrom_norm == target_chrom and pos == args.mutation_pos:
            mutation_rows.append(
                {
                    "chrom": rec["chrom"],
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gt": gt,
                    "phased": "yes" if is_phased(gt) else "no",
                }
            )

        # Goal B: heterozygous SNPs (biallelic SNVs only)
        if not is_biallelic_snp(ref, alt):
            continue
        if not is_het(gt):
            continue

        row = {
            "chrom": rec["chrom"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "gt": gt,
            "phase": "phased" if is_phased(gt) else "unphased",
            "variant_id": vid,
            "info": info,
        }

        if is_phased(gt):
            het_phased_rows.append(row)
        else:
            het_unphased_rows.append(row)

    # Write outputs
    write_csv(
        args.mutation_out,
        mutation_rows,
        ["chrom", "pos", "ref", "alt", "gt", "phased"],
    )
    write_csv(
        args.het_phased_out,
        het_phased_rows,
        ["chrom", "pos", "ref", "alt", "gt", "phase", "variant_id", "info"],
    )
    write_csv(
        args.het_unphased_out,
        het_unphased_rows,
        ["chrom", "pos", "ref", "alt", "gt", "phase", "variant_id", "info"],
    )

    print(f"Mutation check rows: {len(mutation_rows)} -> {args.mutation_out}")
    print(f"Phased het SNPs:    {len(het_phased_rows)} -> {args.het_phased_out}")
    print(f"Unphased het SNPs:  {len(het_unphased_rows)} -> {args.het_unphased_out}")


if __name__ == "__main__":
    main()

