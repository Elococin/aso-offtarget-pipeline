# ASO Off-Target Screening Pipeline

Fast, conservative screening for potential ASO off-target hits against the human transcriptome.

## Quick Start

```bash
# Clone and install
git clone https://github.com/violet-research-institute/aso-offtarget-pipeline
cd aso-offtarget-pipeline
pip install -r requirements.txt  # Optional dependencies only

# Run (works immediately with mock data)
python main.py

# For production speed (15x faster)
pypy3 main.py
```

## What It Does

Scans ASO sequences against RefSeq transcripts to find potential off-targets using:
- **Hamming distance** (substitution-only, no indels)
- **≤2 mismatch threshold**
- **Sliding window** across all transcripts

Results saved to `results.csv` with gene annotations.

## Input Files

**ASO sequences** (`data/aso_sequences.txt`):
```
ASO_001 ATCGATCGATCGATCGATCG
ASO_002 GCTAGCTAGCTAGCTAGCTA
```

**Transcripts**: Uses `data/grch38_refseq_transcripts.fa` if available, falls back to mock data.

## Performance

- **Mock data**: ~seconds (testing only)
- **Full RefSeq (CPython)**: ~5-12 min/ASO
- **Full RefSeq (PyPy)**: ~45 sec/ASO (recommended)

See `validation/VALIDATION_SUMMARY.md` for benchmarks.

## Validation

Internal validation docs for reviewers:
- `validation/ASO_OffTarget_Pipeline_Validation_Report.md` - Methodology overview
- `validation/results_syt1_5aso_test.csv` - Test results

## Project Structure

```
├── main.py                  # Pipeline entry point
├── src/                     # Core library (scan, edit_distance, io, annotate)
├── scripts/                 # Analysis utilities (VCF parsing, inventory builder)
├── data/                    # Inputs + reference (large files gitignored)
├── validation/              # Reviewer packet
└── docs/                    # Analysis reports
```

## Key Limitations

- **Approximation**: Hamming distance is a rough proxy for binding; doesn't model thermodynamics
- **Conservative**: 2-mismatch threshold may miss some real off-targets
- **RefSeq only**: Doesn't screen against whole genome or pseudogenes
- **No clinical use**: Internal research tool, requires expert review

## Optional Features

**Gemini annotation** (adds biological context to hits):
```bash
export GEMINI_API_KEY="your-key"
python annotate_with_gemini.py results.csv annotated_results.csv
```

**Analysis scripts** in `scripts/`:
- `parse_vcf_targets.py` - Extract phased SNPs from VCF
- `join_mutation_phase.py` - Find haplotype-linked SNPs
- `build_aso_inventory_final.py` - Consolidate ASO designs from multiple sources

See individual script docstrings for usage.

## Contact

Nicole Ye (@Elococin) - Violet Research Institute
