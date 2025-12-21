# ASO Off-Target Screening Pipeline

## ⚠️ Important Disclaimer

**This is a mock / conservative screening tool designed for performance validation, NOT clinical decision-making.**

This pipeline is a rapid prototype (48-hour development) intended for:
- Technical validation of computational approaches
- Understanding pipeline scalability
- Educational purposes

**DO NOT use this tool for:**
- Clinical decision-making
- Regulatory submissions
- Patient safety assessments

## Overview

This pipeline performs conservative in silico screening for potential ASO (Antisense Oligonucleotide) off-target hits. It flags sequences with substitution-only distance ≤ 2 mismatches against functional transcribed regions (exons, introns, annotated lncRNAs).

## Key Assumptions and Limitations

1. **Off-target similarity metric**: Off-target binding is approximated using a substitution-only distance (Hamming distance), counting base mismatches between equal-length ASO and transcript windows. Insertions and deletions are ignored, as indels strongly disrupt ASO hybridization and are unlikely to produce meaningful binding.

2. **Mock Data**: Uses synthetic/public genomic data - NOT biologically validated
3. **Simple Algorithm**: Basic substitution-only distance - no sophisticated ASO cutting prediction
4. **Conservative Threshold**: Substitution-only distance ≤ 2 mismatches (may miss some hits, but reduces false positives)
5. **No Clinical Interpretation**: Outputs raw matches - requires human operator review
6. **Scalability**: Designed for up to 40,000 ASOs, but tested with small mock datasets

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python main.py
```

The pipeline will:
1. Read ASO sequences from `data/aso_sequences.txt`
2. Read transcript sequences from `data/mock_transcripts.fa`
3. Scan for off-target hits (substitution-only distance ≤ 2 mismatches)
4. Output results to `results.csv`

## Input Format

### ASO Sequences (`data/aso_sequences.txt`)
```
ASO_001 ATCGATCGATCGATCGATCG
ASO_002 GCTAGCTAGCTAGCTAGCTA
```

Format: `ASO_ID SEQUENCE` (space-separated, one per line)

### Transcript Sequences (`data/mock_transcripts.fa`)
```
>GENE1|exon
ATCGATCGATCGATCGATCGATCGATCG
>GENE2|intron
GCTAGCTAGCTAGCTAGCTAGCTAGCTA
>GENE3|lncRNA
TTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

Format: FASTA with headers encoding `GENE_NAME|region_type`

## Output Format

`results.csv` contains:
- `aso_id`: ASO identifier
- `gene_name`: Gene name from transcript header
- `region_type`: Type of region (exon/intron/lncRNA)
- `transcript_id`: Full transcript identifier (header)
- `hit_position`: Position in transcript where hit starts (0-indexed)
- `edit_distance`: Computed mismatch count (substitution-only distance, 0-2)

## Algorithm

For each ASO:
1. Slide a window of exactly ASO length across each transcript sequence
2. Compute substitution-only distance (mismatch count) between ASO and each window
3. Flag hits where mismatch count ≤ 2
4. Annotate with gene name and region type from FASTA header

Substitution-only distance (Hamming distance) counts only base mismatches between equal-length sequences. Insertions and deletions are ignored, as indels strongly disrupt ASO hybridization.

## Project Structure

```
aso_offtarget_pipeline/
├── README.md              # This file
├── requirements.txt       # Python dependencies
├── main.py               # Main pipeline entry point
├── src/
│   ├── __init__.py       # Package initialization
│   ├── io.py             # Input/output utilities
│   ├── edit_distance.py  # Edit distance computation
│   ├── scan.py           # Sliding window and hit detection
│   └── annotate.py       # FASTA header parsing
├── data/
│   ├── aso_sequences.txt # Input ASO sequences
│   └── mock_transcripts.fa # Mock transcript sequences
└── results.csv           # Output results (generated)
```

## Performance Notes

- Designed for up to 40,000 ASOs
- ASO length: ~20 bp
- Search space: functional transcribed regions only
- No GPU required
- No external bioinformatics tools required

## Next Steps for Production

If this were to be production-ready, consider:
1. Biologically validated reference genomes
2. More sophisticated ASO binding prediction models
3. Clinical annotation databases
4. Performance optimization (k-mer indexing, parallelization)
5. Regulatory validation

