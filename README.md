# ASO Off-Target Screening Pipeline

## Disclaimer

**This is a conservative screening tool designed for rapid performance validation and methodological prototyping, NOT clinical decision-making.**

This pipeline is intended for:
- Technical validation of computational approaches
- Understanding pipeline scalability
- Educational purposes

**DO NOT use this tool for:**
- Clinical decision-making
- Regulatory submissions
- Patient safety assessments

## Overview

This pipeline performs conservative in silico screening for potential ASO (Antisense Oligonucleotide) off-target hits. It flags sequences with substitution-only distance ≤ 2 mismatches against functional transcribed regions (exons, introns, annotated lncRNAs).

## Reference Data (hg38 / RefSeq)

This repository does not include the full hg38 / RefSeq transcript FASTA due to size and redistribution constraints.

The pipeline is designed to run against real hg38 / RefSeq transcript data when the reference FASTA is available locally. If the hg38 file is not found, the pipeline automatically falls back to the small mock transcript dataset included in `data/` for testing and demonstration purposes.

## Key Assumptions and Limitations

1. **Off-target similarity metric**: Off-target binding is approximated using a substitution-only distance (Hamming distance), counting base mismatches between equal-length ASO and transcript windows. Insertions and deletions are ignored, as indels strongly disrupt ASO hybridization and are unlikely to produce meaningful binding.

2. **Mock data**: The repository includes small synthetic/public datasets for testing; real hg38 / RefSeq data is used when available locally.
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

## Gemini-Based Biological Annotation (Optional)

The pipeline includes an optional post-processing step that uses Google's Gemini API to add conservative biological consequence annotations to off-target hits.

**How it works:**
- Gemini acts as a cautious computational biology assistant
- For each hit, generates one conservative sentence (~20-25 words) describing potential biological consequences
- Assumes monoallelic disruption and uses cautious language ("may", "could", "potentially")
- Focuses on general biological role or known mechanisms
- Explicitly states uncertainty for poorly characterized genes
- This annotation step runs strictly downstream of hit detection and does not affect off-target calling

**Usage:**
```bash
export GEMINI_API_KEY="your-api-key-here"
python3 annotate_with_gemini.py results.csv annotated_results.csv
```

**Output:**
The annotated CSV adds two columns:
- `allelic_status`: Always "monoallelic (assumed)"
- `gemini_annotation`: Conservative biological consequence sentence, or "NA" if API unavailable

**Note:** This feature requires a Gemini API key and is optional. The core pipeline works without it.

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
