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

## Quick Start

```bash
# Clone and run (core pipeline requires no external dependencies)
git clone <repository-url>
cd aso_offtarget_pipeline
python main.py
```

Results are written to `results.csv`. The pipeline works immediately with included mock data.

### Optional Dependencies

For optional features, install:
```bash
pip install -r requirements.txt
```

This installs:
- `google-genai>=0.2.0` - For Gemini annotation (optional)
- `pandas>=1.3.0` and `openpyxl>=3.0.0` - For ASO inventory builder (optional)

**Requirements:** Python 3.7+ (standard library only for core pipeline)

## Reference Data

### Included Mock Data

The repository includes small mock data files in `data/`:
- `aso_sequences.txt` - Sample ASO sequences
- `mock_transcripts.fa` - Mock transcript dataset
- `md5sum.txt` - Data checksums

These files allow the pipeline to run immediately without downloads.

### Real Reference Data (Optional)

For production screening, download GRCh38 RefSeq transcripts:
1. Visit [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets)
2. Search for "Homo sapiens" genome assembly GCF_000001405.40 (GRCh38.p14)
3. Download transcript FASTA
4. Place at: `data/grch38_refseq_transcripts.fa`

The pipeline automatically uses `data/grch38_refseq_transcripts.fa` if found, otherwise falls back to mock data.

**Mock data is sufficient for:** Testing, workflow understanding, development, debugging.

**Real RefSeq data is recommended for:** Production screening, comprehensive off-target analysis, research applications.

## Usage

### Basic Pipeline

```bash
python main.py
```

The pipeline:
1. Reads ASO sequences from `data/aso_sequences.txt`
2. Reads transcript sequences from `data/grch38_refseq_transcripts.fa` (if available) or `data/mock_transcripts.fa` (fallback)
3. Scans for off-target hits (substitution-only distance ≤ 2 mismatches)
4. Outputs results to `results.csv`

### Input Format

**ASO Sequences** (`data/aso_sequences.txt`):
```
ASO_001 ATCGATCGATCGATCGATCG
ASO_002 GCTAGCTAGCTAGCTAGCTA
```
Format: `ASO_ID SEQUENCE` (space-separated, one per line)

**Transcript Sequences:**

Mock format (`data/mock_transcripts.fa`):
```
>GENE1|exon
ATCGATCGATCGATCGATCGATCGATCG
```

RefSeq format (`data/grch38_refseq_transcripts.fa`):
```
>NM_000014.6 Homo sapiens alpha-2-macroglobulin (A2M), transcript variant 1, mRNA
ATCGATCGATCGATCGATCGATCGATCG...
```

The pipeline automatically parses both formats.

### Output Format

`results.csv` columns:
- `aso_id` - ASO identifier
- `aso_sequence` - ASO nucleotide sequence (5' to 3')
- `transcript_id` - RefSeq transcript ID (e.g., `NM_000014.6`) or full header for mock data
- `gene_symbol` - Gene symbol (e.g., `A2M`) or `NA` if not available
- `transcript_type` - Transcript type (`mRNA`, `non-coding RNA`, `lncRNA`) or `NA`
- `match_start` - Start position of match in transcript (0-indexed)
- `match_end` - End position of match in transcript (0-indexed, exclusive)
- `matched_sequence` - Actual sequence from transcript that matched the ASO
- `edit_distance` - Mismatch count (substitution-only distance, 0-2)

Example:
```csv
aso_id,aso_sequence,transcript_id,gene_symbol,transcript_type,match_start,match_end,matched_sequence,edit_distance
ASO_001,ATCGATCGATCGATCGATCG,NM_000014.6,A2M,mRNA,100,120,ATCGATCGATCGATCGATCG,0
```

## Optional Features

### Gemini-Based Biological Annotation

Optional post-processing step using Google's Gemini API to add conservative biological consequence annotations to off-target hits.

**How it works:**
- Generates one conservative sentence (~20-25 words) per hit describing potential biological consequences
- Assumes monoallelic disruption and uses cautious language ("may", "could", "potentially")
- Focuses on general biological role or known mechanisms
- Explicitly states uncertainty for poorly characterized genes
- Runs strictly downstream of hit detection and does not affect off-target calling

**Usage:**
```bash
export GEMINI_API_KEY="your-api-key-here"
python annotate_with_gemini.py results.csv annotated_results.csv
```

**Output:** Adds two columns:
- `allelic_status` - Always `"monoallelic (assumed)"`
- `gemini_annotation` - Conservative biological consequence sentence, or `"NA"` if API unavailable

**Note:** Requires a Gemini API key. The core pipeline works without it.

### Allele-Specific Analysis Scripts

Additional scripts for allele-specific ASO design analysis:

#### 1. VCF Target Extraction

Parse a phased VCF to extract pathogenic mutations and heterozygous SNPs:

```bash
python scripts/parse_vcf_targets.py \
  --vcf /path/to/phased.vcf \
  --chrom 12 \
  --mutation-pos 79448958
```

**Arguments:**
- `--vcf` - Path to phased VCF file (provide your own)
- `--chrom` - Chromosome number (e.g., 12) - replace with your chromosome
- `--mutation-pos` - Genomic position (GRCh38 coordinates, e.g., 79448958) - replace with your mutation position

**Outputs:**
- `mutation_check.csv` - Confirmed mutation details
  - Columns: `chrom`, `pos`, `ref`, `alt`, `gt` (genotype), `variant_id`
- `het_phased_snps.csv` - Phased heterozygous SNPs (GT contains `|`)
  - Columns: `chrom`, `pos`, `ref`, `alt`, `gt`, `phase`, `variant_id`, `info`
- `het_unphased_snps.csv` - Unphased heterozygous SNPs (GT is `0/1`)
  - Columns: `chrom`, `pos`, `ref`, `alt`, `gt`, `variant_id`, `info`

**Note:** 
- Provide your own VCF file and mutation coordinates. Output files are created in the current directory.
- Only processes biallelic SNVs (indels excluded).
- Phased SNPs have `|` in GT field (e.g., `0|1`, `1|0`); unphased have `/` (e.g., `0/1`).

#### 2. Haplotype Phase Joining

Join mutation info with phased SNPs to identify SNPs on the same haplotype:

```bash
python scripts/join_mutation_phase.py \
  --mutation mutation_check.csv \
  --phased het_phased_snps.csv \
  --output output_same_haplotype_snps.csv
```

**Input files:**
- `mutation_check.csv` - Output from step 1 (mutation details)
- `het_phased_snps.csv` - Output from step 1 (phased SNPs)

**Required columns:**
- Both files must contain `gt` (genotype) column with phased genotypes (e.g., `0|1`, `1|0`)
- Script determines which haplotype carries the pathogenic ALT based on GT:
  - `GT 0|1` → ALT on right haplotype
  - `GT 1|0` → ALT on left haplotype

**Output:** Filtered CSV containing only SNPs on the same haplotype as the mutation.

**Note:** Provide your own input files from step 1.

#### 3. Candidate Intersection

Intersect haplotype-correct SNPs with a curated candidate list:

```bash
python scripts/intersect_syt1_haplotype_snps.py \
  --hap het_same_haplotype_snps.csv \
  --syt your_candidate_list.csv \
  --output output_category2_final.csv
```

**Input files:**
- `het_same_haplotype_snps.csv` - Output from step 2 (haplotype-correct SNPs)
- `your_candidate_list.csv` - Your curated candidate SNP list

**Required columns:**
- Both files must contain chromosome and position columns (script accepts: `chrom`/`chr`/`chromosome`, `pos`/`position`/`start`)
- Both files must contain `ref` and `alt` columns (or `ref_allele`/`alt_allele`)
- Script automatically normalizes chromosome formatting (e.g., "12" vs "chr12")
- Only biallelic SNPs are kept (indels excluded)

**Output:** Intersection results with category labels and source tracking.

**Note:** Input files are not included in the repository (patient-specific data). You must provide your own candidate lists.

#### 4. Mutation-Based Off-Target Analysis

Extract genomic window around a mutation and run off-target analysis:

```bash
python scripts/run_mutation_offtarget.py \
  --mutation-check mutation_check.csv \
  --genome /path/to/chromosome.fa \
  --transcripts data/grch38_refseq_transcripts.fa \
  --window-size 25
```

**Input files:**
- `mutation_check.csv` - Output from step 1 (mutation details)
- Genome FASTA - Path to chromosome-specific FASTA file (e.g., `chr12.fa`)
- Transcripts - Path to RefSeq transcript FASTA (or use `data/grch38_refseq_transcripts.fa` if available)

**Outputs:**
- `mutation_sequences.csv` - Mutant and WT sequences
- `mutation_offtarget_hits.csv` - Off-target hits for mutant sequence

**Note:** You must provide your own mutation check file and genome FASTA file. The genome file is not included in the repository.

#### 5. ASO Inventory Builder

Build a canonical ASO inventory from multiple sources (CSV and Excel files):

```bash
# Requires: pip install pandas openpyxl
python scripts/build_aso_inventory_final.py
```

**Input files:**
The script looks for ASO sequence files in the repository root directory:
- CSV files: Currently looks for `syt1.csv` (rename your file or modify line 210 in the script)
- Excel files: Any `.xlsx` files matching `Syt1*.xlsx` pattern (modify line 207 in the script to match your pattern)

**File format requirements:**
- CSV files: Should contain columns with ASO sequences in IDT/MOE notation (e.g., `/52MOErT/*/i2MOErA/...`)
  - Script auto-detects sequence columns (looks for `idt`, `sequence`, `seq`, etc.)
- Excel files: Should contain ASO sequences in any column (script auto-detects sequence columns)
- Sequences are automatically extracted from IDT/MOE format by extracting all A, T, C, G characters
- Handles both IDT format (`/52MOErT/...`) and LNA gapmer format (`*T*T*G*T*` or `+A*+T*+A`)

**Required modifications (if needed):**
- To change CSV filename: Edit line 210 in `scripts/build_aso_inventory_final.py`
- To change Excel pattern: Edit line 207 in `scripts/build_aso_inventory_final.py` (change `"Syt1*.xlsx"` to your pattern)

**Output:**
- `canonical_syt1_aso_table.csv` - Deduplicated ASO table with source tracking
  - Columns: `ASO_ID`, `Sequence_5to3`, `Length_nt`, `Source_File`, `Original_Row_ID_or_Name`
  - ASO IDs assigned sequentially (SYT1_ASO_001, SYT1_ASO_002, etc.)

**Dependencies:** Requires `pandas` and `openpyxl`. Install with: `pip install pandas openpyxl`

**Note:** 
- Input files are not included in the repository (patient-specific data). You must provide your own ASO sequence files.
- Script will skip missing files with a warning and process any files that are found.

## Algorithm

For each ASO:
1. Slide a window of exactly ASO length across each transcript sequence
2. Compute substitution-only distance (mismatch count) between ASO and each window
3. Flag hits where mismatch count ≤ 2
4. Annotate with gene symbol and transcript type from FASTA header

**Substitution-only distance (Hamming distance)** counts only base mismatches between equal-length sequences. Insertions and deletions are ignored, as indels strongly disrupt ASO hybridization.

## Key Assumptions and Limitations

1. **Off-target similarity metric**: Off-target binding is approximated using a substitution-only distance (Hamming distance), counting base mismatches between equal-length ASO and transcript windows. Insertions and deletions are ignored, as indels strongly disrupt ASO hybridization and are unlikely to produce meaningful binding.

2. **Mock data**: The repository includes small synthetic/public datasets for testing; real hg38 / RefSeq data is used when available locally.

3. **Simple Algorithm**: Basic substitution-only distance - no sophisticated ASO cutting prediction

4. **Conservative Threshold**: Substitution-only distance ≤ 2 mismatches (may miss some hits, but reduces false positives)

5. **No Clinical Interpretation**: Outputs raw matches - requires human operator review

6. **Scalability**: Designed for up to 40,000 ASOs, but tested with small mock datasets

## Project Structure

```
aso_offtarget_pipeline/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── main.py                      # Main pipeline entry point
├── annotate_with_gemini.py     # Optional Gemini annotation CLI
├── src/
│   ├── __init__.py              # Package initialization
│   ├── io.py                    # Input/output utilities
│   ├── edit_distance.py         # Hamming distance computation
│   ├── scan.py                  # Sliding window and hit detection
│   ├── annotate.py              # FASTA header parsing
│   └── gemini_annotate.py       # Gemini API integration
├── scripts/
│   ├── parse_vcf_targets.py     # VCF parsing for allele-specific targets
│   ├── join_mutation_phase.py  # Haplotype phase joining
│   ├── intersect_syt1_haplotype_snps.py  # Candidate intersection
│   ├── run_mutation_offtarget.py  # Mutation-based off-target analysis
│   └── build_aso_inventory_final.py  # Canonical ASO inventory builder
├── data/
│   ├── aso_sequences.txt        # Input ASO sequences (included)
│   ├── mock_transcripts.fa      # Mock transcript sequences (included)
│   ├── md5sum.txt               # Data checksums
│   ├── README.md                # Data directory documentation
│   └── grch38_refseq_transcripts.fa  # Real RefSeq data (user downloads)
└── report/
    └── syt1_allele_specific_offtarget_report.md  # Example analysis report (user-generated)
```

## Performance Notes

- Designed for up to 40,000 ASOs
- ASO length: ~20 bp
- Search space: functional transcribed regions only
- No GPU required
- No external bioinformatics tools required
- Runtime: ~1-2 minutes per ASO against full RefSeq (estimated)

## Troubleshooting

### "RefSeq file not found" warning

**This is expected and normal!** The pipeline automatically uses `data/mock_transcripts.fa` as a fallback. This is fine for testing and allows the pipeline to work immediately without downloads.

**To use real data (optional):** Download GRCh38 RefSeq transcripts and place at `data/grch38_refseq_transcripts.fa`. See "Reference Data" section above for download instructions.

### "No hits found"

**Possible causes:**
- ASO sequences are too short or contain invalid characters
- Edit distance threshold (≤2) is too strict
- Transcript sequences are shorter than ASO length

**Check:**
- Verify ASO sequences in `data/aso_sequences.txt` are valid (A, T, C, G only)
- Ensure ASO length is ~20 bp
- Review transcript sequences in your FASTA file

### Gemini annotation returns "NA"

**Possible causes:**
- `GEMINI_API_KEY` environment variable not set
- API key is invalid or expired
- Network connectivity issues
- API rate limits exceeded

**Check:**
```bash
echo $GEMINI_API_KEY  # Should show your API key
python -c "from google import genai; print('OK')"  # Test import
```

### Import errors

**For optional features:**
- **Gemini annotation:** `pip install google-genai`
- **ASO inventory builder:** `pip install pandas openpyxl`
- **Core pipeline:** No external dependencies needed - check Python version instead

**Python version:**
This pipeline requires Python 3.7+. Check your version:
```bash
python --version  # Should show 3.7 or higher
python3 --version  # Try python3 if python doesn't work
```

**Standard library modules:**
If you see `ModuleNotFoundError` for standard library modules, this indicates a Python installation issue. The core pipeline uses only standard library modules (`csv`, `sys`, `os`, `pathlib`, `re`, `collections`). These should be available in any Python 3.7+ installation.

## Contributing

This is a research prototype. For questions or issues, please open an issue in the repository.

## License

MIT License - see LICENSE file for details.


## AuthorGitHub: 
Nicole Ye @Elococin
