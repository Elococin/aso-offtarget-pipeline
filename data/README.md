# Data Directory

## Included (for testing / review)

- `aso_sequences.txt` — default ASO input list used by `main.py`
- `mock_transcripts.fa` — small transcript FASTA fixture so the pipeline runs without downloads
- `md5sum.txt` — checksums for included fixtures

## External (download separately)

- `grch38_refseq_transcripts.fa` — full RefSeq transcriptome FASTA (large)
- `chr12.fa` — chromosome reference FASTA used by allele-specific scripts (large)
- Patient-specific FASTA/VCF inputs — do not add if they may contain sensitive data

The repo ignores large reference files by default (see `.gitignore`).

---

# NCBI Datasets

https://www.ncbi.nlm.nih.gov/datasets

This zip archive contains an NCBI Datasets Data Package.

NCBI Datasets Data Packages can include sequence, annotation and other data files, and metadata in one or more data report files.
Data report files are in JSON Lines format.

---
## FAQs
### Where is the data I requested?

Your data is in the subdirectory `ncbi_dataset/data/` contained within this zip archive.

### I still can't find my data, can you help?

We have identified a bug affecting Mac Safari users. When downloading data from the NCBI Datasets web interface, you may see only this README file after the download has completed (while other files appear to be missing).
As a workaround to prevent this issue from recurring, we recommend disabling automatic zip archive extraction in Safari until Apple releases a bug fix.
For more information, visit:
https://www.ncbi.nlm.nih.gov/datasets/docs/reference-docs/mac-zip-bug/

### How do I work with JSON Lines data reports?

Visit our JSON Lines data report documentation page:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/working-with-jsonl-data-reports/

### What is NCBI Datasets?

NCBI Datasets is a resource that lets you easily gather data from across NCBI databases. Find and download gene, transcript, protein and genome sequences, annotation and metadata.

### Where can I find NCBI Datasets documentation?

Visit the NCBI Datasets documentation pages:
https://www.ncbi.nlm.nih.gov/datasets/docs/

---

National Center for Biotechnology Information
National Library of Medicine
info@ncbi.nlm.nih.gov
