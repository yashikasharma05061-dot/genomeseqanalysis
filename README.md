# Genome Sequence Analysis

A small Python pipeline using BioPython that aligns a reference DNA sequence to a sample, finds point mutations, and annotates their effect on the protein.

## What it does

1. **Aligns** the reference and sample sequences (global pairwise alignment, `Bio.Align.PairwiseAligner`)
2. **Identifies** positions where the two sequences differ
3. **Annotates** each variant by translating the affected codon in both sequences and classifying it as:
   - **synonymous**: the DNA changes but the amino acid does not
   - **missense**: the amino acid changes
   - **nonsense**: the change creates a stop codon

## Example output

```
Alignment score: 34.0

Variants found:
Position 20: G>A | codon 7: CGC>CAC | R>H (missense)
Position 27: G>A | codon 9: GGG>GGA | G>G (synonymous)
```

## How to run

```
git clone https://github.com/yashikasharma05061-dot/genomeseqanalysis.git
cd genomeseqanalysis
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python src/pipeline.py
```

## Limitations

- Currently runs on short example sequences written directly in the script
- Handles substitutions only, not insertions or deletions
- Assumes the sequence is a coding region starting at the first base

## Planned improvements

- Read sequences from FASTA files
- Handle insertions and deletions
- Add tests and a command-line interface
