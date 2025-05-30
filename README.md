# Query-Match

**Approximate‐substring matching for DNA sequences**

## Overview

This utility locates all positions where any length-*n* fragment of a “query” DNA string aligns to a length-*n* fragment of a “text” DNA string with at most *k* mismatches. An efficient *l*-mer filtration step (with *l* = ⌊*n*/(*k*+1)⌋) is used to rapidly discard impossible alignments before verifying each candidate.

## Biological Relevance

Identifying approximate occurrences of short DNA motifs in larger genomic sequences is central to many analyses: 
- **Transcription factor binding sites** (e.g. locating candidate estrogen‐receptor motifs in regulatory regions)  
- **Short tandem repeat** (microsatellite) discovery  
- **Primer/probe specificity** checks in PCR and hybridization assays  
- **CRISPR off-target prediction**  

By reporting every fragment of the query that matches within a given tolerance, this tool can highlight potential regulatory elements, repetitive elements, or near-identical sequence regions with biological or clinical significance.

## Requirements

- **Python** 3.6 or higher  
- No external dependencies beyond the Python standard library  

## Installation

1. Download or clone this repository.  
2. Mark the script as executable (on Linux/macOS):
   ```bash
   chmod +x query_match.py

## Usage
Simply run the script without arguments to see a self-contained demonstration:

(_Set directory to project directory to run via bash_)

    python3 query_match.py

It will:

1. Generate a random “text” DNA sequence.
2. Embed a mutated fragment in a random “query” sequence.
3. Slide an n-length window over both strings, using l-mer filtration to find all matches with ≤ k mismatches.
4. Print the matching coordinate pairs (1-based indices).

## Example Output

    Searching for all 15-mers in query vs. text with ≤2 mismatches…
    
    Found matches at (query_start, text_start):
      q[12:27] ≈ t[78:93] → 2 mismatches
      q[33:48] ≈ t[150:165] → 1 mismatches


## Integrating into Other Workflows

To use custom sequences, import and call the query_match(q, t, n, k) function from within your own Python code:

    from query_match import query_match
    
    q = "ACGTACGT…"
    t = "TTGACGTA…"
    matches = query_match(q, t, n=15, k=2)
    # matches is a list of (query_start, text_start) pairs (1-based)

