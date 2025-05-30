# Query-Match
This utility locates all positions where any length-*n* fragment of a “query” DNA string aligns to a length-*n* fragment of a “text” DNA string with at most *k* mismatches. An efficient *l*-mer filtration step (with *l* = ⌊*n*/(*k*+1)⌋) is used to rapidly discard impossible alignments before verifying each candidate.
