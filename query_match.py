#!/usr/bin/env python3
"""
query_match.py

Query Matching

Given:
  – A DNA text t of length m
  – A DNA query q of length p
  – An integer n ≤ p (length of the substrings to match)
  – A maximum Hamming distance k

Output:
  All pairs (i, j) such that the n‐mer q[i:i+n] and t[j:j+n] differ in at most k positions.

We use l-mer filtration where l = floor(n / (k+1)), per Theorem 9.1.
"""

import random
import string

def generate_random_dna(length: int) -> str:
    """Generate a random DNA string of given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def query_match(q: str, t: str, n: int, k: int):
    """
    Find all (i, j) with 0 ≤ i ≤ len(q)-n, 0 ≤ j ≤ len(t)-n,
    such that Hamming(q[i:i+n], t[j:j+n]) ≤ k.
    Uses l-mer filtration: any two strings within k mismatches
    must share an exact l-mer for l = n//(k+1).
    """
    p, m = len(q), len(t)
    if n > p or n > m:
        return []

    # choose l so that any n-mer with ≤k mismatches shares at least one common l-mer
    l = n // (k + 1)
    matches = set()

    # Precompute all l-mers for each query‐n-mer
    # For each starting position i in q, build a hash of its l-mers
    for i in range(p - n + 1):
        block = q[i:i+n]
        lmer_index = {}
        # map each l-mer in block to its offset r within the block
        for r in range(n - l + 1):
            km = block[r:r+l]
            lmer_index.setdefault(km, []).append(r)

        # scan text for any of these l-mers
        for j in range(m - l + 1):
            km = t[j:j+l]
            if km in lmer_index:
                for r in lmer_index[km]:
                    # candidate full‐length alignment starts at text position j-r
                    j0 = j - r
                    if j0 < 0 or j0 + n > m:
                        continue
                    # verify full n-mer with Hamming distance
                    segment_t = t[j0:j0+n]
                    segment_q = block
                    mismatches = sum(1 for a,b in zip(segment_q, segment_t) if a != b)
                    if mismatches <= k:
                        # record 1‐based positions
                        matches.add((i+1, j0+1))

    return sorted(matches)


if __name__ == "__main__":
    # Example usage
    m = 200                 # length of text
    p = 60                  # length of query
    n = 15                  # length of each match
    k = 2                   # max mismatches allowed

    # generate random text, pick a query (here we embed part of text + noise)
    text = generate_random_dna(m)
    # to guarantee some hits, pick a random substring and mutate up to k positions
    pos = random.randint(0, m - n)
    query_seed = list(text[pos:pos+n])
    for _ in range(k):
        idx = random.randrange(n)
        query_seed[idx] = random.choice([x for x in 'ACGT' if x != query_seed[idx]])
    # make the full query longer than n by flanking random DNA
    query = generate_random_dna(p - n)[:(pos % (p-n+1))] + ''.join(query_seed) \
            + generate_random_dna(p - n)[:(p-n) - (pos % (p-n+1))]

    print(f"Text (len={m}): {text}")
    print(f"Query (len={p}): {query}")
    print(f"Searching for all {n}-mers in query vs. text with ≤{k} mismatches...\n")

    hits = query_match(query, text, n, k)
    if hits:
        print("Found matches at (query_start, text_start):")
        for qi, tj in hits:
            print(f"  q[{qi}:{qi+n}] ≈ t[{tj}:{tj+n}]  "
                  f"→ {sum(1 for a,b in zip(query[qi-1:qi-1+n], text[tj-1:tj-1+n]) if a!=b)} mismatches")
    else:
        print("No matches found.")
