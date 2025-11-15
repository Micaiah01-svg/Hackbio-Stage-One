import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np


# =============================================================
#                        CODON TABLE
# =============================================================

CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_',
    'TGC':'C','TGT':'C','TGA':'_','TGG':'W',
}


# =============================================================
#              DNA CLEANING & START CODON SEARCH
# =============================================================

def clean_dna(seq):
    """Cleans DNA sequence by keeping only valid nucleotides."""
    return "".join([base for base in seq.upper() if base in "ATCG"])


def find_start(seq):
    """Finds first ATG start codon."""
    idx = seq.find("ATG")
    return idx if idx != -1 else None


# =============================================================
#                 DNA → PROTEIN TRANSLATION
# =============================================================

def translate_dna(sequence):
    seq = clean_dna(sequence)
    start = find_start(seq)

    if start is None:
        return ""

    protein = ""

    for i in range(start, len(seq), 3):
        codon = seq[i:i+3]

        if len(codon) < 3:
            break

        aa = CODON_TABLE.get(codon)

        if aa is None:
            continue

        if aa == "_":  # stop codon
            break

        protein += aa

    return protein


# =============================================================
#                      HAMMING DISTANCE
# =============================================================

def pad_strings(a, b, pad="*"):
    max_len = max(len(a), len(b))
    return a.ljust(max_len, pad), b.ljust(max_len, pad)


def hamming_distance(s1, s2):
    a, b = pad_strings(s1.lower(), s2.lower())
    return sum(1 for x, y in zip(a, b) if x != y)


# =============================================================
#                      UNIT TESTS
# =============================================================

def run_tests():
    print("Running tests...")

    assert translate_dna("ATGACCTGA") == "MT"
    assert translate_dna("CCCCCC") == ""
    assert hamming_distance("abc", "abc") == 0
    assert hamming_distance("micaiah", "michy") > 0

    print("All tests passed!\n")


# =============================================================
#                        MAIN PROGRAM
# =============================================================

if __name__ == "__main__":
    run_tests()

    dna = "ATGACCTGACTGAATAG"
    protein = translate_dna(dna)
    print("DNA Sequence:", dna)
    print("Protein Translation:", protein)

    slack = "micaiah"
    twitter = "michy"
    distance = hamming_distance(slack, twitter)
    print("\nSlack Username:", slack)
    print("Twitter Handle:", twitter)
    print("Hamming Distance:", distance)
DNA Sequence: ATGACCTGACTGAATAG
Protein Translation: MT

Slack Username: micaiah
Twitter Handle: michy
Hamming Distance: 4
