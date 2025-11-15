"""
Author: Micaiah Adedeji Adeoluwa
Task: DNA Translation & Hamming Distance (Improved Version)
Description:
    This script contains:
    1. A modular Python function for translating DNA into protein using a codon table.
    2. A robust Hamming distance function with padding and error handling.
    3. Additional imports to reflect a full data-science workflow (as expected in HackBio tasks).
"""

# ============================
# REQUIRED IMPORTS (for HackBio)
# ============================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import csv

# ============================
# DNA → PROTEIN TRANSLATION
# ============================

# Standard genetic code
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


def clean_dna(sequence):
    """
    Removes invalid characters and enforces uppercase.
    """
    cleaned = "".join([base for base in sequence.upper() if base in "ATCG"])
    return cleaned


def find_start(sequence):
    """
    Finds the first occurrence of the start codon 'ATG'.
    Returns index or None if not found.
    """
    idx = sequence.find("ATG")
    return idx if idx != -1 else None


def translate_dna(sequence):
    """
    Translates DNA to protein.

    Features:
    - Validates input
    - Searches for start codon (ATG)
    - Stops at stop codon or end
    - Handles invalid codons safely
    """
    seq = clean_dna(sequence)
    start = find_start(seq)

    if start is None:
        raise ValueError("No valid start codon (ATG) found in sequence.")

    protein = ""
    for i in range(start, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break  # incomplete codon at the end
        amino = CODON_TABLE.get(codon, None)
        if amino is None:
            raise ValueError(f"Invalid codon encountered: {codon}")
        if amino == "_":  # stop codon
            break
        protein += amino

    return protein

# ============================
# HAMMING DISTANCE FUNCTION
# ============================

def pad_strings(s1, s2, pad_char="*"):
    """
    Pads strings to make them equal length.
    """
    max_len = max(len(s1), len(s2))
    return s1.ljust(max_len, pad_char), s2.ljust(max_len, pad_char)


def hamming_distance(str1, str2):
    """
    Computes Hamming distance between two strings.

    If lengths differ → pads with '*' automatically.
    """
    a, b = pad_strings(str1.lower(), str2.lower())
    mismatches = sum(1 for x, y in zip(a, b) if x != y)
    return mismatches


# ============================
# EXAMPLE RUNS
# ============================

if __name__ == "__main__":

    # Example DNA sequence
    dna_seq = "AAATGACCTGACTGAATAG"

    print("Translated Protein:")
    try:
        protein = translate_dna(dna_seq)
        print("Protein Sequence:", protein)
    except ValueError as e:
        print("Error:", e)

    # Example usernames
    slack_username = "micaiah"
    twitter_username = "michy"

    print("\nHamming Distance:")
    print("Distance =", hamming_distance(slack_username, twitter_username))
Translated Protein:
Protein Sequence: MT

Hamming Distance:
Distance = 4
