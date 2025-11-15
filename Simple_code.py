import pandas as pd
import numpy as np

# --- DNA to Protein Translation Function ---
def dna_to_protein(dna_seq):
    """
    Translates a DNA sequence into a protein sequence.

    Assumptions:
    - DNA is a coding sequence (no introns).
    - Translation starts at the first nucleotide.
    - Codons not divisible by 3 at the end are ignored.
    - Unknown codons are translated as 'X'.

    Parameters:
    dna_seq (str): DNA sequence containing 'A', 'T', 'C', 'G'

    Returns:
    str: Protein sequence

    Raises:
    ValueError: If input is empty.
    """
    if not dna_seq:
        raise ValueError("DNA sequence must not be empty")

    dna_seq = dna_seq.upper()

    # Warn if sequence length is not divisible by 3
    if len(dna_seq) % 3 != 0:
        print("Warning: Sequence length not divisible by 3; last incomplete codon ignored.")

    # Standard codon table
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    protein_seq = ""
    # Translate codons into amino acids
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein_seq += codon_table.get(codon, 'X')  # 'X' for unknown codons

    return protein_seq


# --- Hamming Distance Function (for usernames) ---
def hamming_distance(str1, str2):
    """
    Calculates Hamming distance between two strings.
    Pads shorter string with spaces for comparison.

    Parameters:
    str1, str2 (str): Strings to compare

    Returns:
    int: Number of mismatched positions
    """
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


# --- Example DNA Sequences ---
dna_sequences = {
    "Seq1": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
    "Seq2": "ATGGCCATTGTAATGGGCCGCTGAAGGGCGCCCGATAG",
    "Seq3": "ATGGCCATTGTAATGGAACGCTGAAAGGGTGCCCGATAG"
}

# Warn if DNA sequences have inconsistent lengths
seq_lengths = [len(seq) for seq in dna_sequences.values()]
if len(set(seq_lengths)) != 1:
    print("Warning: DNA sequences have inconsistent lengths. Pairwise comparisons may be affected.")

# Translate DNA sequences to proteins
protein_sequences = {name: dna_to_protein(seq) for name, seq in dna_sequences.items()}

# Display DNA and Protein sequences in a pandas DataFrame
df_proteins = pd.DataFrame({
    "DNA Sequence Name": list(protein_sequences.keys()),
    "Protein Sequence": list(protein_sequences.values())
})
print("Protein Sequences Table:")
print(df_proteins, "\n")


# --- Hamming Distance Between Slack and Twitter Usernames ---
slack_username = "josoga"
twitter_username = "joseph"
username_distance = hamming_distance(slack_username, twitter_username)

# Display in pandas DataFrame
df_usernames = pd.DataFrame({
    "Slack Username": [slack_username],
    "Twitter Username": [twitter_username],
    "Hamming Distance": [username_distance]
})
print("Hamming Distance Table:")
print(df_usernames, "\n")


# --- Unit Tests ---
def test_dna_to_protein():
    assert dna_to_protein("ATG") == "M"
    assert dna_to_protein("ATGAAATAG") == "MK_"
    assert dna_to_protein("ATGNNN") == "MX"
    print("dna_to_protein tests passed.")

def test_hamming_distance():
    assert hamming_distance("AAAA", "AAAT") == 1
    assert hamming_distance("AAA", "AAAA") == 1
    assert hamming_distance("GTC", "GTC") == 0
    assert hamming_distance("josoga", "joseph") == 3
    print("hamming_distance tests passed.")

# Run unit tests
test_dna_to_protein()
test_hamming_distance()
Protein Sequences Table:
  DNA Sequence Name Protein Sequence
0              Seq1    MAIVMGR_KGAR_
1              Seq2     MAIVMGR_RAPD
2              Seq3    MAIVMER_KGAR_ 

Hamming Distance Table:
  Slack Username Twitter Username  Hamming Distance
0         josoga           joseph                 3 

dna_to_protein tests passed.
hamming_distance tests passed.
