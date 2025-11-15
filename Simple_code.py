import pandas as pd
import numpy as np

# --- Functions ---
def dna_to_protein(dna_seq):
    """
    Translates a DNA sequence into a protein sequence.
    
    Assumptions:
    - DNA is coding sequence (no introns).
    - Translation starts at the first nucleotide.
    - Codons not divisible by 3 at the end are ignored.
    - Unknown codons are translated as 'X'.
    """
    if not dna_seq:
        raise ValueError("DNA sequence must not be empty")
    if any(nuc not in "ATCGatcg" for nuc in dna_seq):
        raise ValueError("DNA sequence contains invalid characters")
    
    dna_seq = dna_seq.upper()
    
    if len(dna_seq) % 3 != 0:
        print("Warning: Sequence length not divisible by 3; last incomplete codon ignored.")
    
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    protein_seq = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein_seq += codon_table.get(codon, 'X')
    
    return protein_seq


def hamming_distance(str1, str2):
    """
    Calculates Hamming distance between two strings.
    Pads the shorter string with spaces.
    """
    if not str1 or not str2:
        raise ValueError("Input strings must not be empty")
    
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)
    
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


# --- Demo DNA Sequences ---
sequences = {
    "Seq1": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
    "Seq2": "ATGGCCATTGTAATGGGCCGCTGAAGGGCGCCCGATAG",
    "Seq3": "ATGGCCATTGTAATGGAACGCTGAAAGGGTGCCCGATAG"
}

# Translate DNA sequences to proteins
proteins = {name: dna_to_protein(seq) for name, seq in sequences.items()}

# Compute pairwise Hamming distances (DNA)
seq_names = list(sequences.keys())
dna_dist_matrix = np.zeros((len(seq_names), len(seq_names)))
protein_dist_matrix = np.zeros((len(seq_names), len(seq_names)))

for i, s1 in enumerate(seq_names):
    for j, s2 in enumerate(seq_names):
        dna_dist_matrix[i, j] = hamming_distance(sequences[s1], sequences[s2])
        protein_dist_matrix[i, j] = hamming_distance(proteins[s1], proteins[s2])

# Convert to DataFrames for nicer display
df_dna_dist = pd.DataFrame(dna_dist_matrix, index=seq_names, columns=seq_names)
df_prot_dist = pd.DataFrame(protein_dist_matrix, index=seq_names, columns=seq_names)
df_proteins = pd.DataFrame(list(proteins.items()), columns=["Sequence", "Protein"])

# --- Demo Hamming distance between usernames ---
usernames = {
    "Slack": "josoga",
    "Twitter": "joseph"
}

username_dist = hamming_distance(usernames["Slack"], usernames["Twitter"])

# --- Output Results ---
print("Protein Sequences:")
print(df_proteins.to_string(index=False))

print("\nPairwise Hamming Distances (DNA sequences):")
print(df_dna_dist)

print("\nPairwise Hamming Distances (Protein sequences):")
print(df_prot_dist)

print(f"\nHamming Distance between Slack and Twitter usernames: {username_dist}")
Warning: Sequence length not divisible by 3; last incomplete codon ignored.
Protein Sequences:
Sequence       Protein
    Seq1 MAIVMGR_KGAR_
    Seq2  MAIVMGR_RAPD
    Seq3 MAIVMER_KGAR_

Pairwise Hamming Distances (DNA sequences):
      Seq1  Seq2  Seq3
Seq1   0.0  10.0   2.0
Seq2  10.0   0.0  12.0
Seq3   2.0  12.0   0.0

Pairwise Hamming Distances (Protein sequences):
      Seq1  Seq2  Seq3
Seq1   0.0   5.0   1.0
Seq2   5.0   0.0   6.0
Seq3   1.0   6.0   0.0

Hamming Distance between Slack and Twitter usernames: 3

