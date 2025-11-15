import pandas as pd
import itertools
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# --- Standard Codon Table ---
CODON_TABLE = {
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

# --- Ambiguous Nucleotide Map ---
AMBIGUOUS_MAP = {
    'N': ['A','T','C','G'],  # N can be any nucleotide
    'R': ['A','G'],           # purine
    'Y': ['C','T'],           # pyrimidine
    'S': ['G','C'],
    'W': ['A','T'],
    'K': ['G','T'],
    'M': ['A','C'],
    'B': ['C','G','T'],
    'D': ['A','G','T'],
    'H': ['A','C','T'],
    'V': ['A','C','G']
}

# --- DNA to Protein Translation Function ---
def dna_to_protein(dna_seq, max_seq_len=10_000):
    """
    Translate DNA sequence into protein sequence, considering ambiguous nucleotides.
    
    Parameters:
    dna_seq (str): DNA sequence with standard or ambiguous nucleotides
    max_seq_len (int): Maximum allowed sequence length to avoid performance issues
    
    Returns:
    str: Protein sequence (ambiguous codons listed as combinations or 'X')
    """
    if not dna_seq:
        raise ValueError("DNA sequence must not be empty")
    if len(dna_seq) > max_seq_len:
        raise ValueError(f"Sequence too long (>{max_seq_len} bp), processing aborted.")
    
    dna_seq = dna_seq.upper()
    
    # Check for completely invalid nucleotides
    valid_bases = set("ATCG" + "".join(AMBIGUOUS_MAP.keys()))
    invalid_bases = set(dna_seq) - valid_bases
    if invalid_bases:
        raise ValueError(f"Invalid nucleotides detected: {invalid_bases}")
    
    protein_seq = ""
    for i in range(0, len(dna_seq)-2, 3):
        codon = dna_seq[i:i+3]
        if any(base in AMBIGUOUS_MAP for base in codon):
            # Expand ambiguous codons to all possible combinations
            possible_bases = [AMBIGUOUS_MAP.get(b, [b]) for b in codon]
            codon_combinations = [''.join(p) for p in itertools.product(*possible_bases)]
            amino_acids = set(CODON_TABLE.get(c, 'X') for c in codon_combinations)
            if len(amino_acids) == 1:
                protein_seq += amino_acids.pop()
            else:
                protein_seq += 'X'  # ambiguous translation
            logging.info(f"Ambiguous codon '{codon}' translated as '{protein_seq[-1]}'")
        else:
            protein_seq += CODON_TABLE.get(codon, 'X')
    
    return protein_seq


# --- Hamming Distance Function ---
def hamming_distance(str1, str2):
    """
    Calculates Hamming distance between two strings.
    Pads shorter string with spaces.
    """
    max_len = max(len(str1), len(str2))
    str1, str2 = str1.ljust(max_len), str2.ljust(max_len)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


# --- Translate Multiple Sequences ---
def translate_sequences(seq_dict):
    """
    Translates a dictionary of DNA sequences into proteins.
    """
    seq_lengths = [len(seq) for seq in seq_dict.values()]
    if len(set(seq_lengths)) != 1:
        logging.warning("Sequences have inconsistent lengths.")
    return {name: dna_to_protein(seq) for name, seq in seq_dict.items()}


# --- Example DNA Sequences ---
dna_sequences = {
    "Seq1": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
    "Seq2": "ATGGCCATTGTAATGGGCCGCTGAAGGGCGCCCGATAG",
    "Seq3": "ATGNNNATTGTAATGGNCCGCTGAAAGGGTGCCCGATAG"  # ambiguous codons
}

protein_sequences = translate_sequences(dna_sequences)

# Display in pandas
df_proteins = pd.DataFrame({
    "DNA Sequence Name": list(protein_sequences.keys()),
    "Protein Sequence": list(protein_sequences.values())
})
print("Protein Sequences Table:")
print(df_proteins, "\n")


# --- Hamming Distance Example ---
slack_username = "josoga"
twitter_username = "joseph"
distance = hamming_distance(slack_username, twitter_username)

df_usernames = pd.DataFrame({
    "Slack Username": [slack_username],
    "Twitter Username": [twitter_username],
    "Hamming Distance": [distance]
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

test_dna_to_protein()
test_hamming_distance()

Protein Sequences Table:
  DNA Sequence Name Protein Sequence
0              Seq1    MAIVMGR_KGAR_
1              Seq2     MAIVMGR_RAPD
2              Seq3    MXIVMXR_KGAR_ 

Hamming Distance Table:
  Slack Username Twitter Username  Hamming Distance
0         josoga           joseph                 3 

dna_to_protein tests passed.
hamming_distance tests passed.
