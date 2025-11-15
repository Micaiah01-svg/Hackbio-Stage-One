# DNA to Protein Translation
def dna_to_protein(dna_seq):
    """
    Translates a DNA sequence into a protein sequence.
    
    Parameters:
    dna_seq (str): DNA sequence consisting of 'A', 'T', 'C', 'G'
    
    Returns:
    str: Protein sequence using the standard genetic code
    """
    # Validate input
    if not dna_seq:
        raise ValueError("DNA sequence must not be empty")
    if any(nuc not in "ATCGatcg" for nuc in dna_seq):
        raise ValueError("DNA sequence contains invalid characters")
    
    dna_seq = dna_seq.upper()
    
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
    
    # Translate codons
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein_seq += codon_table.get(codon, 'X')  # X for unknown codon
    
    return protein_seq


# Hamming Distance
def hamming_distance(str1, str2):
    """
    Calculates the Hamming distance between two strings.
    Pads the shorter string with spaces if lengths are unequal.
    
    Parameters:
    str1 (str): First string (e.g., Slack username)
    str2 (str): Second string (e.g., Twitter/X handle)
    
    Returns:
    int: Hamming distance
    """
    if not str1 or not str2:
        raise ValueError("Input strings must not be empty")
    
    # Pad shorter string with spaces
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)
    
    # Calculate mismatches
    mismatches = sum(c1 != c2 for c1, c2 in zip(str1, str2))
    return mismatches


# Example Usage:
dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
print("Protein:", dna_to_protein(dna_seq))

slack_username = "josoga"
twitter_handle = "joseph"
print("Hamming Distance:", hamming_distance(slack_username, twitter_handle))
Protein: MAIVMGR_KGAR_
Hamming Distance: 3

