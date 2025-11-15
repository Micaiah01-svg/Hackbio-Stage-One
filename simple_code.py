{
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/Micaiah01-svg/Hackbio-Stage-One/blob/main/simple_code.py\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "import pandas as pd\n",
        "import numpy as np\n",
        "\n",
        "# --- Functions ---\n",
        "def dna_to_protein(dna_seq):\n",
        "    \"\"\"\n",
        "    Translates a DNA sequence into a protein sequence.\n",
        "\n",
        "    Assumptions:\n",
        "    - DNA is coding sequence (no introns).\n",
        "    - Translation starts at the first nucleotide.\n",
        "    - Codons not divisible by 3 at the end are ignored.\n",
        "    - Unknown codons are translated as 'X'.\n",
        "    \"\"\"\n",
        "    if not dna_seq:\n",
        "        raise ValueError(\"DNA sequence must not be empty\")\n",
        "    if any(nuc not in \"ATCGatcg\" for nuc in dna_seq):\n",
        "        raise ValueError(\"DNA sequence contains invalid characters\")\n",
        "\n",
        "    dna_seq = dna_seq.upper()\n",
        "\n",
        "    if len(dna_seq) % 3 != 0:\n",
        "        print(\"Warning: Sequence length not divisible by 3; last incomplete codon ignored.\")\n",
        "\n",
        "    codon_table = {\n",
        "        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',\n",
        "        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',\n",
        "        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',\n",
        "        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',\n",
        "        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',\n",
        "        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',\n",
        "        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',\n",
        "        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',\n",
        "    }\n",
        "\n",
        "    protein_seq = \"\"\n",
        "    for i in range(0, len(dna_seq) - 2, 3):\n",
        "        codon = dna_seq[i:i+3]\n",
        "        protein_seq += codon_table.get(codon, 'X')\n",
        "\n",
        "    return protein_seq\n",
        "\n",
        "\n",
        "def hamming_distance(str1, str2):\n",
        "    \"\"\"\n",
        "    Calculates Hamming distance between two strings.\n",
        "    Pads the shorter string with spaces.\n",
        "    \"\"\"\n",
        "    if not str1 or not str2:\n",
        "        raise ValueError(\"Input strings must not be empty\")\n",
        "\n",
        "    max_len = max(len(str1), len(str2))\n",
        "    str1 = str1.ljust(max_len)\n",
        "    str2 = str2.ljust(max_len)\n",
        "\n",
        "    return sum(c1 != c2 for c1, c2 in zip(str1, str2))\n",
        "\n",
        "\n",
        "# --- Demo DNA Sequences ---\n",
        "sequences = {\n",
        "    \"Seq1\": \"ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\",\n",
        "    \"Seq2\": \"ATGGCCATTGTAATGGGCCGCTGAAGGGCGCCCGATAG\",\n",
        "    \"Seq3\": \"ATGGCCATTGTAATGGAACGCTGAAAGGGTGCCCGATAG\"\n",
        "}\n",
        "\n",
        "# Translate DNA sequences to proteins\n",
        "proteins = {name: dna_to_protein(seq) for name, seq in sequences.items()}\n",
        "\n",
        "# Compute pairwise Hamming distances (DNA)\n",
        "seq_names = list(sequences.keys())\n",
        "dna_dist_matrix = np.zeros((len(seq_names), len(seq_names)))\n",
        "protein_dist_matrix = np.zeros((len(seq_names), len(seq_names)))\n",
        "\n",
        "for i, s1 in enumerate(seq_names):\n",
        "    for j, s2 in enumerate(seq_names):\n",
        "        dna_dist_matrix[i, j] = hamming_distance(sequences[s1], sequences[s2])\n",
        "        protein_dist_matrix[i, j] = hamming_distance(proteins[s1], proteins[s2])\n",
        "\n",
        "# Convert to DataFrames for nicer display\n",
        "df_dna_dist = pd.DataFrame(dna_dist_matrix, index=seq_names, columns=seq_names)\n",
        "df_prot_dist = pd.DataFrame(protein_dist_matrix, index=seq_names, columns=seq_names)\n",
        "df_proteins = pd.DataFrame(list(proteins.items()), columns=[\"Sequence\", \"Protein\"])\n",
        "\n",
        "# --- Demo Hamming distance between usernames ---\n",
        "usernames = {\n",
        "    \"Slack\": \"josoga\",\n",
        "    \"Twitter\": \"joseph\"\n",
        "}\n",
        "\n",
        "username_dist = hamming_distance(usernames[\"Slack\"], usernames[\"Twitter\"])\n",
        "\n",
        "# --- Output Results ---\n",
        "print(\"Protein Sequences:\")\n",
        "print(df_proteins.to_string(index=False))\n",
        "\n",
        "print(\"\\nPairwise Hamming Distances (DNA sequences):\")\n",
        "print(df_dna_dist)\n",
        "\n",
        "print(\"\\nPairwise Hamming Distances (Protein sequences):\")\n",
        "print(df_prot_dist)\n",
        "\n",
        "print(f\"\\nHamming Distance between Slack and Twitter usernames: {username_dist}\")\n"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "2Mcav6aVa9PU",
        "outputId": "be82585c-68dd-4716-c70d-7976d4d2f293"
      },
      "execution_count": 14,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "Warning: Sequence length not divisible by 3; last incomplete codon ignored.\n",
            "Protein Sequences:\n",
            "Sequence       Protein\n",
            "    Seq1 MAIVMGR_KGAR_\n",
            "    Seq2  MAIVMGR_RAPD\n",
            "    Seq3 MAIVMER_KGAR_\n",
            "\n",
            "Pairwise Hamming Distances (DNA sequences):\n",
            "      Seq1  Seq2  Seq3\n",
            "Seq1   0.0  10.0   2.0\n",
            "Seq2  10.0   0.0  12.0\n",
            "Seq3   2.0  12.0   0.0\n",
            "\n",
            "Pairwise Hamming Distances (Protein sequences):\n",
            "      Seq1  Seq2  Seq3\n",
            "Seq1   0.0   5.0   1.0\n",
            "Seq2   5.0   0.0   6.0\n",
            "Seq3   1.0   6.0   0.0\n",
            "\n",
            "Hamming Distance between Slack and Twitter usernames: 3\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [],
      "metadata": {
        "id": "jjLbSR-EbonU"
      },
      "execution_count": null,
      "outputs": []
    }
  ],
  "metadata": {
    "colab": {
      "name": "Welcome to Colab",
      "provenance": [],
      "include_colab_link": true
    },
    "kernelspec": {
      "display_name": "Python 3",
      "name": "python3"
    }
  },
  "nbformat": 4,
  "nbformat_minor": 0
}