from Bio import SeqIO

# ✅ Load sequences from a FASTA file
def load_fasta_sequences(fasta_path):
    return list(SeqIO.parse(fasta_path, "fasta"))

# ✅ Standard genetic code (codon → amino acid)
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# ✅ Check if two codons are synonymous
def is_synonymous(codon1, codon2):
    aa1 = codon_table.get(codon1.upper(), None)
    aa2 = codon_table.get(codon2.upper(), None)
    if aa1 is None or aa2 is None:
        return False  # Invalid codon
    return aa1 == aa2

# ✅ Optional: Validate codon format
def is_valid_codon(codon):
    return codon.upper() in codon_table

# ✅ Optional: Translate a codon to amino acid
def translate_codon(codon):
    return codon_table.get(codon.upper(), 'X')  # 'X' for unknown