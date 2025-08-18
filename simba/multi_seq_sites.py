from itertools import combinations
import pandas as pd
from simba.utils import codon_table, is_synonymous
from simba.utils import load_fasta_sequences

def multi_seq_codon_map(sequences):
    seq_len_codons = len(sequences[0].seq) // 3
    data = [
        {"Codon_Position": i+1,
         "Syn_Sites": 0.0,
         "Nonsyn_Sites": 0.0,
         "dN_obs": 0,
         "dS_obs": 0,
         "Comparisons": 0}
        for i in range(seq_len_codons)
    ]
    
    for seq1, seq2 in combinations(sequences, 2):
        codons1 = [str(seq1.seq)[i:i+3] for i in range(0, len(seq1.seq), 3)]
        codons2 = [str(seq2.seq)[i:i+3] for i in range(0, len(seq2.seq), 3)]
        
        for idx, (c1, c2) in enumerate(zip(codons1, codons2)):
            if len(c1) != 3 or len(c2) != 3 or '-' in c1+c2:
                continue
            aa1 = codon_table.get(c1.upper(), None)
            if aa1 is None or aa1 == '*':
                continue

            syn_count, nonsyn_count = 0, 0
            for pos in range(3):
                for nt in ['A','T','C','G']:
                    if nt != c1[pos]:
                        mutant = c1[:pos] + nt + c1[pos+1:]
                        aa_mut = codon_table.get(mutant.upper(), None)
                        if aa_mut and aa_mut != '*':
                            if aa_mut == aa1:
                                syn_count += 1
                            else:
                                nonsyn_count += 1
            total_changes = syn_count + nonsyn_count
            if total_changes > 0:
                data[idx]["Syn_Sites"] += syn_count / total_changes
                data[idx]["Nonsyn_Sites"] += nonsyn_count / total_changes

            if c1 != c2:
                if is_synonymous(c1, c2):
                    data[idx]["dS_obs"] += 1
                else:
                    data[idx]["dN_obs"] += 1
            
            data[idx]["Comparisons"] += 1

    for row in data:
        if row["Comparisons"] > 0:
            row["Syn_Sites"] /= row["Comparisons"]
            row["Nonsyn_Sites"] /= row["Comparisons"]
    
    return pd.DataFrame(data)