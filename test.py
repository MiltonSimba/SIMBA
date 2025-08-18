from simba.utils import codon_table, is_synonymous

print(codon_table['ATG'])  # Should print 'M'
print(is_synonymous('GAA', 'GAG'))  # Should print True
print(is_synonymous('GAA', 'GAC'))  # Should print False
