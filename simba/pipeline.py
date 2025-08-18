from simba.report import (
    add_hotspot_flag,
    generate_mutation_inventory,
    plot_tracks,
    plot_domain_boxplot
)
from simba.utils import load_fasta_sequences
from simba.multi_seq_sites import multi_seq_codon_map

def run_pipeline(fasta_path, annotation_path, threshold=1.5, output_prefix="simba_output"):
     sequences = load_fasta_sequences(fasta_path)
     ref_seq = sequences[0]

     df_pop = multi_seq_codon_map(sequences)
     df_pop = add_hotspot_flag(df_pop, threshold)

     plot_tracks(df_pop, annotation_path, output_path=f"{output_prefix}_tracks.png")
     plot_domain_boxplot(df_pop, annotation_path, output_path=f"{output_prefix}_boxplot.png")

     hotspot_table = generate_mutation_inventory(df_pop, annotation_path, ref_seq, threshold)
     hotspot_table.to_csv(f"{output_prefix}_hotspots.csv", index=False)

     print(f"✅ SIMBA analysis complete. Outputs saved with prefix: {output_prefix}")