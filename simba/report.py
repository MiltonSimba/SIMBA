# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import kruskal
# from Bio.Seq import Seq
# from Bio.Data import CodonTable

# codon_table_std = CodonTable.unambiguous_dna_by_name["Standard"]

# def add_hotspot_flag(df_codon, threshold=1.0, eps=1e-6):
#     df = df_codon.copy()
#     df["dN_rate"] = df["dN_obs"] / (df["Nonsyn_Sites"] + eps)
#     df["dS_rate"] = df["dS_obs"] / (df["Syn_Sites"] + eps)
#     df["dNdS"] = df["dN_rate"] / (df["dS_rate"] + eps)
#     df["Hotspot"] = df["dNdS"] >= threshold
#     df["Selection"] = pd.cut(
#         df["dNdS"],
#         bins=[-float("inf"), 0.5, threshold, float("inf")],
#         labels=["Negative", "Neutral", "Positive"]
#     )
#     return df

# def map_domains(df_codon, annotation_file):
#     ann_df = pd.read_csv(annotation_file)
#     domain_map = []
#     for pos in df_codon["Codon_Position"]:
#         domains = ann_df[
#             (ann_df["Start_Codon"] <= pos) & (ann_df["End_Codon"] >= pos)
#         ]["Feature"].tolist()
#         domain_map.append("; ".join(domains) if domains else "Inter-domain")
#     df_codon["Domain"] = domain_map
#     return df_codon

# def generate_mutation_inventory(df_codon, annotation_file, ref_seq, threshold=1.0, eps=1e-6):
#     df = add_hotspot_flag(df_codon, threshold, eps)
#     df = map_domains(df, annotation_file)

#     ref_codons = [str(ref_seq.seq)[i:i+3] for i in range(0, len(ref_seq.seq), 3)]
#     df["Ref_Codon"] = [c if len(c) == 3 else "NNN" for c in ref_codons]
#     df["Ref_AA"] = [
#         Seq(c).translate(table=codon_table_std) if "N" not in c else "X"
#         for c in df["Ref_Codon"]
#     ]

#     hotspots = df[df["Selection"] == "Positive"].copy()
#     return hotspots[[
#         "Codon_Position", "Ref_Codon", "Ref_AA",
#         "dN_obs", "dS_obs", "Syn_Sites", "Nonsyn_Sites",
#         "dNdS", "Domain", "Selection"
#     ]].sort_values(by="dNdS", ascending=False)

# def plot_tracks(df_codon, annotation_file, output_path="pop_codon_tracks_hotspots.png"):
#     df_codon = map_domains(df_codon, annotation_file)
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14,8),
#                                    sharex=True, gridspec_kw={"height_ratios": [3, 0.6]})

#     ax1.bar(df_codon["Codon_Position"], df_codon["dN_obs"], color='red', alpha=0.6, label="Total dN")
#     ax1.bar(df_codon["Codon_Position"], df_codon["dS_obs"], bottom=df_codon["dN_obs"], 
#             color='blue', alpha=0.6, label="Total dS")
#     ax1.plot(df_codon["Codon_Position"], df_codon["Syn_Sites"], color='cyan', marker='o', label="Mean Syn Sites")
#     ax1.plot(df_codon["Codon_Position"], df_codon["Nonsyn_Sites"], color='orange', marker='o', label="Mean Nonsyn Sites")
#     for idx, row in df_codon[df_codon["Hotspot"]].iterrows():
#         ax1.axvline(row["Codon_Position"], color='gold', linestyle='--', alpha=0.8)
#     ax1.set_ylabel("Counts / Sites")
#     ax1.set_title("Population-level Per-Codon dN/dS with Domain Tracks & Hotspots")
#     ax1.legend()

#     ann_df = pd.read_csv(annotation_file)
#     ax2.set_ylim(0, 1)
#     ax2.axis("off")
#     for _, row in ann_df.iterrows():
#         ax2.axvspan(row["Start_Codon"], row["End_Codon"], color='grey', alpha=0.3)
#         ax2.text((row["Start_Codon"] + row["End_Codon"]) / 2, 0.5, 
#                  row["Feature"], ha='center', va='center', fontsize=8)
#     ax2.set_xlabel("Codon Position")

#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300)
#     plt.show()

# def plot_domain_boxplot(df_codon, annotation_file, output_path="domain_dn_ds_boxplots.png"):
#     df_codon = map_domains(df_codon, annotation_file)
#     pval = kruskal(*[
#         df_codon[df_codon["Domain"] == d]["dNdS"].dropna()
#         for d in df_codon["Domain"].unique()
#     ])[1]

#     plt.figure(figsize=(10,6))
#     df_codon.boxplot(column="dNdS", by="Domain", grid=False, vert=True)
#     plt.axhline(y=1, color='red', linestyle='--', alpha=0.7)
#     plt.xticks(rotation=45, ha='right')
#     plt.ylabel("dN/dS Ratio")
#     plt.title(f"Domain-wise dN/dS Distribution\nKruskal-Wallis p = {pval:.2e}")
#     plt.suptitle("")
#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300)
#     plt.show()


# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.stats import kruskal
# from Bio.Seq import Seq
# from Bio.Data import CodonTable

# codon_table_std = CodonTable.unambiguous_dna_by_name["Standard"]

# def add_hotspot_flag(df_codon, threshold=1.0, eps=1e-6):
#     df = df_codon.copy()
#     df["dN_rate"] = df["dN_obs"] / (df["Nonsyn_Sites"] + eps)
#     df["dS_rate"] = df["dS_obs"] / (df["Syn_Sites"] + eps)
#     df["dNdS"] = df["dN_rate"] / (df["dS_rate"] + eps)
#     df["Hotspot"] = df["dNdS"] >= threshold
#     df["Selection"] = pd.cut(
#         df["dNdS"],
#         bins=[-float("inf"), 0.5, threshold, float("inf")],
#         labels=["Negative", "Neutral", "Positive"]
#     )
#     return df

# def map_domains(df_codon, annotation_file):
#     ann_df = pd.read_csv(annotation_file)
#     domain_map = []
#     for pos in df_codon["Codon_Position"]:
#         domains = ann_df[
#             (ann_df["Start_Codon"] <= pos) & (ann_df["End_Codon"] >= pos)
#         ]["Feature"].tolist()
#         domain_map.append("; ".join(domains) if domains else "Inter-domain")
#     df_codon["Domain"] = domain_map
#     return df_codon

# def generate_mutation_inventory(df_codon, annotation_file, ref_seq, threshold=1.0, eps=1e-6):
#     df = add_hotspot_flag(df_codon, threshold, eps)
#     df = map_domains(df, annotation_file)

#     ref_codons = [str(ref_seq.seq)[i:i+3] for i in range(0, len(ref_seq.seq), 3)]
#     df["Ref_Codon"] = [c if len(c) == 3 else "NNN" for c in ref_codons]
#     df["Ref_AA"] = [
#         Seq(c).translate(table=codon_table_std) if "N" not in c else "X"
#         for c in df["Ref_Codon"]
#     ]

#     hotspots = df[df["Selection"] == "Positive"].copy()
#     return hotspots[[
#         "Codon_Position", "Ref_Codon", "Ref_AA",
#         "dN_obs", "dS_obs", "Syn_Sites", "Nonsyn_Sites",
#         "dNdS", "Domain", "Selection"
#     ]].sort_values(by="dNdS", ascending=False)

# def plot_tracks(df_codon, annotation_file, output_path="pop_codon_tracks_hotspots.png"):
#     df_codon = map_domains(df_codon, annotation_file)
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14,8),
#                                    sharex=True, gridspec_kw={"height_ratios": [3, 0.6]})

#     # Top panel: bars + lines + hotspot markers
#     ax1.bar(df_codon["Codon_Position"], df_codon["dN_obs"], color='red', alpha=0.6, label="Total dN")
#     ax1.bar(df_codon["Codon_Position"], df_codon["dS_obs"], bottom=df_codon["dN_obs"], 
#             color='blue', alpha=0.6, label="Total dS")
#     ax1.plot(df_codon["Codon_Position"], df_codon["Syn_Sites"], color='cyan', marker='o', label="Mean Syn Sites")
#     ax1.plot(df_codon["Codon_Position"], df_codon["Nonsyn_Sites"], color='orange', marker='o', label="Mean Nonsyn Sites")
#     for idx, row in df_codon[df_codon["Hotspot"]].iterrows():
#         ax1.axvline(row["Codon_Position"], color='gold', linestyle='--', alpha=0.8)
#     ax1.set_ylabel("Counts / Sites")
#     ax1.set_title("Population-level Per-Codon dN/dS with Domain Tracks & Hotspots")
#     ax1.legend()

#     # Bottom panel: domain track
#     ann_df = pd.read_csv(annotation_file)

#     # Short label mapping
#     short_labels = {
#         'Signal peptide': 'Sig pep',
#         'Subdomain 1': 'SD1',
#         'Subdomain 2': 'SD2',
#         'Fusion peptide': 'FP',
#         'Heptad repeat 1 (HR1)': 'HR1',
#         'Heptad repeat 2 (HR2)': 'HR2',
#         'Central helix': 'CH',
#         'Connector domain': 'CD',
#         'Transmembrane domain': 'TM',
#         'Cytoplasmic tail': 'CT'
#     }
#     ann_df['Plot_Label'] = ann_df['Feature'].map(lambda x: short_labels.get(x, x))

#     # Colour palette for alternating domains
#     domain_colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#2b8cbe']

#     ax2.set_ylim(0, 1)
#     ax2.axis("off")
#     for i, row in ann_df.iterrows():
#         ax2.axvspan(row["Start_Codon"], row["End_Codon"],
#                     color=domain_colors[i % len(domain_colors)], alpha=0.3)
#         ax2.text((row["Start_Codon"] + row["End_Codon"]) / 2, 0.5, 
#                  row["Plot_Label"], ha='center', va='center', fontsize=8)
#     ax2.set_xlabel("Codon Position")

#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300)
#     plt.show()

# def plot_domain_boxplot(df_codon, annotation_file, output_path="domain_dn_ds_boxplots.png"):
#     df_codon = map_domains(df_codon, annotation_file)
#     pval = kruskal(*[
#         df_codon[df_codon["Domain"] == d]["dNdS"].dropna()
#         for d in df_codon["Domain"].unique()
#     ])[1]

#     plt.figure(figsize=(10,6))
#     df_codon.boxplot(column="dNdS", by="Domain", grid=False, vert=True)
#     plt.axhline(y=1, color='red', linestyle='--', alpha=0.7)
#     plt.xticks(rotation=45, ha='right')
#     plt.ylabel("dN/dS Ratio")
#     plt.title(f"Domain-wise dN/dS Distribution\nKruskal-Wallis p = {pval:.2e}")
#     plt.suptitle("")
#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300)
#     plt.show()   


# import os
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd

# # ==============================
# # Domain-wise dN/dS Plot
# # ==============================
# def plot_domainwise_dn_ds(df, p_value=None, output_path=None,
#                           threshold_line=1.0,
#                           figsize=(10, 6)):
#     """
#     Generates a domain-wise dN/dS plot with:
#       - Semi-transparent boxplots showing median + IQR
#       - Overlaid raw data points
#       - Optional horizontal threshold line
#       - Optional Kruskal–Wallis p-value in title
#     """
#     plt.figure(figsize=figsize)

#     # Summary stats
#     sns.boxplot(
#         data=df, x="Domain", y="dN_dS",
#         whis=1.5, showfliers=False,
#         boxprops={'facecolor': 'lightgray', 'alpha': 0.5},
#         medianprops={'color': 'blue', 'linewidth': 2}
#     )

#     # Raw points
#     sns.stripplot(
#         data=df, x="Domain", y="dN_dS",
#         jitter=True, color='black', alpha=0.7
#     )

#     # Threshold
#     if threshold_line is not None:
#         plt.axhline(threshold_line, color='red', linestyle='--', linewidth=1)

#     # Labels/title
#     plt.xticks(rotation=45, ha='right')
#     title = "Domain-wise dN/dS Distribution"
#     if p_value is not None:
#         title += f"\nKruskal–Wallis p = {p_value:.2e}"
#     plt.title(title, fontsize=14)
#     plt.ylabel("dN/dS Ratio")

#     plt.tight_layout()

#     if output_path:
#         plt.savefig(output_path, dpi=300)
#         plt.close()
#     else:
#         plt.show()


# # ==============================
# # Report Generation
# # ==============================
# def generate_report(domainwise_df, p_value, output_dir):
#     """
#     Creates all final report plots and saves them to the output directory.
#     """
#     os.makedirs(output_dir, exist_ok=True)

#     # Domain-wise plot
#     domain_plot_path = os.path.join(output_dir, "domainwise_dn_ds.png")
#     plot_domainwise_dn_ds(
#         df=domainwise_df,
#         p_value=p_value,
#         output_path=domain_plot_path
#     )

#     # Additional plots / report sections can be added here
#     # e.g., per-site dN/dS plots, histograms, summary tables, etc.


# # ==============================
# # Utility: Save Summary Table
# # ==============================
# def save_summary_table(df, output_path):
#     """
#     Saves a given DataFrame as a CSV file for the report.
#     """
#     df.to_csv(output_path, index=False)


# # ==============================
# # Main Entry for Testing
# # ==============================
# if __name__ == "__main__":
#     # Example: Dummy run for testing layout
#     example_df = pd.DataFrame({
#         "Domain": ["NTD", "RBD", "HR1", "NTD", "RBD", "HR1"],
#         "dN_dS": [0.5, 2.1, 1.2, 0.8, 3.5, 0.9]
#     })
#     generate_report(example_df, p_value=2.03e-4, output_dir="test_report")
#     save_summary_table(example_df, "test_report/domainwise_table.csv")
#     print("Test report generated.")




import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from Bio.Seq import Seq
from Bio.Data import CodonTable

codon_table_std = CodonTable.unambiguous_dna_by_name["Standard"]

def add_hotspot_flag(df_codon, threshold=1.0, eps=1e-6):
    df = df_codon.copy()
    df["dN_rate"] = df["dN_obs"] / (df["Nonsyn_Sites"] + eps)
    df["dS_rate"] = df["dS_obs"] / (df["Syn_Sites"] + eps)
    df["dNdS"] = df["dN_rate"] / (df["dS_rate"] + eps)
    df["Hotspot"] = df["dNdS"] >= threshold
    df["Selection"] = pd.cut(
        df["dNdS"],
        bins=[-float("inf"), 0.5, threshold, float("inf")],
        labels=["Negative", "Neutral", "Positive"]
    )
    return df

def map_domains(df_codon, annotation_file):
    ann_df = pd.read_csv(annotation_file)
    domain_map = []
    for pos in df_codon["Codon_Position"]:
        domains = ann_df[
            (ann_df["Start_Codon"] <= pos) & (ann_df["End_Codon"] >= pos)
        ]["Feature"].tolist()
        domain_map.append("; ".join(domains) if domains else "Inter-domain")
    df_codon["Domain"] = domain_map
    return df_codon

def generate_mutation_inventory(df_codon, annotation_file, ref_seq, threshold=1.0, eps=1e-6):
    df = add_hotspot_flag(df_codon, threshold, eps)
    df = map_domains(df, annotation_file)

    ref_codons = [str(ref_seq.seq)[i:i+3] for i in range(0, len(ref_seq.seq), 3)]
    df["Ref_Codon"] = [c if len(c) == 3 else "NNN" for c in ref_codons]
    df["Ref_AA"] = [
        Seq(c).translate(table=codon_table_std) if "N" not in c else "X"
        for c in df["Ref_Codon"]
    ]

    hotspots = df[df["Selection"] == "Positive"].copy()
    return hotspots[[
        "Codon_Position", "Ref_Codon", "Ref_AA",
        "dN_obs", "dS_obs", "Syn_Sites", "Nonsyn_Sites",
        "dNdS", "Domain", "Selection"
    ]].sort_values(by="dNdS", ascending=False)

def plot_tracks(df_codon, annotation_file, output_path="pop_codon_tracks_hotspots.png"):
    df_codon = map_domains(df_codon, annotation_file)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14,8),
                                   sharex=True, gridspec_kw={"height_ratios": [3, 0.6]})

    # Top panel: bars + lines + hotspot markers
    ax1.bar(df_codon["Codon_Position"], df_codon["dN_obs"], color='red', alpha=0.6, label="Total dN")
    ax1.bar(df_codon["Codon_Position"], df_codon["dS_obs"], bottom=df_codon["dN_obs"], 
            color='blue', alpha=0.6, label="Total dS")
    ax1.plot(df_codon["Codon_Position"], df_codon["Syn_Sites"], color='cyan', marker='o', label="Mean Syn Sites")
    ax1.plot(df_codon["Codon_Position"], df_codon["Nonsyn_Sites"], color='orange', marker='o', label="Mean Nonsyn Sites")
    for idx, row in df_codon[df_codon["Hotspot"]].iterrows():
        ax1.axvline(row["Codon_Position"], color='gold', linestyle='--', alpha=0.8)
    ax1.set_ylabel("Counts / Sites")
    ax1.set_title("Population-level Per-Codon dN/dS with Domain Tracks & Hotspots")
    ax1.legend()

    # Bottom panel: domain track
    ann_df = pd.read_csv(annotation_file)

    # Short label mapping
    short_labels = {
        'Signal peptide': 'Sig pep',
        'Subdomain 1': 'SD1',
        'Subdomain 2': 'SD2',
        'Fusion peptide': 'FP',
        'Heptad repeat 1 (HR1)': 'HR1',
        'Heptad repeat 2 (HR2)': 'HR2',
        'Central helix': 'CH',
        'Connector domain': 'CD',
        'Transmembrane domain': 'TM',
        'Cytoplasmic tail': 'CT'
    }
    ann_df['Plot_Label'] = ann_df['Feature'].map(lambda x: short_labels.get(x, x))

    # Colour palette for alternating domains
    domain_colors = ['#f0f9e8', '#bae4bc', '#7bccc4', '#2b8cbe']

    ax2.set_ylim(0, 1)
    ax2.axis("off")
    for i, row in ann_df.iterrows():
        ax2.axvspan(row["Start_Codon"], row["End_Codon"],
                    color=domain_colors[i % len(domain_colors)], alpha=0.3)
        ax2.text((row["Start_Codon"] + row["End_Codon"]) / 2, 0.5, 
                 row["Plot_Label"], ha='center', va='center', fontsize=8)
    ax2.set_xlabel("Codon Position")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

def plot_domain_boxplot(df_codon, annotation_file, output_path="domain_dn_ds_boxplots.png"):
    """
    Updated: now uses seaborn for semi-transparent boxplots + overlaid raw data points.
    """
    df_codon = map_domains(df_codon, annotation_file)
    pval = kruskal(*[
        df_codon[df_codon["Domain"] == d]["dNdS"].dropna()
        for d in df_codon["Domain"].unique()
    ])[1]

    plt.figure(figsize=(10,6))

    # Boxplot layer
    sns.boxplot(
        data=df_codon, x="Domain", y="dNdS",
        whis=1.5, showfliers=False,
        boxprops={'facecolor': 'lightgray', 'alpha': 0.5},
        medianprops={'color': 'blue', 'linewidth': 2}
    )

    # Raw points
    sns.stripplot(
        data=df_codon, x="Domain", y="dNdS",
        jitter=True, color='black', alpha=0.7
    )

    plt.axhline(y=1, color='red', linestyle='--', alpha=0.7)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("dN/dS Ratio")
    plt.title(f"Domain-wise dN/dS Distribution\nKruskal-Wallis p = {pval:.2e}")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
