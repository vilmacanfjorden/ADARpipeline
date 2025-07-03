import pandas as pd
import matplotlib.pyplot as plt
import venn

def extract_position_and_gene(file_path):
    with open(file_path, 'r') as file:
        return [(line.split('\t')[1], line.split('\t')[5]) for line in file]

def get_position_to_gene_map(*files):
    return [{pos: gene for pos, gene in extract_position_and_gene(file)} for file in files]

def get_common_positions(*files):
    position_sets = [set(pos for pos, _ in extract_position_and_gene(file)) for file in files]
    return set.intersection(*position_sets)

def plot_venn_and_get_treatment_specific_positions(treatment_positions, dmso_positions, treatment_label):
    plt.figure(figsize=(8, 8))
    venn.venn({
        f"3'UTR {treatment_label}": treatment_positions,
        "3'UTR DMSO": dmso_positions
    })
    plt.title(f"Venn Diagram: {treatment_label} vs DMSO")
    plt.show()
    return treatment_positions - dmso_positions

def extract_ad_af_from_vcf(vcf_file, positions_of_interest, gene_map):
    positions_set = set(positions_of_interest)
    all_records = []

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = cols[1]
            vcf_pos = f"{chrom}:{pos}-{pos}"

            if vcf_pos not in positions_set:
                continue

            format_fields = cols[8].split(":")
            sample_fields = cols[9:]

            try:
                ad_index = format_fields.index("AD")
            except ValueError:
                continue  # Skip if AD not found

            gene = gene_map.get(vcf_pos, "Unknown")
            record = {
                "Position": vcf_pos,
                "GeneName": gene
            }

            for i, sample in enumerate(sample_fields):
                values = sample.split(":")
                try:
                    ad = values[ad_index]
                    ref_count, alt_count = map(int, ad.split(",")[:2])
                    total = ref_count + alt_count
                    af = alt_count / total if total > 0 else None
                except (ValueError, IndexError):
                    ref_count = alt_count = af = None

                record[f"AD_REF_{i+1}"] = ref_count
                record[f"AD_ALT_{i+1}"] = alt_count
                record[f"AF_{i+1}"] = af

            all_records.append(record)

    return pd.DataFrame(all_records)


# Main 

treatment = "pano"
cell_line = "MLS-Avory"

# Input file paths

treatment_files = [
    f"MLS-1765/MLS-1765_{treatment}_all3_r3_3utr.txt",
    f"MLS-402/MLS-402_{treatment}_all3_r3_3utr.txt",
    f"MLS-Avory/MLS-Avory_{treatment}_all3_r3_3utr.txt"
]

dmso_files = [
    "MLS-1765/MLS-1765_all3_r3_3utr.txt",
    "MLS-402_all3_r3_3utr.txt",
    "MLS-Avory/MLS-Avory_all3_r3_3utr.txt"
]


# Generate maps and get unique positions
treat_maps = get_position_to_gene_map(*treatment_files)
dmso_maps = get_position_to_gene_map(*dmso_files)

merged_treatment_map = {k: v for d in treat_maps for k, v in d.items()}
merged_dmso_map = {k: v for d in dmso_maps for k, v in d.items()}

common_treat_positions = get_common_positions(*treatment_files)
common_dmso_positions = get_common_positions(*dmso_files)

unique_treatment_positions = plot_venn_and_get_treatment_specific_positions(
    common_treat_positions,
    common_dmso_positions,
    treatment
)

# Format VCF positions
vcf_positions = [f"chr{pos}" if not pos.startswith("chr") else pos for pos in unique_treatment_positions]
gene_list = {f"chr{pos}" if not pos.startswith("chr") else pos: merged_treatment_map[pos] for pos in unique_treatment_positions}

# VCF file path
treatment_vcf_file = f"{cell_line}/{cell_line}_{treatment}_RNA_AG_TC_DP5_ALL3.vcf" 
no_treatment_vcf_file = f"{cell_line}/{cell_line}_RNA_AG_TC_DP5_ALL3.vcf"

# Extract AD/AF
treatment_df = extract_ad_af_from_vcf(treatment_vcf_file, vcf_positions, gene_list)
no_treatment_df = extract_ad_af_from_vcf(no_treatment_vcf_file, vcf_positions, gene_list)
print(no_treatment_df)
# Merge and compute AF ratios
merged_df = pd.merge(treatment_df, no_treatment_df, on=["Position", "GeneName"],how="left", suffixes=("_T", "_NT"))

# Calculate max AF in treatment and no-treatment samples
treatment_af_cols = [col for col in merged_df.columns if col.startswith("AF_") and col.endswith("_T")]
no_treatment_af_cols = [col for col in merged_df.columns if col.startswith("AF_") and col.endswith("_NT")]

merged_df["Max_AF_T"] = merged_df[treatment_af_cols].max(axis=1, skipna=True)
merged_df["Max_AF_NT"] = merged_df[no_treatment_af_cols].max(axis=1, skipna=True)

# Compute max AF ratio
merged_df["Max_AF_Ratio"] = merged_df.apply(
    lambda row: row["Max_AF_T"] / row["Max_AF_NT"]
    if pd.notnull(row["Max_AF_T"]) and pd.notnull(row["Max_AF_NT"]) and row["Max_AF_NT"] != 0
    else None,
    axis=1
)

# Save output
# Replace NaN with 0 (can also use .fillna(0, inplace=True) if you prefer)
output_path = f"250703_{treatment}_{cell_line}_AF_Ratios_AD_3replicates_3utr.csv"
merged_df_filled = merged_df.fillna(0)

# Save to CSV with formatted floats
merged_df_filled.to_csv(
    output_path,
    index=False,
    float_format="%.6f"
)
#output_path = f"250626_{treatment}_{cell_line}_AF_Ratios_AD.csv"
merged_df_filled.to_csv(output_path, index=False, sep=";")

print("\nAF Ratios (Treatment / No-Treatment):")
print(merged_df_filled)

# Save to CSV if needed
#merged_af_df.to_csv(f"250613_{treatment}_{cell_line}_AF_Ratios_AD.csv", index=False)
