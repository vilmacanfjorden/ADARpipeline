#!/bin/bash -l
#$ -cwd
#$ -N runVC_MLS-Avory_azd
#$ -o $JOB_NAME_$TASK_ID.stdout
#$ -e $JOB_NAME_$TASK_ID.stderr
#$ -q [queue name]
#$ -pe mpi 1 

module load bcftools
module load htslib

# Define paths
DNA_VCF="path/to/DNAfile"
# DBSNP_VCF="dbsnp_150.vcf.gz"  # Uncomment and update path if SNP filtering desired
OUTDIR="outputdir"
TMP_ISEC_DIR="${OUTDIR}/tmp_isec_RNA_AZD_MLS-Avory"


# Treated RNA VCFs (AZD) called by Mutect2
RNA_VCFs=(
  "replicate1"
  "replicate2"
  "replicate3"
  "replicate4"
)

# Untreated RNA VCFs (DMSO) called by Mutect2
UNTREATED_RNA_VCFs=(
  "replicate1"
  "replicate1"
  "replicate3"
  "replicate4"
)

# Intersect DNA and treated RNA to get RNA-only variants
mkdir -p "$TMP_ISEC_DIR"

for RNA_VCF in "${RNA_VCFs[@]}"; do
  samplename=$(basename "$RNA_VCF" .mutect.vcf.gz)
  outdir="${TMP_ISEC_DIR}/${samplename}"
  mkdir -p "$outdir"
  bcftools isec -p "$outdir" -n=1 "$RNA_VCF" "$DNA_VCF"
  bgzip -f "$outdir/0000.vcf"
  bcftools index -f "$outdir/0000.vcf.gz"
done

# Merge RNA-only VCFs from treated replicates
INTERSECTED_VCFs=()
for RNA_VCF in "${RNA_VCFs[@]}"; do
  samplename=$(basename "$RNA_VCF" .mutect.vcf.gz)
  INTERSECTED_VCFs+=("${TMP_ISEC_DIR}/${samplename}/0000.vcf.gz")
done

MERGED_VCF="${OUTDIR}/MLS-Avory_azd_merged_RNA_updated_r3.vcf.gz"
bcftools merge --force-samples "${INTERSECTED_VCFs[@]}" -Oz -o "$MERGED_VCF"
bcftools index -f "$MERGED_VCF"

# Keep only A>G and T>C changes (ADAR candidates)
AGTC_VCF="${OUTDIR}/MLS-Avory_azd_RNA_AG_TC_updated_r3.vcf.gz"
bcftools view -i '(REF="A" & ALT="G") | (REF="T" & ALT="C")' "$MERGED_VCF" -Oz -o "$AGTC_VCF"
bcftools index -f "$AGTC_VCF"

# Filter for consistent presence and depth > 5 in all at least 3 replicates
FILTERED_ALL_VCF="${OUTDIR}/MLS-Avory_azd_RNA_AG_TC_DP5_ALL3.vcf.gz"
bcftools view -i 'COUNT(GT!="./." & FORMAT/DP>5)>=3' "$AGTC_VCF" -Oz -o "$FILTERED_ALL_VCF"
bcftools index -f "$FILTERED_ALL_VCF"

# Process untreated DMSO replicates to exclude variants found in untreated samples (can be to conservative)
TMP_ISEC_DMSO_DIR="${OUTDIR}/tmp_isec_RNA_DMSO_MLS-Avory"
mkdir -p "$TMP_ISEC_DMSO_DIR"

UNTREATED_INTERSECTED_VCFs=()
for DMSO_VCF in "${UNTREATED_RNA_VCFs[@]}"; do
  samplename=$(basename "$DMSO_VCF" .mutect.vcf.gz)
  outdir="${TMP_ISEC_DMSO_DIR}/${samplename}"
  mkdir -p "$outdir"
  bcftools isec -p "$outdir" -n=1 "$DMSO_VCF" "$DNA_VCF"
  bgzip -f "${outdir}/0000.vcf"
  bcftools index -f "${outdir}/0000.vcf.gz"
  UNTREATED_INTERSECTED_VCFs+=("${outdir}/0000.vcf.gz")
done

UNTREATED_MERGED_VCF="${OUTDIR}/MLS-Avory_dmso_merged_RNA.vcf.gz"
bcftools merge --force-samples "${UNTREATED_INTERSECTED_VCFs[@]}" -Oz -o "$UNTREATED_MERGED_VCF"
bcftools index -f "$UNTREATED_MERGED_VCF"

# Subtract untreated variants from treated filtered variants
TREATMENT_SPECIFIC_VCF="${OUTDIR}/MLS-Avory_azd_DP5_ALL3_noDMSO.vcf.gz"
#echo "Overlap count:"
#bcftools isec -n=2 -w1 -Oz -p tmp_overlap "$FILTERED_ALL_VCF" "$UNTREATED_MERGED_VCF"
#bcftools view -H tmp_overlap/0000.vcf.gz | wc -l

bcftools isec -C -p MLS-Avory_tmp_isec_out_azd_nodmso_r3 "$FILTERED_ALL_VCF" "$UNTREATED_MERGED_VCF" -Oz -o "$TREATMENT_SPECIFIC_VCF" # Can be to conservative
#bcftools index -f "$TREATMENT_SPECIFIC_VCF"
mv MLS-Avory_tmp_isec_out_azd_nodmso_r3/0000.vcf.gz MLS-Avory_azd_nodmso_treatment_specific_variants_r3.vcf.gz
bcftools index -f MLS-Avory_azd_nodmso_treatment_specific_variants_r3.vcf.gz

# Step 6: Optional — Remove known common SNPs (uncomment and update DBSNP_VCF path if available)
# FINAL_VCF="${OUTDIR}/MLS-Avory_azd_DP5_ALL4_noDMSO_noSNPs.vcf.gz"
# bcftools isec -C "$TREATMENT_SPECIFIC_VCF" "$DBSNP_VCF" -Oz -o "$FINAL_VCF"
# bcftools index -f "$FINAL_VCF"

# Final output message
echo "Final filtered treatment-specific RNA editing candidates (DMSO removed) saved in:"
echo "$TREATMENT_SPECIFIC_VCF"

