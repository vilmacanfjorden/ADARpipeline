#!/bin/bash
#$ -cwd
#$ -N runVC_HT1080_FUS
#$ -o $JOB_NAME_$TASK_ID.stdout
#$ -e $JOB_NAME_$TASK_ID.stderr
#$ -q [queue name]
#$ -pe mpi 1 

# Load modules
module load bcftools
module load htslib


DNA_VCF="path/to/DNAfile"
OUTDIR="outputdir"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "$OUTDIR" "$TMPDIR"

# FUS RNA replicates (4)
FUS_RNA_VCFs=(
  "replicate1"
  "replicate2",
  "replicate3",
  "replicate4",
)

# WT RNA replicates (3)
WT_RNA_VCFs=(
  "replicate1",
  "replicate1",
  "replicate3"
)
# Intersect RNA and DNA to get RNA-specific variants 
INTERSECTED_FUS_VCFS=()
INTERSECTED_WT_VCFS=()

# FUS
for RNA_VCF in "${FUS_RNA_VCFs[@]}"; do
  sample=$(basename "$RNA_VCF" .vcf.gz)
  outdir="${TMPDIR}/isec_${sample}"
  mkdir -p "$outdir"
  bcftools isec -p "$outdir" -n=1 "$RNA_VCF" "$DNA_VCF"
  bgzip -f "$outdir/0000.vcf"
  bcftools index -f "$outdir/0000.vcf.gz"
  INTERSECTED_FUS_VCFS+=("${outdir}/0000.vcf.gz")
done

# WT
for RNA_VCF in "${WT_RNA_VCFs[@]}"; do
  sample=$(basename "$RNA_VCF" .vcf.gz)
  outdir="${TMPDIR}/isec_${sample}"
  mkdir -p "$outdir"
  bcftools isec -p "$outdir" -n=1 "$RNA_VCF" "$DNA_VCF"
  bgzip -f "$outdir/0000.vcf"
  bcftools index -f "$outdir/0000.vcf.gz"
  INTERSECTED_WT_VCFS+=("${outdir}/0000.vcf.gz")
done

# Merge VCFs per group
MERGED_FUS_VCF="${OUTDIR}/FUS_merged.vcf.gz"
MERGED_WT_VCF="${OUTDIR}/WT_merged.vcf.gz"

bcftools merge --force-samples "${INTERSECTED_FUS_VCFS[@]}" -Oz -o "$MERGED_FUS_VCF"
bcftools index -f "$MERGED_FUS_VCF"

bcftools merge --force-samples "${INTERSECTED_WT_VCFS[@]}" -Oz -o "$MERGED_WT_VCF"
bcftools index -f "$MERGED_WT_VCF"

# Filter A>G and T>C variants 
FUS_AGTC="${OUTDIR}/FUS_AG_TC.vcf.gz"
WT_AGTC="${OUTDIR}/WT_AG_TC.vcf.gz"

bcftools view -i '(REF="A" & ALT="G") | (REF="T" & ALT="C")' "$MERGED_FUS_VCF" -Oz -o "$FUS_AGTC"
bcftools index -f "$FUS_AGTC"

bcftools view -i '(REF="A" & ALT="G") | (REF="T" & ALT="C")' "$MERGED_WT_VCF" -Oz -o "$WT_AGTC"
bcftools index -f "$WT_AGTC"

# Filter for DP > 5 and present in at least 3 replicates 
FUS_FILTERED="${OUTDIR}/FUS_AG_TC_DP5_in3.vcf.gz"
WT_FILTERED="${OUTDIR}/WT_AG_TC_DP5_in3.vcf.gz"

bcftools view -i 'COUNT(GT!="./." & FORMAT/DP>5)>=3' "$FUS_AGTC" -Oz -o "$FUS_FILTERED"
bcftools index -f "$FUS_FILTERED"

bcftools view -i 'COUNT(GT!="./." & FORMAT/DP>5)>=3' "$WT_AGTC" -Oz -o "$WT_FILTERED"
bcftools index -f "$WT_FILTERED"

# Subtract WT from FUS to get FUS-specific variants
bcftools isec -C -p "${OUTDIR}/isec_FUS_only" "$FUS_FILTERED" "$WT_FILTERED"
bgzip -f "${OUTDIR}/isec_FUS_only/0000.vcf"
bcftools index -f "${OUTDIR}/isec_FUS_only/0000.vcf.gz"

FINAL_OUTPUT="${OUTDIR}/FUS_specific_variants.vcf.gz"
mv "${OUTDIR}/isec_FUS_only/0000.vcf.gz" "$FINAL_OUTPUT"
bcftools index -f "$FINAL_OUTPUT"

# Done 
echo "Pipeline completed successfully!"
echo "Final FUS-specific variants (A>G, T>C, DP>5, in ≥3 replicates): $FINAL_OUTPUT"
