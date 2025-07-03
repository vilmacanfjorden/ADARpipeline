#!/bin/bash -l
#$ -cwd
#$ -N runVC_MLS-402
#$ -o $JOB_NAME_$TASK_ID.stdout
#$ -e $JOB_NAME_$TASK_ID.stderr
#$ -q [queue name]
#$ -pe mpi 1 

module load bcftools
module load htslib

# Input VCF files
DNA_VCF="DNAfile"
RNA_VCFs=(
  "replicate1"
  "replicate2"
  "replicate3"
  "replicate4"
)

# Create directories for intersections
mkdir -p isec_RNA_MLS-402

# Intersect DNA VCF with each RNA VCF (keep RNA-only variants)
for RNA_VCF in "${RNA_VCFs[@]}"; do
  samplename=$(basename "$RNA_VCF" .mutect.vcf.gz)
  bcftools isec -p "isec_RNA_MLS-402/$samplename" -n=1 "$RNA_VCF" "$DNA_VCF"
done

# Compress and index the intersected VCFs
for dir in isec_RNA_MLS-402/*/; do
  bgzip -f "${dir}0000.vcf"
  bcftools index -f "${dir}0000.vcf.gz"
done

# Define an array for the intersected VCFs
INTERSECTED_RNA_VCFs=(
  "isec_RNA_MLS-402/MLS-402-R1/0000.vcf.gz"
  "isec_RNA_MLS-402/MLS-402-R2/0000.vcf.gz"
  "isec_RNA_MLS-402/MLS-402-R3/0000.vcf.gz"
  "isec_RNA_MLS-402/MLS-402-R4/0000.vcf.gz"
)

# Merge RNA replicates into a multisample VCF
bcftools merge --force-samples "${INTERSECTED_RNA_VCFs[@]}" -Oz -o MLS-402_merged_RNA_r3.vcf.gz
bcftools index -f MLS-402_merged_RNA_r3.vcf.gz

# Filter variants with depth > 5 in ALL 4 replicates and genotype present
bcftools view -i 'COUNT(GT!="." & FORMAT/DP>5)>=3' MLS-402_merged_RNA_r3.vcf.gz -Oz -o MLS-402_RNA_DP5_ALL3.vcf.gz
bcftools index -f MLS-402_RNA_DP5_ALL3.vcf.gz

# Filter for A->G and T->C variants (ADAR editing candidates)
bcftools view -i '(REF="A" & ALT="G") | (REF="T" & ALT="C")' MLS-402_RNA_DP5_ALL3.vcf.gz -Oz -o MLS-402_RNA_AG_TC_DP5_ALL3.vcf.gz
bcftools index -f MLS-402_RNA_AG_TC_DP5_ALL3.vcf.gz

echo "Filtered RNA variants (present with DP>5 in all replicates and A>G/T>C) saved in MLS-402_RNA_AG_TC_DP5_ALL3.vcf.gz"

