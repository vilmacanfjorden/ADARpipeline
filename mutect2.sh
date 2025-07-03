#!/bin/bash
#$ -cwd
#$ -N mutect_MLS
#$ -o $JOB_NAME_$TASK_ID.stdout
#$ -e $JOB_NAME_$TASK_ID.stderr
#$ -q [queue name]
#$ -pe mpi 20

# Load required modules
module load samtools
module load picard-tools/2.9.0
module load java
module load gatk

# Paths to resources
Reference_fasta="/path/to/hg38.fa"
realignment_file1="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# Input and output directories
BAM_DIR="STARoutput/"
OUTPUT_DIR="STARoutput/"
PT_HOME="/path/to/picard/2.9.0/"

# Loop through all BAM files in the input directory
for BAM_FILE in "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam; do
  if [[ -f "$BAM_FILE" ]]; then
    echo "Processing $BAM_FILE"
    
    # Extract sample name
    samplename=$(basename "$BAM_FILE" | sed 's/_Aligned.sortedByCoord.out.bam//g')

    # MarkDuplicates
    java -Xmx10g -jar $PT_HOME/picard.jar MarkDuplicates \
      INPUT="$BAM_FILE" \
      METRICS_FILE="$OUTPUT_DIR/${samplename}_Aligned.out.bam.metrics" \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=LENIENT \
      CREATE_INDEX=true \
      OUTPUT="$OUTPUT_DIR/${samplename}_md.bam"

    # SplitNCigarReads
    /path/to/gatk/4.6.1.0/gatk SplitNCigarReads \
      -R "$Reference_fasta" \
      -I "$OUTPUT_DIR/${samplename}_md.bam" \
      -O "$OUTPUT_DIR/${samplename}.bam"

    # AddOrReplaceReadGroups
    java -Xmx4g -jar $PT_HOME/picard.jar AddOrReplaceReadGroups \
      I="$OUTPUT_DIR/${samplename}.bam" \
      O="$OUTPUT_DIR/${samplename}.rg.bam" \
      RGID=4 \
      RGLB="$samplename" \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=20

    # BaseRecalibrator
    /path/to/gatk/4.6.1.0/gatk BaseRecalibrator \
      -R "$Reference_fasta" \
      -I "$OUTPUT_DIR/${samplename}.rg.bam" \
      --known-sites "$realignment_file1" \
      -O "$OUTPUT_DIR/${samplename}.recal_data.table"

    # ApplyBQSR
    /path/to/gatk/4.6.1.0/gatk ApplyBQSR \
      -R "$Reference_fasta" \
      -I "$OUTPUT_DIR/${samplename}.rg.bam" \
      --bqsr-recal-file "$OUTPUT_DIR/${samplename}.recal_data.table" \
      -O "$OUTPUT_DIR/${samplename}_recal.bam"

     # Optional: Run Mutect2 for variant calling (comment if not needed)
     /path/to/gatk/4.6.1.0/gatk --java-options "-Xmx4g" Mutect2 \
       -R "$Reference_fasta" \
       -I "$OUTPUT_DIR/${samplename}_recal.bam" \
       --tumor-sample "$samplename" \
       -O "$OUTPUT_DIR/${samplename}.mutect.vcf.gz"

  else
    echo "No BAM files found in $BAM_DIR"
    break
  fi
done
