#!/bin/bash -l
#$ -cwd
#$ -o run_STAR.stdout
#$ -e run_STAR.stderr
#$ -q [queue name]
#$ -pe mpi 40


module load star/2.5.2b
module load samtools

# paths
STAR_index_path=/path/to/STARindex
GTF_file=/path/to/gtf/gencode.v34.annotation.gtf
Reference_fasta=path/to/hg38/hg38.fa
fastq_path=MLS/
output_path=STARoutput/

# List of directories (sample names) with fastq files 
samplelist=(
)

mkdir STARoutput

for i in "${samplelist[@]}"
do
     STAR --genomeDir $STAR_index_path --runThreadN 40 --readFilesIn  $fastq_path/$i"_R1_001.fastq.gz" $fastq_path/$i"_R2_001.fastq.gz" --readFilesCommand zcat --outFileNamePrefix $output_path/$i"_" --outSAMtype BAM SortedByCoordinate
    samtools index $output_path/$i"_Aligned.out.bam"
done

