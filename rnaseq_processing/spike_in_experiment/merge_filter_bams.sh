#!/bin/bash
#SBATCH -q express -p erprp
#SBATCH --mem=100G
#SBATCH --time=4-01:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=go7535@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err

module load samtools
if [[ ! -f ${var}_merged_clean.bam ]]; then

	###Merge Reads###
	samtools merge ${var}_merged.bam ${var}_run1_sorted.bam ${var}_run2_sorted.bam 
	samtools index ${var}_merged.bam
	samtools view -c ${var}_merged.bam > ${var}_merged_sorted_count.txt

	###Quality Filter###
	samtools view -b1 -q10 ${var}_merged.bam > ${var}_merged_quality.bam
	samtools index ${var}_merged_quality.bam
	samtools view -c ${var}_merged_quality.bam > ${var}_merged_quality_count.txt


	###Deduplication###
	samtools rmdup ${var}_merged_quality.bam ${var}_merged_clean.bam
	samtools index ${var}_merged_clean.bam
	samtools view -c ${var}_merged_clean.bam > ${var}_merged_clean_count.txt
fi