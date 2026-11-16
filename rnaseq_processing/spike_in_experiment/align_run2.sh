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


module load hisat2/2.0.4
module load samtools/1.11

echo "Processing sample: ${var}"
# Skip if final output already exists
if [ -f "${var}_run2_clean.bam" ]; then
    echo "${var}: clean BAM already exists. Skipping."
    exit 0
fi


genomeindex=/wsu/el7/groups/piquelab/refData/hisat2/grch38_snp_tran/genome_snp_tran
###Align Reads###-
hisat2 -p 8 -x ${genomeindex} -1 /rs/rs_grp_ibdeqtl/GxM_Validation/seq_data/run2_trimmed_fastqs/FL-GG-1s-FL012-FC02-${var}_L001_R1_001_val_1.fq.gz \
                              -2 /rs/rs_grp_ibdeqtl/GxM_Validation/seq_data/run2_trimmed_fastqs/FL-GG-1s-FL012-FC02-${var}_L001_R2_001_val_2.fq.gz \
      2> ${var}_run2_aligned.bam.e | samtools view -b1 - > ${var}_run2_aligned.bam



###Sort Reads###
samtools sort -@ 4 -T tmp_${var}_run2_aligned.bam -o ${var}_run2_sorted.bam ${var}_run2_aligned.bam
samtools index ${var}_run2_sorted.bam
samtools view -c ${var}_run2_sorted.bam > ${var}_run2_sorted_count.txt


###Quality Filter###
samtools view -b1 -q10 ${var}_run2_sorted.bam > ${var}_run2_quality.bam
samtools index ${var}_run2_quality.bam
samtools view -c ${var}_run2_quality.bam > ${var}_run2_quality_count.txt


###Deduplication###
samtools rmdup ${var}_run2_quality.bam ${var}_run2_clean.bam
samtools index ${var}_run2_clean.bam
samtools view -c ${var}_run2_clean.bam > ${var}_run2_clean_count.txt
