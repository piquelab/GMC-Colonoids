#!/bin/bash
#SBATCH -q express -p erprp
#SBATCH --mem=100G
#SBATCH --time=4-00:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o slurm_%j.out
#SBATCH -e slurm_%j.err


# Purpose: add gene info to SNPs tested for ASE
# Made by AF  8/1/19
# Shreya Nirmalan - modified for GMCcol project 1/30/2026
# Update: build ASE_SNPs.txt from QuASAR output (gz)

set -u
set -v

WORKDIR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/ASE"
cd "$WORKDIR"

PILEUP_GZ="Colon-87-C_allele_counts_long.txt.gz"
ANNO="/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.no_head.bam_sort.tsv.gz"

# ------------------------------
# 1) Get SNPs tested for ASE
#    Create ASE_SNPs.txt with 6 columns (BED-like):
#    chr  pos0  pos1  ref  alt  rsID
# ------------------------------

zcat "$PILEUP_GZ" \
  | awk 'BEGIN{FS=","; OFS="\t"} NR>1 {print $1, $2-1, $2, $4, $5, $3}' \
  | sort -k1,1 -k2,2n \
  | uniq \
  > ASE_SNPs.txt

wc -l ASE_SNPs.txt

# ------------------------------
# 2) Intersect with gene annotations
# ------------------------------
module load bedtools

bedtools intersect -sorted -u \
  -a ASE_SNPs.txt \
  -b "$ANNO" | wc -l

bedtools intersect -sorted -wa -wb \
  -a ASE_SNPs.txt \
  -b "$ANNO" \
  > ASE_SNPs.genes.txt

cut -f 1-6,14,16,17 ASE_SNPs.genes.txt \
  | sort -k1,1 -k2,2n | uniq \
  > ASE_SNPs.genes.uniq.txt