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

WORKDIR="/rs/rs_grp_gmccol/GMCcol_ASE_analysis/per_individual_joint_pileup/keep_zero_ALT_sites_analysis/stringent_SNP_filter_analysis/ASE/QuASAR_results"
cd "$WORKDIR"

QUASAR_GZ="all_individuals.Quasar_output.ASE.txt.gz"
ANNO="/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.no_head.bam_sort.tsv.gz"

# ------------------------------
# 1) Get SNPs tested for ASE
#    Create ASE_SNPs.txt with 6 columns (BED-like):
#    chr  pos0  pos1  ref  alt  rsID
# ------------------------------

zcat "$QUASAR_GZ" \
  | awk 'BEGIN{FS=OFS="\t"} NR>1 {print $2, $3-1, $3, $5, $6, $4}' \
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
