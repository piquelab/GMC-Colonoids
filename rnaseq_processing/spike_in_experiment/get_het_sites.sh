#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=100G
#SBATCH --time=4-00:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o logs/het_%j.out
#SBATCH -e logs/het_%j.err

set -euo pipefail

ind=$1
WORKDIR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/pileups"
cd "$WORKDIR"

vcfFileMain="/wsu/home/groups/piquelab/OurData/genewiz/30-792835734/MB_imputation_redo/locoglimpse2/ref.ac1.vcf.gz"
dbSnpFile="/nfs/rprdata/RefSnps/dbSNP144.hg38/snp144Common.bed.gz"
GP_THRESH="0.999"

module load bcftools/1.19
module load htslib


TMP_VCF_CHR="hets/${ind}.tmp.sample.chr.vcf.gz"
TMP_VCF_NOCHR="hets/${ind}.tmp.sample.nochr.vcf.gz"
GENO_VCF="hets/${ind}.filtered_genotype_hets.nochr.vcf.gz"
GENO_VCF_RENAMED="hets/${ind}.filtered_genotype_hets.nochr.DNA.vcf.gz"
HET_BED="hets/${ind}.filtered_genotype_hets.sites.bed.gz"

STRIP_CHR_MAP="hets/strip_chr_map.txt"
RENAME_SAMPLES="hets/${ind}.rename_samples.txt"



# 1) subset to one individual

bcftools view \
  -s "$ind" \
  -Oz -o "$TMP_VCF_CHR" \
  "$vcfFileMain"

tabix -f -p vcf "$TMP_VCF_CHR"


# 2) strip chr names to match no-chr pileup/dbSNP

bcftools annotate \
  --rename-chrs "$STRIP_CHR_MAP" \
  -Oz -o "$TMP_VCF_NOCHR" \
  "$TMP_VCF_CHR"

tabix -f -p vcf "$TMP_VCF_NOCHR"


# 3) filter: GP het threshold, biallelic SNPs, dbSNP-common

bcftools view \
  -i "FMT/GP[0:1] > ${GP_THRESH}" \
  -m2 -M2 \
  -v snps \
  -T "$dbSnpFile" \
  -Oz -o "$GENO_VCF" \
  "$TMP_VCF_NOCHR"

tabix -f -p vcf "$GENO_VCF"


# 4) rename genotype sample for merge
printf "%s\t%s_DNA\n" "$ind" "$ind" > "$RENAME_SAMPLES"

bcftools reheader \
  -s "$RENAME_SAMPLES" \
  -o "$GENO_VCF_RENAMED" \
  "$GENO_VCF"

tabix -f -p vcf "$GENO_VCF_RENAMED"


# 5) make BED for pileup targeting from final genotype VCF
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' "$GENO_VCF" \
  | bgzip > "$HET_BED"

tabix -f -p bed "$HET_BED"

echo "Wrote: $GENO_VCF_RENAMED"
echo "Wrote: $HET_BED"
echo "Number of targeted SNPs: $(bcftools view -H "$GENO_VCF" | wc -l)"