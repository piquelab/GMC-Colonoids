#!/bin/bash
#SBATCH -q express -p erprp
#SBATCH --job-name=genotype_pileup_ind
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o logs/pileup_%A_%a.out
#SBATCH -e logs/pileup_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=go7535@wayne.edu

set -euo pipefail

WORKDIR="/rs/rs_grp_gmccol/GMCcol_ASE_analysis/per_individual_joint_pileup/keep_zero_ALT_sites_analysis/pileups"
ASEDIR="/rs/rs_grp_gmccol/GMCcol_ASE_analysis/per_individual_joint_pileup/keep_zero_ALT_sites_analysis/ASE_input"
BAMDIR="/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/bams"
fastaRef="/nfs/rprdata/Anthony/data/HISAT2_Index/hg38_noChr.fa"
dbSnpFile="/nfs/rprdata/RefSnps/dbSNP144.hg38/snp144Common.bed.gz"

ncores="${SLURM_CPUS_PER_TASK:-8}"

cd "$WORKDIR"

module load bcftools/1.19
module load htslib
module load R

IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" individuals.txt)
echo "[$(date)] Individual: $IND"


# 1) Build filtered genotype VCF and bed file of credible het sites

bash get_het_sites.sh "$IND"

GENO_VCF="hets/${IND}.filtered_genotype_hets.nochr.DNA.vcf.gz"
HET_BED="hets/${IND}.filtered_genotype_hets.sites.bed.gz"

[[ -f "$HET_BED" ]] || { echo "ERROR: missing het BED: $HET_BED" >&2; exit 1; }

# 2) Build list of bam files from same individual/cell line

BAM_LIST="bam_files_${IND}.txt"

awk -v ind="$IND" -v bamdir="$BAMDIR" 'BEGIN{FS=OFS="\t"} $2==ind {print bamdir "/" $1 "_clean.bam"}' \
  sample_individual_map.txt > "$BAM_LIST"

awk 'NF' "$BAM_LIST" | while read -r bam; do
  [[ -f "$bam" ]] && echo "$bam"
done > "${BAM_LIST}.tmp"
mv "${BAM_LIST}.tmp" "$BAM_LIST"

echo "[$(date)] BAMs found for $IND: $(wc -l < "$BAM_LIST")"
if [[ $(wc -l < "$BAM_LIST") -eq 0 ]]; then
  echo "ERROR: No BAMs found for $IND after existence check" >&2
  exit 1
fi

# 3) Generate pileup VCF across all samples in individual -> multi sample pileup
VCF_RAW="${IND}.pileup.vcf.gz"

if [[ ! -f "$VCF_RAW" ]]; then
  echo "[$(date)] Generating pileup: $VCF_RAW"

  bcftools mpileup -d 100000 --min-MQ 30 --min-BQ 30 \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP \
    -f "$fastaRef" -T "$HET_BED" \
    --threads "$ncores" \
    --bam-list "$BAM_LIST" \
  | bcftools call -m -Oz --threads "$ncores" -o "$VCF_RAW"

  tabix -f -p vcf "$VCF_RAW"
fi

# 4️) Merge genotype + pileup to annotate ALT allele information
VCF_MERGED="${IND}.genotype_plus_pileup.merged.vcf.gz"
if [[ ! -f "$VCF_MERGED" ]]; then
  echo "[$(date)] Merging pileup with genotype information: $VCF_MERGED"
  bcftools merge \
    -m none \
    -Oz \
    -o "$VCF_MERGED" \
    "$GENO_VCF" \
    "$VCF_RAW"

  tabix -f -p vcf "$VCF_MERGED"
fi

# 5️) Keep only biallelic SNPs
VCF_CLEAN="${IND}.merged.clean.snps.vcf.gz"
if [[ ! -f "$VCF_CLEAN" ]]; then
  echo "[$(date)] Filtering for biallelic SNPs: $VCF_CLEAN"
  bcftools view \
    -m2 -M2 \
    -v snps \
    -Oz -o "$VCF_CLEAN" \
    "$VCF_MERGED"

  tabix -f -p vcf "$VCF_CLEAN"
fi

# 6) Annotate rsIDs from dbSNP
VCF_RSID="${IND}.merged.clean.snps.rsid.vcf.gz"
if [[ ! -f "$VCF_RSID" ]]; then
  echo "[$(date)] Annotating rsID from dbSNP: $VCF_RSID"
  bcftools annotate \
    -a "$dbSnpFile" \
    -c CHROM,FROM,TO,ID \
    -Oz -o "$VCF_RSID" \
    "$VCF_CLEAN"

  tabix -f -p vcf "$VCF_RSID"
fi

# 7️) Export RNA Allele Depth (AD) matrix
 echo "[$(date)] Making AD matrix for $IND"
 
 AD_WIDE="${IND}.AD.wide.tsv.gz"
if [[ ! -f "$AD_WIDE" ]]; then
    bcftools query -l "$VCF_RAW" > "${IND}.samples.full.txt"
    
    awk -F/ '{print $NF}' "${IND}.samples.full.txt" \
      | sed 's/_clean\.bam$//' \
      > "${IND}.samples.txt"
    
    RNA_SAMPLES_CSV=$(paste -sd, "${IND}.samples.full.txt")
    
    
    
    bcftools query \
      -s "$RNA_SAMPLES_CSV" \
      -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%AD]\n' \
      "$VCF_RSID" \
    | bgzip > "$AD_WIDE"
fi

# 8) Allelic Imbalance pre-processing - make long table of row per SNP/sample
echo "[$(date)] Running AI processing $IND"
 
Rscript ai_processing.R \
  "$IND" \
  "$WORKDIR/sample_individual_map.txt" \
  "$WORKDIR" \
  "$ASEDIR" 

echo "[$(date)] Finished $IND"

# 9) Pileup QC

echo "Raw pileup rows:"
bcftools view -H "$VCF_RAW" | wc -l

echo "Dot ALT count in raw pileup:"
bcftools view -H "$VCF_RAW" | awk '$5=="."' | wc -l

echo "Merged VCF rows:"
bcftools view -H "$VCF_MERGED" | wc -l

echo "Dot ALT count in merged VCF:"
bcftools view -H "$VCF_MERGED" | awk '$5=="."' | wc -l


echo "Merged biallelic SNP VCF rows:"
bcftools view -H "$VCF_CLEAN" | wc -l

echo "Dot ALT count in merged biallelic SNP VCF:"
bcftools view -H "$VCF_CLEAN" | awk '$5=="."' | wc -l
