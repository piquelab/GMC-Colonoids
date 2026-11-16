#!/bin/bash
#SBATCH -q express -p erprp
#SBATCH --job-name=GxM_spikein_pileup_ind
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -o logs/GxM_spikein_%A_%a.out
#SBATCH -e logs/GxM_spikein_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=go7535@wayne.edu

set -euo pipefail
module unload gnu7
module load gnu9
module load bcftools/1.19
module load htslib
module load samtools
module load R

IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" individuals.txt)

WORKDIR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/pileups"
ASEDIR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/ASE"
BAMDIR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/bams"
fastaRef="/nfs/rprdata/Anthony/data/HISAT2_Index/hg38_noChr.fa"
dbSnpFile="/nfs/rprdata/RefSnps/dbSNP144.hg38/snp144Common.bed.gz"
RSCRIPT="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/pileups/ai_processing.R"

ncores="${SLURM_CPUS_PER_TASK:-8}"
MAX_DEPTH=100000

cd "$WORKDIR"

echo "[$(date)] Individual: $IND"

# 1) Build filtered genotype VCF and bed file of credible het sites

bash get_het_sites.sh "$IND"

GENO_VCF="hets/${IND}.filtered_genotype_hets.nochr.DNA.vcf.gz"
HET_BED="hets/${IND}.filtered_genotype_hets.sites.bed.gz"

[[ -f "$GENO_VCF" ]] || { echo "ERROR: missing genotype VCF: $GENO_VCF" >&2; exit 1; }
[[ -f "$HET_BED" ]]  || { echo "ERROR: missing het BED: $HET_BED" >&2; exit 1; }

# 2) Build list of bam files from same individual/cell line

# Build BAM list for individual
BAM_LIST="bam_files_${IND}.txt"
BAD_LOG="logs/${IND}.bad_bams.txt"

awk -v ind="$IND" -v bamdir="$BAMDIR" 'BEGIN{FS=OFS="\t"} $2==ind {print bamdir "/" $1 "_merged_clean.bam"}' \
  sample_individual_map.txt > "${BAM_LIST}.tmp"

: > "$BAM_LIST"
: > "$BAD_LOG"

# Keep only existing + valid BAMs, log bad/missing ones
while read -r bam; do
  [[ -z "$bam" ]] && continue

  if [[ -f "$bam" ]] && samtools quickcheck "$bam" 2>/dev/null; then
    echo "$bam" >> "$BAM_LIST"
  else
    echo "$bam" >> "$BAD_LOG"
  fi
done < "${BAM_LIST}.tmp"

rm -f "${BAM_LIST}.tmp"

echo "[$(date)] Good BAMs for $IND: $(wc -l < "$BAM_LIST")"
echo "[$(date)] Skipped BAMs logged to: $BAD_LOG"

if [[ ! -s "$BAM_LIST" ]]; then
  echo "ERROR: No usable BAMs found for $IND" >&2
  exit 1
fi

# 3) Generate pileup VCF across all samples in individual -> multi sample pileup
VCF_RAW="${IND}.pileup.vcf.gz"

if [[ ! -f "$VCF_RAW" ]]; then
  echo "[$(date)] Generating pileup: $VCF_RAW"

  bcftools mpileup -Ou -d "$MAX_DEPTH" --min-MQ 30 --min-BQ 30 \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP \
    -f "$fastaRef" -T "$HET_BED" \
    --threads "$ncores" \
    --bam-list "$BAM_LIST" \
  | bcftools call -m -Ou --threads "$ncores" \
  | bcftools view -Oz --threads "$ncores" -o "$VCF_RAW"

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
 
Rscript "$RSCRIPT" \
  "$IND" \
  "$WORKDIR" \
  "$ASEDIR" 

echo "[$(date)] Finished $IND"
