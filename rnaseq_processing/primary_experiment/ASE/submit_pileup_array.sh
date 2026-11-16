#!/bin/bash
set -euo pipefail

WORKDIR="/rs/rs_grp_gmccol/GMCcol_ASE_analysis/per_individual_joint_pileup/keep_zero_ALT_sites_analysis/pileups"
COVAR="/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates/GMCcol_deep_covariates.txt"

cd "$WORKDIR"

mkdir -p hets logs

# sample → individual map
tail -n +2 "$COVAR" | cut -f2,11 > sample_individual_map.txt

cut -f2 sample_individual_map.txt | sort -u > individuals.txt
N=$(wc -l < individuals.txt)

STRIP_CHR_MAP="hets/strip_chr_map.txt"
###############################################
# chr -> no-chr map
###############################################
if [[ ! -f "$STRIP_CHR_MAP" ]]; then
  cat > "$STRIP_CHR_MAP" <<'EOF'
chr1 1
chr2 2
chr3 3
chr4 4
chr5 5
chr6 6
chr7 7
chr8 8
chr9 9
chr10 10
chr11 11
chr12 12
chr13 13
chr14 14
chr15 15
chr16 16
chr17 17
chr18 18
chr19 19
chr20 20
chr21 21
chr22 22
chrX X
chrY Y
chrM MT
chrMT MT
EOF
fi

chmod +x run_pileup_individual.sh  get_het_sites.sh

echo "Submitting array with $N individuals"
sbatch --array=1-"$N" run_pileup_individual.sh
