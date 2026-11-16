#!/bin/bash
set -euo pipefail

WORKDIR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/pileups"
COVAR="/rs/rs_grp_ibdeqtl/GxM_Validation/IBD_GxM_Validation/IBDGXM_validation_covariates_updated_10-14-25.txt"

cd "$WORKDIR"

mkdir -p logs hets

# sample → individual map
tail -n +2 "$COVAR" | cut -f2,16 > sample_individual_map.txt

cut -f2 sample_individual_map.txt | sort -u > individuals.txt
N=$(wc -l < individuals.txt)

STRIP_CHR_MAP="hets/strip_chr_map.txt"

# chr -> no-chr map
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


chmod +x run_individual_multisample_pileup.sh get_het_sites.sh

echo "Submitting array with $N individuals"
sbatch --array=1-"$N" run_individual_multisample_pileup.sh
