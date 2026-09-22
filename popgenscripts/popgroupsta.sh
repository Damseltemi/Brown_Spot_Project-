#!/bin/bash
#SBATCH --job-name=sydowia_diversity
#SBATCH --output=sydowia_diversity_%j.out
#SBATCH --error=sydowia_diversity_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=60:00:00
#SBATCH -p jrw0107_std,general

set -euo pipefail

module load bcftools
module load vcftools

# -------- PATHS --------
BASE="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results"
POP_DIR="/scratch/acicola_project/fastas/sydowiaspp/vcf_results/fst_pops"
OUTDIR="${BASE}/diversity_statsnew"
WINDOW=5000

mkdir -p "$OUTDIR"

# Main VCFs
VCF_PI_TD="${BASE}/filtered_variants.vcf.gz"   # for pi & TajimaD
VCF_FST="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/hetero_homozy/analysis_snps.vcf.gz"         # QC-filtered SNPs for FST

# -------- π AND TAJIMA'S D PER GROUP --------
# Group names must match your pop_*.txt files in $POP_DIR

# 1) North_LowModerate
vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_North_LowModerate.txt" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_North_LowModerate_${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_North_LowModerate.txt" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_North_LowModerate_${WINDOW}bp"

# 2) South_LowModerate
vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_South_LowModerate.txt" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_South_LowModerate_${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_South_LowModerate.txt" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_South_LowModerate_${WINDOW}bp"

# 3) North_High
vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_North_High.txt" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_North_High_${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_North_High.txt" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_North_High_${WINDOW}bp"

# 4) South_High
vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_South_High.txt" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_South_High_${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/pop_South_High.txt" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_South_High_${WINDOW}bp"
# 5) OAWEST
vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/OAWEST.txt" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_OAWEST_${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/OAWEST.txt" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_OAWEST_${WINDOW}bp"

  # 6) South
vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/South" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_South${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/South" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_South${WINDOW}bp"

 # 7) North

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/North" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_North_${WINDOW}bp"

vcftools \
  --gzvcf "$VCF_PI_TD" \
  --keep "${POP_DIR}/North" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajima_North_${WINDOW}bp"

# -------- FST: PAIRWISE COMPARISONS --------
# You must run FST once per pair of populations.

# Low/Moderate:  South region
vcftools \
  --gzvcf "$VCF_FST" \
  --weir-fst-pop "${POP_DIR}/pop_North_LowModerate.txt" \
  --weir-fst-pop "${POP_DIR}/pop_South_LowModerate.txt" \
  --fst-window-size $WINDOW \
  --fst-window-step $WINDOW \
  --out "${OUTDIR}/fst_South_LowModHigh_${WINDOW}bp"

# High severity: North 
vcftools \
  --gzvcf "$VCF_FST" \
  --weir-fst-pop "${POP_DIR}/pop_North_High.txt" \
  --weir-fst-pop "${POP_DIR}/pop_South_High.txt" \
  --fst-window-size $WINDOW \
  --fst-window-step $WINDOW \
  --out "${OUTDIR}/fst_North_LowModHigh_${WINDOW}bp"

  # Region
vcftools \
  --gzvcf "$VCF_FST" \
  --weir-fst-pop "${POP_DIR}/North" \
  --weir-fst-pop "${POP_DIR}/South" \
  --fst-window-size $WINDOW \
  --fst-window-step $WINDOW \
  --out "${OUTDIR}/fst_region_${WINDOW}bp"

    # Region
vcftools \
  --gzvcf "$VCF_FST" \
  --weir-fst-pop "${POP_DIR}/North" \
  --weir-fst-pop "${POP_DIR}/OAWEST" \
  --fst-window-size $WINDOW \
  --fst-window-step $WINDOW \
  --out "${OUTDIR}/fst_region_${WINDOW}bp"

    # Region
vcftools \
  --gzvcf "$VCF_FST" \
  --weir-fst-pop "${POP_DIR}/OAWEST" \
  --weir-fst-pop "${POP_DIR}/South" \
  --fst-window-size $WINDOW \
  --fst-window-step $WINDOW \
  --out "${OUTDIR}/fst_region_${WINDOW}bp"

echo "Done ✅  π, Tajima's D, and FST window stats written to: $OUTDIR"