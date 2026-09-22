#!/bin/bash
#SBATCH --job-name=diversity_stats
#SBATCH --output=logs/diversity_stats_%j.out
#SBATCH --error=logs/diversity_stats_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=60:00:00
#SBATCH -p jrw0107_std,general

# ===== Modules / env =====
module load bcftools
module load vcftools

# ===== CONFIG =====
# Original merged VCF (from bcftools mpileup/call)
RAW_VCF="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/merged_variants_fixed.vcf.gz"

# Base working dir for all downstream files
BASE_OUT="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results"

# Renamed + filtered VCFs
RENAMED_VCF="${BASE_OUT}/merged_variants_fixed_renamed.vcf.gz"
FILTERED_VCF="${BASE_OUT}/filtered_variants.vcf.gz"

# Diversity output directory
OUTDIR="${BASE_OUT}/diversity_stats"

# Window size (bp) for all window-based stats
WINDOW=5000

mkdir -p "$BASE_OUT" "$OUTDIR" logs

echo "=== Diversity stats pipeline ==="
echo "Raw VCF:        $RAW_VCF"
echo "Renamed VCF:    $RENAMED_VCF"
echo "Filtered VCF:   $FILTERED_VCF"
echo "Diversity out:  $OUTDIR"
echo "Window size:    $WINDOW bp"
echo

############################
# 1. Build rename.txt file #
############################
RENAME_FILE="${BASE_OUT}/rename.txt"

echo "[1/4] Creating rename file: $RENAME_FILE"
bcftools query -l "$RAW_VCF" \
  | sed 's#.*/##; s/_trimmed\.sorted\.bam$//; s/_combined\.sorted\.bam$//' \
  > "$RENAME_FILE"

echo "First few sample names BEFORE renaming:"
bcftools query -l "$RAW_VCF" | head
echo
echo "First few sample names AFTER renaming:"
head "$RENAME_FILE"
echo

##################################
# 2. Reheader VCF with new names #
##################################
echo "[2/4] Reheadering VCF with new sample names..."
bcftools reheader \
  -s "$RENAME_FILE" \
  -o "$RENAMED_VCF" \
  "$RAW_VCF"

tabix -p vcf "$RENAMED_VCF"

echo "Samples in renamed VCF:"
bcftools query -l "$RENAMED_VCF"
echo

#########################################
# 3. Simple filtering (SNPs, bi-allelic) #
#########################################
echo "[3/4] Filtering VCF to SNPs, bi-allelic..."
bcftools view \
  -m2 -M2 -v snps \
  -Oz \
  -o "$FILTERED_VCF" \
  "$RENAMED_VCF"

tabix -p vcf "$FILTERED_VCF"

echo "Filtered VCF created: $FILTERED_VCF"
echo

##########################################
# 4. Window-based diversity statistics   #
##########################################
VCF_FOR_STATS="$FILTERED_VCF"

echo "[4/4] Calculating diversity statistics from: $VCF_FOR_STATS"

# ---- π (Nucleotide diversity) ----
echo "  • π (nucleotide diversity)..."
vcftools \
  --gzvcf "$VCF_FOR_STATS" \
  --window-pi $WINDOW \
  --out "${OUTDIR}/pi_${WINDOW}bp"

# Creates: ${OUTDIR}/pi_${WINDOW}bp.windowed.pi

# ---- Watterson’s θw (if your vcftools has it) ----
# If this fails with “Unknown option: --window-theta”, just comment this block out.
echo "  • Watterson's theta (if supported by vcftools)..."
vcftools \
  --gzvcf "$VCF_FOR_STATS" \
  --window-theta $WINDOW \
  --out "${OUTDIR}/thetaW_${WINDOW}bp" || echo "    (window-theta not supported in this vcftools build)"

# Creates: ${OUTDIR}/thetaW_${WINDOW}bp.windowed.theta

# ---- Tajima’s D ----
echo "  • Tajima's D..."
vcftools \
  --gzvcf "$VCF_FOR_STATS" \
  --TajimaD $WINDOW \
  --out "${OUTDIR}/tajimaD_${WINDOW}bp"

# Creates: ${OUTDIR}/tajimaD_${WINDOW}bp.Tajima.D

echo
echo "=== Done. Diversity files in: $OUTDIR ==="
