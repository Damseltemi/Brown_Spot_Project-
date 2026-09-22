#!/bin/bash
#SBATCH --job-name=sydowia_popgen
#SBATCH --output=logs/sydowia_popgen_%j.out
#SBATCH --error=logs/sydowia_popgen_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH -p jrw0107_std,general

set -euo pipefail

# === INPUTS ===
VCF="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/filtered_variants.vcf.gz"
OUTDIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/hetero_homozy"

mkdir -p "$OUTDIR" logs

# === LOAD MODULES ===
module load bcftools
module load vcftools
module load plink/1.9
module load htslib      # <--- THIS PROVIDES TABIX

echo "### 1. Extract high-quality biallelic SNPs"
bcftools view \
  -m2 -M2 -v snps \
  -i 'QUAL>30 && DP>4 && F_MISSING<0.2' \
  -Oz -o "$OUTDIR/analysis_snps.vcf.gz" \
  "$VCF"

# Index the VCF
tabix -p vcf "$OUTDIR/analysis_snps.vcf.gz"

echo "### 2. Individual heterozygosity"
vcftools \
  --gzvcf "$OUTDIR/analysis_snps.vcf.gz" \
  --het \
  --out "$OUTDIR/het_individuals"

echo "### 3. Convert to PLINK format"
plink \
  --vcf "$OUTDIR/analysis_snps.vcf.gz" \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out "$OUTDIR/Sydowia"

echo "### 4. Run ROH detection"
plink \
  --bfile "$OUTDIR/Sydowia" \
  --allow-extra-chr \
  --homozyg \
  --homozyg-snp 50 \
  --homozyg-kb 50 \
  --homozyg-density 50 \
  --homozyg-window-snp 50 \
  --homozyg-window-het 1 \
  --homozyg-window-missing 5 \
  --out "$OUTDIR/Sydowia_roh"

echo "Done 🎉 Population genetics completed!"
