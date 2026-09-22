#!/bin/bash

# ============================================================
# Sydowia polyspora Genome-wide SNP QC + LD Pruning + PCA
# ============================================================
#
# INPUT:
#   analysis_snps.vcf.gz
#
# OPTIONAL:
#   analysis_snps.vcf.gz.tbi
#   Sydowia_metadata.tsv
#
# OUTPUT:
#   PCA/
#       01_raw_plink/
#       02_qc/
#       03_ld_pruned/
#       04_pca/
#
# ============================================================

set -euo pipefail

echo "============================================================"
echo " Sydowia polyspora SNP QC + LD Pruning + PCA"
echo "============================================================"

# ============================================================
# 1. USER SETTINGS
# ============================================================

# Directory containing your VCF
WORKDIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/hetero_homozy1"

# Input VCF
VCF="$WORKDIR/analysis_snps.vcf.gz"

# Output directory
OUTDIR="$WORKDIR/PCA"

# Metadata file
# Set to "" if you do not have metadata yet.
METADATA="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/Sydowiaregion.txt"

# ============================================================
# SNP QC PARAMETERS
# ============================================================

# Maximum missingness per SNP
# 0.10 = retain SNPs present in >=90% of samples
MAX_MISSING=0.10

# Minimum minor allele frequency
# With only 25 isolates, 0.02 is approximately equivalent
# to retaining alleles observed at low frequency.
MIN_MAF=0.02

# LD pruning parameters
WINDOW=50
STEP=5
R2=0.20

# Number of PCs to calculate
NPC=10


# ============================================================
# 2. CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p "$OUTDIR"

mkdir -p "$OUTDIR/01_raw_plink"
mkdir -p "$OUTDIR/02_qc"
mkdir -p "$OUTDIR/03_ld_pruned"
mkdir -p "$OUTDIR/04_pca"
mkdir -p "$OUTDIR/logs"

echo ""
echo "Working directory:"
echo "$WORKDIR"

echo ""
echo "Input VCF:"
echo "$VCF"

# ============================================================
# 3. LOAD SOFTWARE
# ============================================================

echo ""
echo "============================================================"
echo "Loading software"
echo "============================================================"

module load bcftools
module load htslib
module load plink/1.9

echo ""
echo "bcftools:"
bcftools --version | head -1

echo ""
echo "PLINK:"
plink --version

# ============================================================
# 4. CHECK INPUT VCF
# ============================================================

echo ""
echo "============================================================"
echo "Checking input VCF"
echo "============================================================"

if [[ ! -f "$VCF" ]]; then
    echo "ERROR: VCF not found:"
    echo "$VCF"
    exit 1
fi

echo "VCF found."

# ============================================================
# 5. INDEX VCF IF NECESSARY
# ============================================================

if [[ ! -f "${VCF}.tbi" ]]; then

    echo ""
    echo "VCF index not found."
    echo "Creating tabix index..."

    tabix -p vcf "$VCF"

else

    echo ""
    echo "VCF index found."

fi

# ============================================================
# 6. BASIC VCF SUMMARY
# ============================================================

echo ""
echo "============================================================"
echo "Basic VCF statistics"
echo "============================================================"

bcftools stats "$VCF" \
    > "$OUTDIR/01_raw_vcf_stats.txt"

echo ""
echo "Number of samples:"
bcftools query -l "$VCF" | wc -l

echo ""
echo "Sample names:"
bcftools query -l "$VCF"

echo ""
echo "Number of SNP records:"
bcftools view -H "$VCF" | wc -l

# ============================================================
# 7. SAVE SAMPLE NAMES
# ============================================================

bcftools query -l "$VCF" \
    > "$OUTDIR/sample_names.txt"

# ============================================================
# 8. CONVERT VCF TO PLINK
# ============================================================

echo ""
echo "============================================================"
echo "Converting VCF to PLINK"
echo "============================================================"

plink \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --out "$OUTDIR/01_raw_plink/Sydowia_raw"

echo ""
echo "Raw PLINK dataset created."

# ============================================================
# 9. REPORT RAW SNP NUMBER
# ============================================================

echo ""
echo "============================================================"
echo "Raw SNP count"
echo "============================================================"

RAW_SNPS=$(wc -l < "$OUTDIR/01_raw_plink/Sydowia_raw.bim")

echo "Raw SNPs retained by PLINK: $RAW_SNPS"

# ============================================================
# 10. SNP QC
# ============================================================

echo ""
echo "============================================================"
echo "Applying SNP QC"
echo "============================================================"

echo ""
echo "Missingness threshold:"
echo "$MAX_MISSING"

echo ""
echo "MAF threshold:"
echo "$MIN_MAF"

plink \
    --bfile "$OUTDIR/01_raw_plink/Sydowia_raw" \
    --allow-extra-chr \
    --geno "$MAX_MISSING" \
    --maf "$MIN_MAF" \
    --make-bed \
    --out "$OUTDIR/02_qc/Sydowia_QC"

echo ""
echo "SNP QC completed."

# ============================================================
# 11. QC SNP COUNT
# ============================================================

QC_SNPS=$(wc -l < "$OUTDIR/02_qc/Sydowia_QC.bim")

echo ""
echo "============================================================"
echo "SNP counts after QC"
echo "============================================================"

echo "Raw SNPs: $RAW_SNPS"
echo "QC SNPs:  $QC_SNPS"

# Calculate percentage retained
QC_PERCENT=$(awk -v q="$QC_SNPS" -v r="$RAW_SNPS" \
    'BEGIN {printf "%.2f", (q/r)*100}')

echo "Percentage retained: ${QC_PERCENT}%"

# ============================================================
# 12. ADDITIONAL MISSINGNESS REPORT
# ============================================================

echo ""
echo "============================================================"
echo "Generating missingness reports"
echo "============================================================"

plink \
    --bfile "$OUTDIR/02_qc/Sydowia_QC" \
    --allow-extra-chr \
    --missing \
    --out "$OUTDIR/02_qc/Sydowia_missingness"

# ============================================================
# 13. ALLELE FREQUENCY REPORT
# ============================================================

echo ""
echo "============================================================"
echo "Generating allele-frequency report"
echo "============================================================"

plink \
    --bfile "$OUTDIR/02_qc/Sydowia_QC" \
    --allow-extra-chr \
    --freq \
    --out "$OUTDIR/02_qc/Sydowia_frequency"

# ============================================================
# 14. LD PRUNING
# ============================================================

echo ""
echo "============================================================"
echo "LD pruning"
echo "============================================================"

echo ""
echo "Parameters:"
echo "Window = $WINDOW SNPs"
echo "Step   = $STEP SNPs"
echo "r2     = $R2"

plink \
    --bfile "$OUTDIR/02_qc/Sydowia_QC" \
    --allow-extra-chr \
    --indep-pairwise "$WINDOW" "$STEP" "$R2" \
    --out "$OUTDIR/03_ld_pruned/Sydowia_LD"

# ============================================================
# 15. REPORT LD-PRUNED SNP NUMBER
# ============================================================

PRUNED_SNPS=$(wc -l < "$OUTDIR/03_ld_pruned/Sydowia_LD.prune.in")

echo ""
echo "============================================================"
echo "LD pruning results"
echo "============================================================"

echo "SNPs before LD pruning: $QC_SNPS"
echo "Approximately independent SNPs: $PRUNED_SNPS"

# ============================================================
# 16. CREATE PCA DATASET
# ============================================================

echo ""
echo "============================================================"
echo "Creating LD-pruned PCA dataset"
echo "============================================================"

plink \
    --bfile "$OUTDIR/02_qc/Sydowia_QC" \
    --extract "$OUTDIR/03_ld_pruned/Sydowia_LD.prune.in" \
    --allow-extra-chr \
    --make-bed \
    --out "$OUTDIR/03_ld_pruned/Sydowia_PCA"

# ============================================================
# 17. RUN PCA
# ============================================================

echo ""
echo "============================================================"
echo "Running PCA"
echo "============================================================"

plink \
    --bfile "$OUTDIR/03_ld_pruned/Sydowia_PCA" \
    --allow-extra-chr \
    --pca "$NPC" \
    --out "$OUTDIR/04_pca/Sydowia_PCA"

echo ""
echo "PCA completed."

# ============================================================
# 18. CALCULATE VARIANCE EXPLAINED
# ============================================================

echo ""
echo "============================================================"
echo "Calculating variance explained"
echo "============================================================"

awk '
{
    sum += $1
}
END {
    print sum > "'"$OUTDIR"'/04_pca/eigenvalue_sum.txt"
}' "$OUTDIR/04_pca/Sydowia_PCA.eigenval"

TOTAL_EIGEN=$(cat "$OUTDIR/04_pca/eigenvalue_sum.txt")

awk -v total="$TOTAL_EIGEN" '
BEGIN {
    print "PC\tEigenvalue\tPercent_variance"
}
{
    printf "PC%d\t%.6f\t%.4f\n", NR, $1, ($1/total)*100
}
' "$OUTDIR/04_pca/Sydowia_PCA.eigenval" \
> "$OUTDIR/04_pca/Sydowia_PCA_variance.txt"

# ============================================================
# 19. COPY SAMPLE NAMES + PCA RESULTS TO A SIMPLE TABLE
# ============================================================

echo ""
echo "============================================================"
echo "Preparing PCA table"
echo "============================================================"

awk '
BEGIN {
    print "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10"
}
{
    print
}
' "$OUTDIR/04_pca/Sydowia_PCA.eigenvec" \
> "$OUTDIR/04_pca/Sydowia_PCA_table.tsv"

# ============================================================
# 20. OPTIONAL METADATA CHECK
# ============================================================

if [[ -f "$METADATA" ]]; then

    echo ""
    echo "============================================================"
    echo "Metadata file found"
    echo "============================================================"

    echo "Metadata:"
    echo "$METADATA"

    echo ""
    echo "Metadata preview:"
    head "$METADATA"

else

    echo ""
    echo "============================================================"
    echo "No metadata file found"
    echo "============================================================"

    echo "PCA will still be completed."
    echo ""
    echo "To interpret PCA by Region or Severity, create:"
    echo "$METADATA"

fi

# ============================================================
# 21. FINAL SUMMARY
# ============================================================

echo ""
echo "============================================================"
echo "ANALYSIS COMPLETE"
echo "============================================================"

echo ""
echo "Input VCF:"
echo "$VCF"

echo ""
echo "Number of samples:"
bcftools query -l "$VCF" | wc -l

echo ""
echo "Raw SNPs:"
echo "$RAW_SNPS"

echo ""
echo "SNPs after QC:"
echo "$QC_SNPS"

echo ""
echo "LD-pruned SNPs:"
echo "$PRUNED_SNPS"

echo ""
echo "PCA results:"
echo "$OUTDIR/04_pca/Sydowia_PCA.eigenvec"

echo ""
echo "Variance explained:"
echo "$OUTDIR/04_pca/Sydowia_PCA_variance.txt"

echo ""
echo "PCA table:"
echo "$OUTDIR/04_pca/Sydowia_PCA_table.tsv"

echo ""
echo "============================================================"
echo "IMPORTANT:"
echo "PCA is being treated as an exploratory analysis of"
echo "genome-wide genetic relationships among recovered isolates."
echo "Heterozygosity and ROH analyses are NOT included."
echo "============================================================"