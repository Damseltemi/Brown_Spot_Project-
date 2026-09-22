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
OUTDIR="${BASE}/diversityfst"
WINDOW=5000

mkdir -p "$OUTDIR"

# VCF for FST (QC-filtered SNPs)
VCF_FST="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood/vcf_results/hetero_homozy/analysis_snps.vcf.gz"

# Population labels (must match keys in POP_FILE below)
POPS=(
  North_High
  North_LowModerate
  South_High
  South_LowModerate
  South
  North
  OAWEST
)

# Map labels -> pop files
declare -A POP_FILE=(
  [North_High]="${POP_DIR}/pop_North_High.txt"
  [North_LowModerate]="${POP_DIR}/pop_North_LowModerate.txt"
  [South_High]="${POP_DIR}/pop_South_High.txt"
  [South_LowModerate]="${POP_DIR}/pop_South_LowModerate.txt"
  [South]="${POP_DIR}/South"
  [North]="${POP_DIR}/North"
  [OAWEST]="${POP_DIR}/OAWEST.txt"
)

# -------- SANITY CHECKS --------
echo "[Check] VCF: $VCF_FST"
[[ -f "$VCF_FST" ]] || { echo "ERROR: Missing VCF_FST: $VCF_FST"; exit 1; }
[[ -f "${VCF_FST}.tbi" ]] || { echo "ERROR: Missing index ${VCF_FST}.tbi (run: tabix -p vcf $VCF_FST)"; exit 1; }

echo "[Check] Population files:"
for p in "${POPS[@]}"; do
  f="${POP_FILE[$p]:-}"
  if [[ -z "$f" ]]; then
    echo "ERROR: No POP_FILE mapping for '$p'"
    exit 1
  fi
  if [[ ! -s "$f" ]]; then
    echo "ERROR: Pop file missing or empty for '$p': $f"
    exit 1
  fi
  echo "  - $p -> $f ($(wc -l < "$f") samples)"
done

# Optional: ensure sample IDs exist in VCF (prevents silent empty pops)
echo "[Check] Verifying pop IDs exist in VCF header..."
VCF_SAMPLES="$(mktemp)"
bcftools query -l "$VCF_FST" > "$VCF_SAMPLES"

for p in "${POPS[@]}"; do
  f="${POP_FILE[$p]}"
  missing_count=$(grep -Fvx -f "$VCF_SAMPLES" "$f" | wc -l || true)
  if [[ "$missing_count" -gt 0 ]]; then
    echo "ERROR: $missing_count sample IDs in $f not found in VCF samples."
    echo "First few missing:"
    grep -Fvx -f "$VCF_SAMPLES" "$f" | head
    rm -f "$VCF_SAMPLES"
    exit 1
  fi
done
rm -f "$VCF_SAMPLES"

echo "Starting pairwise FST calculations..."
echo "Window size = $WINDOW bp; step = $WINDOW bp"

# Loop over all unique pairs i < j (includes ALL groups in POPS)
for (( i=0; i<${#POPS[@]}-1; i++ )); do
  for (( j=i+1; j<${#POPS[@]}; j++ )); do

    POPA="${POPS[$i]}"
    POPB="${POPS[$j]}"

    FILEA="${POP_FILE[$POPA]}"
    FILEB="${POP_FILE[$POPB]}"

    # Clean prefix (vcftools writes multiple files with this prefix)
    OUT_PREFIX="${OUTDIR}/fst_${POPA}_vs_${POPB}_${WINDOW}bp"

    # Skip if already completed (vcftools produces *.windowed.weir.fst for windowed)
    if [[ -s "${OUT_PREFIX}.windowed.weir.fst" ]]; then
      echo "SKIP (exists): ${POPA} vs ${POPB} -> ${OUT_PREFIX}.windowed.weir.fst"
      continue
    fi

    echo "Running FST: ${POPA} vs ${POPB} -> ${OUT_PREFIX}"

    vcftools \
      --gzvcf "$VCF_FST" \
      --weir-fst-pop "$FILEA" \
      --weir-fst-pop "$FILEB" \
      --fst-window-size "$WINDOW" \
      --fst-window-step "$WINDOW" \
      --out "$OUT_PREFIX"

    # Hard fail if output not produced (prevents “finished” when nothing ran)
    if [[ ! -s "${OUT_PREFIX}.windowed.weir.fst" ]]; then
      echo "ERROR: No output produced for ${POPA} vs ${POPB} (expected ${OUT_PREFIX}.windowed.weir.fst)"
      exit 1
    fi

  done
done

echo "All pairwise FST calculations completed."
