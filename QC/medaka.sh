#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH --output=medaka_%j.log
#SBATCH --partition=jrw0107_std,general
#SBATCH --cpus-per-task=15
#SBATCH --mem=80G
#SBATCH --time=5-00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=trf0019@auburn.edu

set -euo pipefail

source /home/trf0019/miniforge3/bin/activate medaka_clean

ASSEMBLY_DIR="/scratch/acicola_project/20260114/20260114_demux/combined_fastq/assembly/final_fastas"
READS_DIR="/scratch/acicola_project/20260114/20260114_demux/combined_fastq"
OUTROOT="/scratch/acicola_project/20260114/20260114_demux/assembly/medaka_polished"

THREADS="${SLURM_CPUS_PER_TASK:-15}"
mkdir -p "$OUTROOT"

###############################################################################
# Get valid model list (one per line) and choose a good one for Dorado SUP R10.4.1
###############################################################################
MODELS_FILE="$OUTROOT/medaka_models.$SLURM_JOB_ID.txt"

# Convert "Available: a, b, c" into one-model-per-line:
medaka tools list_models \
  | sed 's/^Available:[[:space:]]*//' \
  | tr ',' '\n' \
  | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
  | grep -v '^$' > "$MODELS_FILE"

echo "Found $(wc -l < "$MODELS_FILE") medaka models"
echo "First 10 models:"
head -n 10 "$MODELS_FILE"

# Prefer newest SUP 400bps R10.4.1 model
preferred_order=(
  "r1041_e82_400bps_sup_v5.2.0"
  "r1041_e82_400bps_sup_v5.0.0"
  "r1041_e82_400bps_sup_v4.3.0"
  "r1041_e82_400bps_sup_v4.2.0"
  "r1041_e82_400bps_sup_v4.1.0"
  "r1041_e82_400bps_sup_v4.0.0"
)

MODEL=""
for m in "${preferred_order[@]}"; do
  if grep -qx "$m" "$MODELS_FILE"; then
    MODEL="$m"
    break
  fi
done

# If none matched exactly, fall back to the first SUP R1041 400bps model
if [[ -z "$MODEL" ]]; then
  MODEL="$(grep -i '^r1041_.*400bps_.*sup' "$MODELS_FILE" | head -n 1 || true)"
fi

if [[ -z "$MODEL" ]]; then
  echo "ERROR: Could not find an r1041*400bps*sup* model."
  echo "Try running: medaka tools list_models | head"
  exit 1
fi

echo "Using medaka model: $MODEL"

###############################################################################
# Loop over assemblies and polish (resumable)
###############################################################################
shopt -s nullglob
for assembly in "$ASSEMBLY_DIR"/*.fasta; do
  sample="$(basename "$assembly" .fasta)"

  # Find matching reads (your files look like: <sample>.fastq or <sample>.fastq.gz)
  reads=""
  if [[ -f "$READS_DIR/${sample}.fastq.gz" ]]; then
    reads="$READS_DIR/${sample}.fastq.gz"
  elif [[ -f "$READS_DIR/${sample}.fastq" ]]; then
    reads="$READS_DIR/${sample}.fastq"
  else
    # also allow a ".unique" naming mismatch between fasta and reads
    if [[ -f "$READS_DIR/${sample}.fastq" ]]; then
      reads="$READS_DIR/${sample}.fastq"
    elif [[ -f "$READS_DIR/${sample}.fastq.gz" ]]; then
      reads="$READS_DIR/${sample}.fastq.gz"
    else
      echo "WARNING: No reads found for sample '$sample' in $READS_DIR. Skipping."
      continue
    fi
  fi

  outdir="$OUTROOT/$sample"
  mkdir -p "$outdir"

  # Resume: skip if consensus exists and non-empty
  if [[ -s "$outdir/consensus.fasta" ]]; then
    echo "SKIP (already polished): $sample"
    continue
  fi

  echo "=================================================="
  echo "Polishing: $sample"
  echo "Assembly : $assembly"
  echo "Reads    : $reads"
  echo "Outdir   : $outdir"
  echo "Model    : $MODEL"
  echo "=================================================="

  medaka_consensus \
    -i "$reads" \
    -d "$assembly" \
    -o "$outdir" \
    -t "$THREADS" \
    -m "$MODEL"

  echo "DONE: $sample -> $outdir/consensus.fasta"
done
shopt -u nullglob

echo "All medaka polishing finished."
