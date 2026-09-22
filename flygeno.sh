#!/usr/bin/env bash
#SBATCH --partition=jrw0107_std,general
#SBATCH --nodes=1                        # Request 1 node
#SBATCH --ntasks-per-node=10             # 20 tasks per node
#SBATCH --mem=50G                        # Request 50 GB of memory
#SBATCH --time=50:00:00                  # Set the time limit to 10 hours
#SBATCH --job-name=flye_assembly         # Job name
#SBATCH --output=flye_output_%j.log      # Output log file
#SBATCH --error=flye_error_%j.log        # Error log file

# Load Flye environment
source /home/trf0019/miniforge3/bin/activate flye_env

# Define paths
INPUT_DIR="/scratch/acicola_project/20260114/20260114_demux/combined_fastq"
OUTPUT_DIR="/scratch/acicola_project/20260114/20260114_demux/combined_fastq/assembly"
SAMPLES_FILE="/scratch/acicola_project/20260114/no_sample_id/20260114_1105_MC-114852_FBC81239_8536e92c/samplesheetccomb.csv"


# Number of threads to use (align with ntasks-per-node)
THREADS=10

# Ensure the output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
  echo "Output directory created: $OUTPUT_DIR"
fi

# **Loop Through Each Sample**
IFS=,
tail -n +2 "$SAMPLES_FILE" | while read -r sample_id alias barcode experiment_id position_id kit; do
SAMPLE_INPUT=$(find "$INPUT_DIR" -type f -name "*${sample_id}*.fastq" | head -n 1)
    SAMPLE_OUTPUT="${OUTPUT_DIR}/${sample_id}"

  # **Check if the input file exists**
  if [ ! -f "$SAMPLE_INPUT" ]; then
    echo "Error: Input file $SAMPLE_INPUT does not exist for sample $sample_id."
    continue  # Make sure this is inside the loop
  fi

  # **Check if assembly already exists and is up-to-date**
  if [[ ! -d "$SAMPLE_OUTPUT" || ! -f "$SAMPLE_OUTPUT/assembly.fasta" || "$SAMPLE_INPUT" -nt "$SAMPLE_OUTPUT/assembly.fasta" ]]; then
    echo "Running Flye for sample $sample_id..."
    flye --nano-hq "$SAMPLE_INPUT" \
         --out-dir "$SAMPLE_OUTPUT" \
         --threads "$THREADS" \
         --scaffold
    echo "Flye assembly completed for sample $sample_id."
  else
    echo "Skipping sample $sample_id: Assembly already exists and is up-to-date."
  fi
done

# Confirm successful execution
echo "Flye genome assembly script completed."
