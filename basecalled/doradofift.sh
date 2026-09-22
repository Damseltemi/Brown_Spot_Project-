#!/usr/bin/env bash
#SBATCH --job-name=dorado_basecalling
#SBATCH --output=logs/dorado_basecall.out
#SBATCH --error=logs/dorado_basecall.err
#SBATCH --partition=gpu2,gpu4
#SBATCH --gres=gpu:tesla:2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40G
#SBATCH --time=2-00:00:00 
#SBATCH --mail-type=END,FAIL                 # Email user when job finishes or fails
#SBATCH --mail-user=trf0019@auburn.edu

# Load necessary modules for CUDA
module load cuda11.0/toolkit/11.8.0

# Set environment variable for Dorado library
export LD_LIBRARY_PATH="/home/trf0019/software_app/dorado-0.7.1-linux-x64/lib:$LD_LIBRARY_PATH"

# Set paths
POD5_DIR="/scratch/acicola_project/20260114/no_sample_id/20260114_1105_MC-114852_FBC81239_8536e92c/pod5_skip"
SAMPLE_SHEET="/scratch/acicola_project/20260114/no_sample_id/20260114_1105_MC-114852_FBC81239_8536e92c/samplesheet2026.csv"
OUTPUT_BAM="/scratch/acicola_project/20260114/20251106.bam"

# Run the basecaller
/home/trf0019/software_app/dorado-0.7.1-linux-x64/bin/dorado \
    basecaller sup "$POD5_DIR" \
    --kit-name "SQK-NBD114-24" \
    --sample-sheet "$SAMPLE_SHEET" > "$OUTPUT_BAM"
