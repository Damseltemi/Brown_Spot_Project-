#!/bin/bash
#SBATCH --job-name=fastqc_job
#SBATCH --output=/scratch/trf0019/fastqc_job.out
#SBATCH --error=/scratch/trf0019/fastqc_job.err
#SBATCH --time=01:00:00
#SBATCH --mem=100G  # Increase the memory allocation

# Load the FastQC module
module load fastqc

# Create the output directory if it doesn't exist
mkdir -p /scratch/trf0019/Acicola_seq

# Directory containing FASTQ files
FASTQ_DIR=/home/trf0019/Acicola_seq/fastq
OUTPUT_DIR=/scratch/trf0019/Acicola_seq

# Increase Java heap space
export JAVA_OPTS="-Xmx12g"

# Find all files ending with fastq.gz in the directory
FILES=(${FASTQ_DIR}/*.fastq.gz)

# Check if files exist and have read permissions
for FILE in "${FILES[@]}"; do
    if [[ -r "${FASTQ_DIR}/${FILE}" ]]; then
        echo "Running FastQC on ${FILE}"
        fastqc -t 5 -o "${OUTPUT_DIR}" "${FASTQ_DIR}/${FILE}"
    else
        echo "Skipping '${FILE}' which didn't exist, or couldn't be read"
    fi
done

echo "FastQC analysis completed for all files."






