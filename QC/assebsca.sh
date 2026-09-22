#!/usr/bin/env bash
#SBATCH --partition=gpu2,gpu4            # Specify the GPU partition
#SBATCH --gres=gpu:tesla:2               # Request 2 Tesla GPUs
#SBATCH --nodes=1                        # Request 1 node
#SBATCH --ntasks-per-node=20             # 20 tasks per node
#SBATCH --mem=50G                        # Request 50 GB of memory
#SBATCH --time=40:00:00                  # Set the time limit to 10 hours
#SBATCH --job-name=var_assessment      # Job name
#SBATCH --output=logs/vanrt.out
#SBATCH --error=logs/vanrt.err

# Set directories
RAW_DIR="/scratch/acicola_project/fastqs"
FILTERED_DIR="/scratch/acicola_project/filtered_reads"
ASSEMBLY_DIR="/scratch/acicola_project/assemblies"
SCAFFOLD_DIR="/scratch/acicola_project/scaffolded"

source 
mkdir -p "$FILTERED_DIR" "$ASSEMBLY_DIR" "$SCAFFOLD_DIR"

# Step 1: Quality Filtering with Filtlong
echo "Running Filtlong filtering..."
for fq in "$RAW_DIR"/*fastq*; do
    base=$(basename "$fq" | sed 's/\.fastq\.gz//;s/\.fastq//')
    filtlong --min_length 1000 --keep_percent 90 "$fq" > "$FILTERED_DIR/${base}_filtered.fastq"
done

# Step 2: Assembly with Flye
echo "Running Flye assembly..."
for fq in "$FILTERED_DIR"/*_filtered.fastq; do
    sample=$(basename "$fq" _filtered.fastq)
    flye --nano-raw "$fq" --out-dir "$ASSEMBLY_DIR/$sample" --genome-size 40m --threads 12
done

# Step 3: Optional Polishing with Medaka
echo "Running Medaka polishing..."
for fq in "$FILTERED_DIR"/*_filtered.fastq; do
    sample=$(basename "$fq" _filtered.fastq)
    medaka_consensus -i "$fq" \
                     -d "$ASSEMBLY_DIR/$sample/assembly.fasta" \
                     -o "$SCAFFOLD_DIR/${sample}_medaka" \
                     -t 8 -m r941_min_high_g303
done

echo "Pipeline completed!"