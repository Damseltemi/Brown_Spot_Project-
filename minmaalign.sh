#!/bin/bash
#SBATCH --job-name=minimap2_alignment
#SBATCH --output=logs/minimap2_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00
#SBATCH --mem=80G
#SBATCH --partition=jrw0107_std,general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=trf0019@auburn.edu

set -euo pipefail

# Activate conda env
source /home/trf0019/miniforge3/bin/activate minimap2_env

# === CONFIGURATION ===
REF_GENOME="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/goodbusco/BS06_240_394-4.scaffolds.fa" 
READS_DIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq"
OUTPUT_DIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentgood"
THREADS=8

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# === STEP 1: Index reference ===
if [[ ! -f "$OUTPUT_DIR/ref.mmi" ]]; then
    echo "📦 Indexing reference genome..."
    minimap2 -d "$OUTPUT_DIR/ref.mmi" "$REF_GENOME"
else
    echo "📦 Using existing index: $OUTPUT_DIR/ref.mmi"
fi

# === STEP 2: Loop through FASTQ and FASTQ.GZ ===
echo "🔍 Searching for FASTQ and FASTQ.GZ files..."

shopt -s nullglob

for READ in "$READS_DIR"/*.fastq "$READS_DIR"/*.fastq.gz; do
    
    # skip loop if no files found
    [[ -e "$READ" ]] || { echo "❌ No FASTQ(.gz) files found in $READS_DIR"; break; }

    # Get name without extension
    BASENAME=$(basename "$READ")
    BASENAME=${BASENAME%.fastq}
    BASENAME=${BASENAME%.fastq.gz}

    SAMPLE_OUT="$OUTPUT_DIR/$BASENAME"

    echo "🔁 Processing sample: $BASENAME"

    # === ALIGN ===
    minimap2 -ax map-ont -t "$THREADS" "$OUTPUT_DIR/ref.mmi" "$READ" > "$SAMPLE_OUT.sam"

    # === SAM → BAM ===
    samtools view -@ "$THREADS" -bS "$SAMPLE_OUT.sam" \
        | samtools sort -@ "$THREADS" -o "$SAMPLE_OUT.sorted.bam"

    samtools index "$SAMPLE_OUT.sorted.bam"

    rm "$SAMPLE_OUT.sam"

    samtools flagstat "$SAMPLE_OUT.sorted.bam" > "$SAMPLE_OUT.stats.txt"

    echo "✅ Completed $BASENAME — Stats: $SAMPLE_OUT.stats.txt"
done

echo "🎉 All alignments saved in: $OUTPUT_DIR"
