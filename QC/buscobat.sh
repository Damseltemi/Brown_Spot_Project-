#!/bin/bash
#SBATCH --job-name=BUSCO_batch
#SBATCH --output=BUSCO_batch_%j.log
#SBATCH --partition=jrw0107_std,general
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=10-00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=trf0019@auburn.edu


# Activate BUSCO 6 environment
source /home/trf0019/miniforge3/bin/activate busco6_env

# Input / output directories
INDIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/goodbusco"
OUTPATH="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/busco_runsnw"

# CPUs: use SLURM_CPUS_PER_TASK if set, otherwise default to 15
CPU="${SLURM_CPUS_PER_TASK:-10}"

# Lineage dataset
LINEAGE="/home/trf0019/busco_downloads/lineages/dothideomycetes_odb10"
# Or, alternative:
# LINEAGE="/home/trf0019/busco_downloads/lineages/ascomycota_odb10"

mkdir -p "$OUTPATH"

# Put temp dir on scratch to avoid home quota problems
export TMPDIR="/scratch/acicola_project/tmp/${SLURM_JOB_ID:-busco_tmp}"
mkdir -p "$TMPDIR"

echo "Input dir:  $INDIR"
echo "Output dir: $OUTPATH"
echo "Lineage:    $LINEAGE"
echo "CPUs:       $CPU"
echo

shopt -s nullglob

# Loop over all .fa and .fasta files in INDIR
for fasta in "$INDIR"/*.fa "$INDIR"/*.fasta; do
  # In case one of the globs matches nothing
  [ -e "$fasta" ] || continue

  base="$(basename "$fasta")"
  sample="${base%.*}"   # strips .fa or .fasta

  echo "=================================================="
  echo "Running BUSCO for: $sample"
  echo "Input:  $fasta"
  echo "Output: ${OUTPATH}/${sample}"
  echo "=================================================="

  busco -i "$fasta" \
        -o "$sample" \
        --out_path "$OUTPATH" \
        -m genome \
        -l "$LINEAGE" \
        -c "$CPU"

  echo "BUSCO finished: $sample"
  echo
done

shopt -u nullglob

echo "All BUSCO runs completed."
