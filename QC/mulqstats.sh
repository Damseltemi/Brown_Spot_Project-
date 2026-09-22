#!/bin/bash
#SBATCH --job-name=assembly_stats_summary
#SBATCH --output=assembly_stats_summary.log
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=24G
#SBATCH --time=01:00:00

set -euo pipefail

################################################################################
# Sydowia assembly stats + BUSCO + MultiQC
################################################################################

# === Activate MultiQC environment ===
source /home/trf0019/miniforge3/bin/activate multiqc_env

# === Paths ===
FASTA_DIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/goodbusco"
STATS_DIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/assembly_statsnew1"
BUSCO_DIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafasta/busco_runsnw"

SUMMARY_FILE="${STATS_DIR}/assembly_stats_combined.txt"
MULTIQC_OUTPUT="${STATS_DIR}/multiqc"
CLEAN_BUSCO_DIR="${STATS_DIR}/busco_multiqc_clean"

ASSEMBLY_STATS_PY="/home/trf0019/Working_script/assembly_stats.py"

mkdir -p "${STATS_DIR}" "${MULTIQC_OUTPUT}" "${CLEAN_BUSCO_DIR}"

################################################################################
# Step 1: Generate assembly stats per FASTA (.fa / .fasta)
################################################################################
echo "📊 Generating individual assembly stats from: ${FASTA_DIR}"

shopt -s nullglob
fasta_files=( "${FASTA_DIR}"/*.fa "${FASTA_DIR}"/*.fasta )
shopt -u nullglob

if [[ ${#fasta_files[@]} -eq 0 ]]; then
  echo "ERROR: No .fa or .fasta files found in ${FASTA_DIR}"
  exit 1
fi

for file in "${fasta_files[@]}"; do
  fname="$(basename "${file}")"
  base="${fname%.fasta}"
  base="${base%.fa}"

  echo "  -> ${fname}"
  python "${ASSEMBLY_STATS_PY}" "${file}" > "${STATS_DIR}/${base}_stats.txt"
done

################################################################################
# Step 2: Combine all stats into one summary file
################################################################################
echo "📋 Combining stats into summary table: ${SUMMARY_FILE}"
echo -e "Sample\tTotal_Length\tScaffolds\tContigs\tScaffold_N50\tContig_N50\tGC_Content" > "${SUMMARY_FILE}"

shopt -s nullglob
stat_files=( "${STATS_DIR}"/*_stats.txt )
shopt -u nullglob

if [[ ${#stat_files[@]} -eq 0 ]]; then
  echo "ERROR: No *_stats.txt files produced in ${STATS_DIR}"
  exit 1
fi

for file in "${stat_files[@]}"; do
  sample="$(basename "${file}" _stats.txt)"

  ##########################################################################
  # 2a. Find the corresponding FASTA and count scaffolds (= number of '>')
  ##########################################################################
  fasta_path=""

  # Exact .fa / .fasta first
  if [[ -f "${FASTA_DIR}/${sample}.fa" ]]; then
    fasta_path="${FASTA_DIR}/${sample}.fa"
  elif [[ -f "${FASTA_DIR}/${sample}.fasta" ]]; then
    fasta_path="${FASTA_DIR}/${sample}.fasta"
  else
    # More relaxed search
    fasta_path="$(find "${FASTA_DIR}" -maxdepth 1 -type f \
      \( -name "${sample}.fa" -o -name "${sample}.fasta" \
         -o -name "${sample}.fa.gz" -o -name "${sample}.fasta.gz" \
         -o -name "${sample}.fna" -o -name "${sample}.fna.gz" \) \
      | head -n 1 || true)"
  fi

  scaffolds="NA"
  if [[ -n "${fasta_path}" ]]; then
    if [[ "${fasta_path}" == *.gz ]]; then
      scaffolds="$(zgrep -c '^>' "${fasta_path}" || echo "NA")"
    else
      scaffolds="$(grep -c '^>' "${fasta_path}" || echo "NA")"
    fi
  else
    echo "WARNING: no FASTA found for sample ${sample}" >&2
  fi

  ##########################################################################
  # 2b. Pull the other stats from the *_stats.txt file
  ##########################################################################
  total_length="$(
    awk 'tolower($0) ~ /total[[:space:]]+length/ {print $NF; exit}' "${file}"
  )"

  contigs="$(
    awk 'tolower($0) ~ /number[[:space:]]+of[[:space:]]+contigs/ {print $NF; exit}' "${file}"
  )"

  contig_n50="$(
    awk 'tolower($0) ~ /(^|[[:space:]])n50([[:space:]]|$)/ {print $NF; exit}' "${file}"
  )"

  # We don't have true scaffold N50 from assembly_stats.py; leave as NA
  scaffold_n50="NA"

  gc="$(
    awk 'tolower($0) ~ /gc[[:space:]]+content/ {print $NF; exit}' "${file}"
  )"

  # Fallbacks
  total_length="${total_length:-NA}"
  contigs="${contigs:-NA}"
  contig_n50="${contig_n50:-NA}"
  gc="${gc:-NA}"

  # If scaffold count failed but contig count is available,
  # use contigs as scaffolds (many of your files are scaffold-level FASTA)
  if [[ "${scaffolds}" == "NA" || -z "${scaffolds}" ]]; then
    if [[ "${contigs}" != "NA" && -n "${contigs}" ]]; then
      scaffolds="${contigs}"
    fi
  fi

  echo -e "${sample}\t${total_length}\t${scaffolds}\t${contigs}\t${scaffold_n50}\t${contig_n50}\t${gc}" \
    >> "${SUMMARY_FILE}"
done

################################################################################
# Step 3: Prepare BUSCO short summaries for MultiQC
################################################################################
echo "🧹 Preparing BUSCO inputs for MultiQC..."

shopt -s nullglob
busco_summaries=( "${BUSCO_DIR}"/short_summary*.txt "${BUSCO_DIR}"/*/short_summary*.txt )
shopt -u nullglob

if [[ ${#busco_summaries[@]} -eq 0 ]]; then
  echo "WARNING: No BUSCO short_summary*.txt files found under ${BUSCO_DIR}"
  echo "         MultiQC will run, but BUSCO module may be empty."
else
  rm -f "${CLEAN_BUSCO_DIR}"/*.txt || true
  for f in "${busco_summaries[@]}"; do
    cp -v "${f}" "${CLEAN_BUSCO_DIR}/" >/dev/null
  done
fi

################################################################################
# Step 4: Run MultiQC
################################################################################
echo "📈 Running MultiQC..."
multiqc -f -o "${MULTIQC_OUTPUT}" "${STATS_DIR}" "${CLEAN_BUSCO_DIR}"

echo "✅ Done!"
echo "   Summary table : ${SUMMARY_FILE}"
echo "   MultiQC report: ${MULTIQC_OUTPUT}/multiqc_report.html"
echo "   Clean BUSCO   : ${CLEAN_BUSCO_DIR}"
