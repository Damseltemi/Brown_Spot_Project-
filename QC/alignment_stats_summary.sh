#!/bin/bash
set -euo pipefail

# Directory that contains your BAM + *.stats.txt files
ALIGN_DIR="/scratch/acicola_project/fastas/sydowiaspp/sydowiafastq/Alignmentsydonew2"

cd "$ALIGN_DIR"

SUMMARY="${ALIGN_DIR}/alignment_stats_summary.tsv1"

echo "Writing summary to: $SUMMARY"

# Header
echo -e "Sample\tTotal_reads\tMapped_reads\tMapped_pct\tProperly_paired_reads\tProperly_paired_pct" \
  > "$SUMMARY"

shopt -s nullglob
stat_files=( *.stats.txt )
shopt -u nullglob

if [[ ${#stat_files[@]} -eq 0 ]]; then
  echo "ERROR: No *.stats.txt files found in $ALIGN_DIR"
  exit 1
fi

for f in "${stat_files[@]}"; do
  sample="${f%.stats.txt}"

  vals=$(
    awk '
      / in total/ {
        total=$1
      }
      / mapped \(/ && $4=="mapped" {
        mapped=$1
        tmp=$5
        gsub("[()%]", "", tmp)   # remove ( and % and )
        mapped_pct=tmp
      }
      / properly paired \(/ {
        pp=$1
        tmp=$5
        gsub("[()%]", "", tmp)
        pp_pct=tmp
      }
      END {
        if (total=="")      total="NA"
        if (mapped=="")     mapped="NA"
        if (mapped_pct=="") mapped_pct="NA"
        if (pp=="")         pp="NA"
        if (pp_pct=="")     pp_pct="NA"
        OFS="\t"
        print total, mapped, mapped_pct, pp, pp_pct
      }
    ' "$f"
  )

  echo -e "${sample}\t${vals}" >> "$SUMMARY"
done

echo "✅ Done."
echo "   Summary table: $SUMMARY"
#!/bin/bash
set -euo pipefail

