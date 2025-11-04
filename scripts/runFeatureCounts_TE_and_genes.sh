#!/usr/bin/env bash
set -euo pipefail

# runFeatureCounts_TE_and_genes.sh
#
# What it does
#   1) TE counts (featureCounts on SAF; -M --fraction; fragments; default unstranded)
#   2) Gene counts (calls your existing runFeatureCounts.sh for GTF; fragments; per-project strandedness)
#   3) Optionally TE sense/antisense matrices (for stranded libs)
#   4) Outputs tidy count matrices and an optional combined matrix (genes + TE)
#
# Inputs (all absolute paths):
#   -i   BAM directory with *.bam (non-recursive)
#   -o   Output base directory (we’ll create subfolders)
#   -g   Gene GTF (e.g. gencode.vM37.primary_assembly.annotation.filtered.gtf)
#   -e   TE SAF (your /data1/shared/ref/.../GRCm39_rmsk_TE_formatted.saf)
#   -S   Gene strandedness for featureCounts: 0|1|2   (0=unstranded, 1=forward, 2=reverse)
#   -t   Threads (default 12)
#   --te-strand   One of: unstranded|sense_antisense  (default: unstranded)
#
# Notes
#   - TE is counted with:  -M --fraction  -p --countReadPairs -B -C
#   - TE default is unstranded (-s 0) to avoid halving signal from antisense TE transcription.
#   - If you set --te-strand sense_antisense AND your library is stranded (e.g., reverse), we also produce TE-sense and TE-antisense matrices by running -s 2 and -s 1 respectively.
#
# Author: Anton Zhelonkin
# Date: 2025-09-12

# Defaults
THREADS=12
TE_STRAND_MODE="unstranded"

usage() {
  cat >&2 <<EOF
Usage: $0 -i BAM_DIR -o OUT_DIR -g GENE_GTF -e TE_SAF -S GENE_STRAND [-t THREADS] [--te-strand unstranded|sense_antisense]

Examples:
  # 13441-DM (reverse-stranded genes; TE unstranded)
  $0 \\
    -i /data1/users/antonz/data/DM_summer_2025/13441-DM/results_mm39_TE/star_salmon/bam_fin \\
    -o /data1/users/antonz/data/DM_summer_2025/13441-DM/results_mm39_TE/star_salmon \\
    -g /data1/users/antonz/data/DM_summer_2025/13441-DM/results_mm39_TE/genome/gencode.vM37.primary_assembly.annotation.filtered.gtf \\
    -e /data1/shared/ref/mouse/Ensembl/mm39/GRCm39_rmsk_TE_formatted.saf \\
    -S 2

  # 13444-DM (unstranded genes; also produce TE sense/antisense if library were stranded)
  $0 -i .../13444-DM/.../bam_fin -o .../13444-DM/.../star_salmon -g .../gencode...gtf -e .../TE_formatted.saf -S 0 --te-strand unstranded
EOF
  exit 1
}

# Parse args
ARGS=()
while (( "$#" )); do
  case "$1" in
    -i) BAM_DIR="$2"; shift 2;;
    -o) OUT_BASE="$2"; shift 2;;
    -g) GENE_GTF="$2"; shift 2;;
    -e) TE_SAF="$2"; shift 2;;
    -S) GENE_STRAND="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --te-strand) TE_STRAND_MODE="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

# Validate
for v in BAM_DIR OUT_BASE GENE_GTF TE_SAF GENE_STRAND; do
  [[ -z "${!v:-}" ]] && echo "Missing required: $v" >&2 && usage
done
[[ -d "$BAM_DIR" ]] || { echo "BAM_DIR not found: $BAM_DIR" >&2; exit 1; }
[[ -f "$GENE_GTF" ]] || { echo "GENE_GTF not found: $GENE_GTF" >&2; exit 1; }
[[ -f "$TE_SAF"  ]] || { echo "TE_SAF not found:  $TE_SAF"  >&2; exit 1; }
[[ "$GENE_STRAND" =~ ^[012]$ ]] || { echo "-S must be 0,1,2" >&2; exit 1; }
[[ "$TE_STRAND_MODE" =~ ^(unstranded|sense_antisense)$ ]] || { echo "--te-strand must be unstranded|sense_antisense" >&2; exit 1; }

# Make output dirs
TE_DIR="${OUT_BASE}/featurecounts_TE"
GENE_DIR="${OUT_BASE}/fc_genes"
mkdir -p "$TE_DIR" "$GENE_DIR"

# Deterministic BAM list
mapfile -t BAM_LIST < <(ls "$BAM_DIR"/*.bam | sort)
if (( ${#BAM_LIST[@]} == 0 )); then
  echo "No BAMs in $BAM_DIR" >&2; exit 1
fi

echo "BAMs (${#BAM_LIST[@]}):"
printf '  %s\n' "${BAM_LIST[@]}"

# Helper: build a space-separated, shell-escaped list
join_bams() {
  local IFS=' '
  echo "${BAM_LIST[*]}"
}

# --------------------------- TE COUNTS (first) ---------------------------
echo "[TE] Counting TEs (fractional multi-mappers; fragments)..."

# Unstranded matrix
TE_RAW="${TE_DIR}/te_counts_raw.txt"
TE_MAT="${TE_DIR}/te_counts_matrix.txt"

featureCounts \
  -M \
  -F SAF -a "$TE_SAF" \
  -o "$TE_RAW" \
  -s 0 \
  -p --countReadPairs -B -C \
  -T "$THREADS" \
  $(join_bams)

# Post-process to a clean matrix with sample names
awk 'BEGIN{FS=OFS="\t"} 
  NR==1{next}
  NR==2{printf "Geneid"; for(i=7;i<=NF;i++){split($i,a,"/"); split(a[length(a)],b,"."); printf "\t" b[1]} printf "\n"; next}
  {printf $1; for(i=7;i<=NF;i++) printf "\t" $i; printf "\n"}' \
  "$TE_RAW" > "$TE_MAT"

echo "[TE] Unstranded matrix -> $TE_MAT"

# Optional: TE sense/antisense (only meaningful if library is stranded)
if [[ "$TE_STRAND_MODE" == "sense_antisense" ]]; then
  echo "[TE] Producing sense/antisense matrices (requires stranded lib)."

  # For reverse-stranded libs: sense = -s 2, antisense = -s 1
  # For forward-stranded libs:  sense = -s 1, antisense = -s 2
  # For unstranded libs: not possible (skip).
  if [[ "$GENE_STRAND" == "0" ]]; then
    echo "[TE] Library is unstranded (-S 0). Skipping sense/antisense." >&2
  else
    if [[ "$GENE_STRAND" == "2" ]]; then SENSE_S=2; ANTISENSE_S=1; else SENSE_S=1; ANTISENSE_S=2; fi

    TE_RAW_S="${TE_DIR}/te_counts_sense_raw.txt"
    TE_MAT_S="${TE_DIR}/te_counts_sense_matrix.txt"
    featureCounts -M --fraction -F SAF -a "$TE_SAF" -o "$TE_RAW_S" -s $SENSE_S -p --countReadPairs -B -C -T "$THREADS" $(join_bams)
    awk 'BEGIN{FS=OFS="\t"} NR==1{next} NR==2{printf "Geneid"; for(i=7;i<=NF;i++){split($i,a,"/"); split(a[length(a)],b,"."); printf "\t" b[1]} printf "\n"; next} {printf $1; for(i=7;i<=NF;i++) printf "\t" $i; printf "\n"}' "$TE_RAW_S" > "$TE_MAT_S"
    echo "[TE] Sense matrix -> $TE_MAT_S"

    TE_RAW_A="${TE_DIR}/te_counts_antisense_raw.txt"
    TE_MAT_A="${TE_DIR}/te_counts_antisense_matrix.txt"
    featureCounts -M --fraction -F SAF -a "$TE_SAF" -o "$TE_RAW_A" -s $ANTISENSE_S -p --countReadPairs -B -C -T "$THREADS" $(join_bams)
    awk 'BEGIN{FS=OFS="\t"} NR==1{next} NR==2{printf "Geneid"; for(i=7;i<=NF;i++){split($i,a,"/"); split(a[length(a)],b,"."); printf "\t" b[1]} printf "\n"; next} {printf $1; for(i=7;i<=NF;i++) printf "\t" $i; printf "\n"}' "$TE_RAW_A" > "$TE_MAT_A"
    echo "[TE] Antisense matrix -> $TE_MAT_A"
  fi
fi

# --------------------------- GENE COUNTS (second) ---------------------------
echo "[GENE] Counting genes via your runFeatureCounts.sh..."

RAW_GENE_DIR="${GENE_DIR}/raw_fc_output"
MAT_GENE_DIR="${GENE_DIR}/count_matrices_fc"
mkdir -p "$RAW_GENE_DIR" "$MAT_GENE_DIR"

# We reuse your script; it already sets:  -p --countReadPairs -B -C  when -p yes
# and does tidy post-processing to sorted_counts_matrix.txt
$(dirname "$0")/runFeatureCounts.sh \
  -i "$BAM_DIR" \
  -o "$GENE_DIR" \
  -a "$GENE_GTF" \
  -s "$GENE_STRAND" \
  -t "$THREADS" \
  -f exon \
  -g gene_id \
  -p yes

GENE_MAT="${MAT_GENE_DIR}/sorted_counts_matrix.txt"
echo "[GENE] Gene matrix -> $GENE_MAT"

# --------------------------- OPTIONAL: combine genes + TE ---------------------------
# Row-bind (headers identical order). We assume sample order is identical because both used sorted BAM lists.
COMBINED="${OUT_BASE}/combined_gene_TE_counts.tsv"
if [[ -f "$GENE_MAT" && -f "$TE_MAT" ]]; then
  { head -n 1 "$GENE_MAT"; tail -n +2 "$GENE_MAT"; tail -n +2 "$TE_MAT"; } > "$COMBINED"
  echo "[COMBINE] Combined matrix -> $COMBINED"
else
  echo "[COMBINE] Skipped (missing one of the matrices)." >&2
fi

echo "All done."
