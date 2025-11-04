#!/bin/bash

# runFeatureCounts.sh
#
# Description:
#   This script runs featureCounts to generate count matrices from BAM files
#   supporting different library strandedness types.
#
# Usage:
#   ./runFeatureCounts.sh -i input_dir -o output_dir -a annotation_file -s strandedness [-t threads] [-f feature_type]
#
# Parameters:
#   -i  Input directory containing BAM files (default: current directory)
#   -o  Output directory for count matrices (required)
#   -a  GTF annotation file (required)
#   -s  Strandedness: 0=unstranded, 1=stranded, 2=reversely stranded (default: 0)
#   -t  Number of threads for featureCounts (default: 8)
#   -f  Feature type to count in GTF (default: exon)
#   -p  Count read pairs instead of reads (default: yes, use empty to disable)
#   -g  GTF attribute to use for grouping (default: gene_id)
#
# If `-p yes``, then "-p --countReadPairs -B -C" are added to arguments list
#
# Output files:
#   - Raw featureCounts output: <output_dir>/raw/counts_matrix.txt
#   - Processed count matrix: <output_dir>/count_matrices/sorted_counts_matrix.txt
#   - Log files: <output_dir>/logs/featureCounts.log
#
# Examples:
#   # Unstranded library
#   ./runFeatureCounts.sh -i ./aligned_bams -o ./counts -a /path/to/annotation.gtf -s 0
#
#   # Reverse stranded library (common for Illumina TruSeq)
#   ./runFeatureCounts.sh -i ./aligned_bams -o ./counts -a /path/to/annotation.gtf -s 2
#
# Author: Anton Zhelonkin
# Date: 2025

# Default parameters
INPUT_DIR="."
THREADS=8
STRANDEDNESS=0
FEATURE_TYPE="exon"
PAIR_MODE="yes"
GTF_ATTRIBUTE="gene_id"

# Usage message function
usage() {
    echo "Usage: $0 -i input_dir -o output_dir -a annotation_file -s strandedness [-t threads] [-f feature_type]"
    echo "  -i  Input directory containing BAM files (default: current directory)"
    echo "  -o  Output directory for count matrices (required)"
    echo "  -a  GTF annotation file (required)"
    echo "  -s  Strandedness: 0=unstranded, 1=stranded, 2=reversely stranded (default: 0)"
    echo "  -t  Number of threads for featureCounts (default: 8)"
    echo "  -f  Feature type to count in GTF (default: exon)"
    echo "  -p  Count read pairs instead of reads (default: yes, use empty to disable)"
    echo "  -g  GTF attribute to use for grouping (default: gene_id)"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:a:s:t:f:p:g:" opt; do
    case ${opt} in
        i )
            INPUT_DIR=$OPTARG
            ;;
        o )
            OUTPUT_BASE_DIR=$OPTARG
            ;;
        a )
            ANNOTATION_FILE=$OPTARG
            ;;
        s )
            STRANDEDNESS=$OPTARG
            ;;
        t )
            THREADS=$OPTARG
            ;;
        f )
            FEATURE_TYPE=$OPTARG
            ;;
        p )
            PAIR_MODE=$OPTARG
            ;;
        g )
            GTF_ATTRIBUTE=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Validate required parameters
if [ -z "$OUTPUT_BASE_DIR" ] || [ -z "$ANNOTATION_FILE" ]; then
    echo "Error: Output directory (-o) and annotation file (-a) are required."
    usage
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist."
    exit 1
fi

# Check if annotation file exists
if [ ! -f "$ANNOTATION_FILE" ]; then
    echo "Error: Annotation file '$ANNOTATION_FILE' does not exist."
    exit 1
fi

# Validate strandedness parameter
if [[ "$STRANDEDNESS" != "0" && "$STRANDEDNESS" != "1" && "$STRANDEDNESS" != "2" ]]; then
    echo "Error: Strandedness (-s) must be 0, 1, or 2."
    usage
fi

# Set up strandedness description for logs
STRAND_DESC=""
case $STRANDEDNESS in
    0) STRAND_DESC="unstranded" ;;
    1) STRAND_DESC="stranded" ;;
    2) STRAND_DESC="reversely stranded" ;;
esac

# Define output directories
RAW_OUTPUT_DIR="${OUTPUT_BASE_DIR}/raw_fc_output"
MATRIX_OUTPUT_DIR="${OUTPUT_BASE_DIR}/count_matrices_fc"
LOG_DIR="${OUTPUT_BASE_DIR}/logs_fc"

# Create output directories
mkdir -p "$RAW_OUTPUT_DIR"
mkdir -p "$MATRIX_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Output files
COUNTS_FILE="$RAW_OUTPUT_DIR/counts_matrix.txt"
SORTED_COUNTS_FILE="$MATRIX_OUTPUT_DIR/sorted_counts_matrix.txt"
LOG_FILE="$LOG_DIR/featureCounts.log"

# Log parameters
echo "------------------------------------------------" > "$LOG_FILE"
echo "featureCounts run: $(date)" >> "$LOG_FILE"
echo "Input directory: $INPUT_DIR" >> "$LOG_FILE"
echo "Output directory: $OUTPUT_BASE_DIR" >> "$LOG_FILE"
echo "Annotation file: $ANNOTATION_FILE" >> "$LOG_FILE"
echo "Strandedness: $STRANDEDNESS ($STRAND_DESC)" >> "$LOG_FILE"
echo "Feature type: $FEATURE_TYPE" >> "$LOG_FILE"
echo "GTF attribute: $GTF_ATTRIBUTE" >> "$LOG_FILE"
echo "Threads: $THREADS" >> "$LOG_FILE"
echo "Pair mode: $PAIR_MODE" >> "$LOG_FILE"
echo "------------------------------------------------" >> "$LOG_FILE"

# Find all BAM files in the input directory (non-recursive) and sort them lexicographically
mapfile -t SORTED_BAM_FILES < <(ls "$INPUT_DIR"/*.bam | sort)

echo "Processing BAM files: ${SORTED_BAM_FILES[*]}" | tee -a "$LOG_FILE"

# Add debugging information
echo "Number of BAM files: ${#SORTED_BAM_FILES[@]}" | tee -a "$LOG_FILE"
for file in "${SORTED_BAM_FILES[@]}"; do
    echo "File: $file" >> "$LOG_FILE"
done

# Check if any BAM files were found
if [ ${#SORTED_BAM_FILES[@]} -eq 0 ]; then
    echo "No BAM files found in $INPUT_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

# Setup paired-end read parameters
PAIRED_ARGS=""
if [ "$PAIR_MODE" == "yes" ]; then
    PAIRED_ARGS="-p --countReadPairs -B -C"
fi

# Run featureCounts on all sorted BAM files in one command
echo "Starting featureCounts with strandedness=$STRANDEDNESS ($STRAND_DESC)" | tee -a "$LOG_FILE"

featureCounts \
    -a "$ANNOTATION_FILE" \
    -o "$COUNTS_FILE" \
    $PAIRED_ARGS \
    -s $STRANDEDNESS \
    -t $FEATURE_TYPE \
    -g $GTF_ATTRIBUTE \
    -T "$THREADS" \
    "${SORTED_BAM_FILES[@]}" \
    >> "$LOG_FILE" 2>&1

if [ $? -ne 0 ]; then
    echo "featureCounts failed. Check the log file at $LOG_FILE for details." | tee -a "$LOG_FILE"
    exit 1
fi

# Reorder columns in the output file
awk 'BEGIN {FS=OFS="\t"}
    NR==1 {next}  # Skip the first line (featureCounts summary)
    NR==2 {
        printf "Geneid";
        for (i=7; i<=NF; i++) {
            split($i, a, "/");  # Split the path
            split(a[length(a)], b, ".");  # Split the filename
            printf "\t%s", b[1];  # Print only the sample name
        }
        printf "\n";
        next;
    }
    {
        printf "%s", $1;
        for (i=7; i<=NF; i++) printf "\t%s", $i;
        printf "\n";
    }' "$COUNTS_FILE" > "$SORTED_COUNTS_FILE"

# Log completion
echo "featureCounts analysis completed with strandedness=$STRANDEDNESS ($STRAND_DESC)" | tee -a "$LOG_FILE"
echo "Results saved to $SORTED_COUNTS_FILE" | tee -a "$LOG_FILE"
echo "Completed processing ${#SORTED_BAM_FILES[@]} BAM files."
