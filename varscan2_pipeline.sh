#!/bin/bash

# =============================================================================
# VarScan2 Somatic Variant and CNV Analysis Pipeline
# Author      : Parameswar
# Description : Paired tumour/normal somatic SNP, INDEL, Germline, LOH, and
#               copy number variation analysis using VarScan2, samtools,
#               bam-readcount, and fpfilter.
#
# QUICK SETUP
# -----------
# 1. Set all paths and parameters in the CONFIGURATION section below.
# 2. Create your pairs file (see FILE FORMAT notes).
# 3. chmod +x varscan_pipeline.sh && ./varscan_pipeline.sh
#
# FILE FORMAT
# -----------
# FILE_PAIRS_LIST is a CSV with one tumour/normal pair per line.
# Column 1 = TUMOUR identifier, Column 2 = NORMAL identifier.
# Each identifier is the BAM filename as it appears in the pairs file.
#
# The sample name is derived by stripping PAIRS_SUFFIX from the identifier.
# The BAM is then located as: $BAM_DIR/<sample_name>$BAM_SUFFIX
#
# Example — current CRC cohort (sample_pairs.csv):
#   50T_CRC_final.bam,50N_CRC_final.bam
#   PAIRS_SUFFIX="_final.bam"
#   BAM_SUFFIX="_final.bam"
#   BAM_DIR="$PWD"
#
# IMPORTANT: Incomplete pairs (empty normal field) are skipped with a warning.
#   Verify all pairs are complete before running. Example of incomplete entry:
#   19T_CRC_final.bam,     ← no matched normal — will be skipped
# =============================================================================

echo "========================================================================"
echo "VarScan2 Somatic Variant and CNV Analysis Pipeline"
echo "========================================================================"

set -o pipefail

# =============================================================================
# CONFIGURATION — edit this section before running
# =============================================================================

# --- Paths -------------------------------------------------------------------
# GENOMEIDX1: the only path you must configure — absolute path to your GRCh38 reference FASTA.
# All other tools (VarScan jar, fpfilter.pl) are bundled in the repo under software/ and scripts/.
GENOMEIDX1=""                           # e.g. /data/ref/Homo_sapiens.GRCh38.dna.toplevel.fna
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOFTWAREDIR="$SCRIPT_DIR/software"     # directory containing VarScan.v2.3.9.jar
SCRIPTSDIR="$SCRIPT_DIR/scripts"       # directory containing fpfilter.pl
VARSCANDIR="$PWD"
SCRIPT_PATH="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

# --- BAM files ---------------------------------------------------------------
BAM_DIR="$PWD"                          # directory containing BAM files
BAM_SUFFIX="_final.bam"                 # suffix of BAM files

# --- Pairs file --------------------------------------------------------------
FILE_PAIRS_LIST="sample_pairs.csv"      # one tumour,normal pair per line (tumour col 1, normal col 2)
PAIRS_SUFFIX="_final.bam"              # strip this from each entry to get sample name
                                        # use ".flagstats" if pairs file lists flagstat filenames

# --- Target regions (WES) ----------------------------------------------------
# Path to capture kit BED file (e.g. Agilent SureSelect, Twist Exome).
# When set, mpileup is restricted to captured regions only, which:
#   - reduces mpileup file sizes by ~95% compared to whole-genome pileup
#   - excludes off-target reads that generate noise in CNV windows
#   - naturally excludes alt contigs, unplaced scaffolds, and decoy sequences
# Leave empty for WGS (full-genome pileup).
TARGET_BED=""

# --- VCF sample naming -------------------------------------------------------
# Optional: plain-text file with one sample name per line (normal first, tumour second).
# Provides correct column labels in VCF headers instead of generic NORMAL/TUMOR.
# Leave empty to use VarScan default generic names.
VCF_SAMPLE_LIST=""

# --- Somatic calling parameters ----------------------------------------------
MIN_COVERAGE=20                         # minimum total site depth
MIN_COVERAGE_NORMAL=10                  # minimum depth in normal sample
MIN_COVERAGE_TUMOR=20                   # minimum depth in tumour sample
                                        # WES: 20x tumour minimum ensures ≥2 supporting reads
                                        # at MIN_VAR_FREQ=0.10, reducing single-read false positives
MIN_BASE_QUAL=20                        # minimum base quality applied by VarScan and mpileup
                                        # (VarScan2 official default: 20)
MIN_VAR_FREQ=0.10                       # minimum variant allele frequency to call
                                        # (VarScan2 official somatic default: 0.10)
MIN_FREQ_FOR_HOM=0.75                   # minimum VAF to classify as homozygous
NORMAL_PURITY=1.0
TUMOR_PURITY=1.0
P_VALUE=0.99                            # general calling threshold (VarScan2 default: 0.99)
                                        # IMPORTANT: setting this value causes VarScan to
                                        # CALCULATE actual Fisher's p-values and accept
                                        # variants where p <= 0.99. If unset or set to 1.0,
                                        # VarScan skips p-value calculation entirely and
                                        # inserts a placeholder 0.98 — per VarScan FAQ.
                                        # Somatic specificity is controlled by SOMATIC_P_VALUE.
SOMATIC_P_VALUE=0.05                    # significance threshold for somatic classification
STRAND_FILTER=1                         # 1 = exclude variants supported on only one strand

# --- processSomatic parameters -----------------------------------------------
MIN_TUMOR_FREQ=0.10                     # min tumour VAF for high-confidence somatic call
MAX_NORMAL_FREQ=0.05                    # max normal VAF for high-confidence somatic call
PROCESS_P_VALUE=0.05                    # p-value cutoff for .hc classification

# --- CNV parameters ----------------------------------------------------------
CNV_MIN_COVERAGE=20                     # minimum depth per window for CNV analysis
CNV_P_VALUE=0.01                        # change-point significance threshold
                                        # (VarScan2 official default: 0.01)
MIN_SEGMENT_SIZE=10                     # minimum windows per CNV segment
MAX_SEGMENT_SIZE=100                    # maximum windows per CNV segment

# --- copyCaller thresholds ---------------------------------------------------
CNV_AMP_THRESHOLD=0.25                  # log2 ratio lower bound for amplification call
                                        # (VarScan2 official default: 0.25)
CNV_DEL_THRESHOLD=0.25                  # log2 ratio upper bound for deletion call
                                        # (VarScan2 official default: 0.25)
CNV_RECENTER_UP=0                       # shift log2 baseline up by this amount before calling
CNV_RECENTER_DOWN=0                     # shift log2 baseline down by this amount before calling
                                        # IMPORTANT: inspect *.copynumber output before running
                                        # copyCaller. In tumours with widespread chromosomal
                                        # instability, the median log2 ratio is offset from 0.
                                        # Set RECENTER_UP or RECENTER_DOWN to correct this
                                        # offset; incorrect baseline produces systematically
                                        # wrong amplification/deletion classifications.

# --- bam-readcount parameters ------------------------------------------------
BRC_MAP_QUAL=10                         # minimum mapping quality
                                        # VarScan FAQ: "A minimum mapping quality of 10 is
                                        # even better" to exclude ambiguously mapped reads.
                                        # MAPQ 1-9 reads have low placement confidence and
                                        # introduce phantom variant signals in segmental
                                        # duplications and paralogous gene regions.
BRC_BASE_QUAL=15                        # minimum base quality (fpfilter docs recommend 15)

# --- Performance -------------------------------------------------------------
MAX_PARALLEL_JOBS=30

# --- Batch processing controls -----------------------------------------------
# These defaults are overridden at runtime by command-line arguments.
# Use --from, --to, --stage, --resume, and --dry-run to control execution.
FROM_STAGE=1        # first stage to run  (1–10)
TO_STAGE=10         # last  stage to run  (1–10)
RESUME=false        # skip stages when SHA256 resume marker and outputs both match
DRY_RUN=false       # print what would run without executing anything

# --- Resume state -------------------------------------------------------------
RESUME_STATEDIR="$VARSCANDIR/.resume_state"

# --- Internal job pool state --------------------------------------------------
JOB_LIMIT=1
JOB_PIDS=()

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

create_directory() {
    local dir_path="$1"
    if [ ! -d "$dir_path" ]; then
        mkdir -p "$dir_path"
        chmod 755 "$dir_path"
        log_message "Created: $dir_path"
    fi
}

check_file_exists() {
    if [ ! -f "$1" ]; then
        log_message "ERROR: Required file not found: $1"
        exit 1
    fi
}

check_command_exists() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log_message "ERROR: Required command not found on PATH: $1"
        exit 1
    fi
}

# Extract sample name from a pairs-file entry by stripping PAIRS_SUFFIX.
get_sample_name() {
    echo "${1%"$PAIRS_SUFFIX"}"
}

# Return the full path to a sample's BAM file.
get_bam_path() {
    echo "$BAM_DIR/$(get_sample_name "$1")${BAM_SUFFIX}"
}

# Verify that a BAM index exists alongside the BAM file.
# bam-readcount requires an index; missing index causes a silent empty output.
# Accepts both <file>.bam.bai and <file>.bai naming conventions.
check_bam_index() {
    local bam="$1"
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ]; then
        log_message "ERROR: BAM index not found for: $bam"
        log_message "  Expected: ${bam}.bai  or  ${bam%.bam}.bai"
        log_message "  Run: samtools index $bam"
        exit 1
    fi
}

# Return 0 (true) if a stage has expected output files on disk.
# Used together with SHA256 resume markers to skip completed stages.
stage_output_exists() {
    case "$1" in
        1)  find "$FLAGSTATDIR"   -name '*.flagstats'         -type f 2>/dev/null | grep -q . ;;
        2)  find "$MPILEUPDIR"    -name '*.mpileup'           -type f 2>/dev/null | grep -q . ;;
        3)  find "$SOMATICDIR"    -name '*.snp.vcf'           -type f 2>/dev/null | grep -q . ;;
        4)  find "$SOMATICDIR"    -name '*.hc.vcf'            -type f 2>/dev/null | grep -q . ;;
        5)  find "$COPYNUMBERDIR" -name '*.copynumber'        -type f 2>/dev/null | grep -q . ;;
        6)  find "$COPYNUMBERDIR" -name '*.copynumber.called' -type f 2>/dev/null | grep -q . ;;
        7)  find "$SNPVARDIR"     -name '*.var'               -type f 2>/dev/null | grep -q . ;;
        8)  find "$READCOUNTDIR"  -name '*.readcount'         -type f 2>/dev/null | grep -q . ;;
        9)  find "$FILTEREDDIR"   -name '*.fpfilter.vcf'      -type f 2>/dev/null | grep -q . ;;
        10) [ -f "$VARSCANDIR/varscan_pipeline_summary.txt" ] ;;
        *)  return 1 ;;
    esac
}

sha256_of_string() {
    printf '%s' "$1" | sha256sum | awk '{print $1}'
}

sha256_of_file_or_missing() {
    local f="$1"
    if [ -f "$f" ]; then
        sha256sum "$f" | awk '{print $1}'
    else
        echo "MISSING"
    fi
}

hash_manifest() {
    local dir="$1"
    local pattern="$2"
    if [ ! -d "$dir" ]; then
        echo "MISSING_DIR"
        return 0
    fi
    local listing
    listing=$(find "$dir" -maxdepth 1 -type f -name "$pattern" -printf '%f|%s|%T@\n' 2>/dev/null | LC_ALL=C sort)
    if [ -z "$listing" ]; then
        echo "EMPTY"
    else
        sha256_of_string "$listing"
    fi
}

hash_pairs_cleaned() {
    if [ ! -f "$FILE_PAIRS_LIST" ]; then
        echo "MISSING"
        return 0
    fi
    local normalized
    normalized=$(while IFS=',' read -r entry1 entry2; do
        entry1=$(echo "$entry1" | tr -d ' \r')
        entry2=$(echo "$entry2" | tr -d ' \r')
        validate_pair "$entry1" "$entry2" || continue
        printf '%s,%s\n' "$entry1" "$entry2"
    done < "$FILE_PAIRS_LIST")
    if [ -z "$normalized" ]; then
        echo "EMPTY"
    else
        sha256_of_string "$normalized"
    fi
}

hash_paired_bam_metadata() {
    if [ ! -f "$FILE_PAIRS_LIST" ]; then
        echo "MISSING"
        return 0
    fi

    local manifest=""
    while IFS=',' read -r entry1 entry2; do
        entry1=$(echo "$entry1" | tr -d ' \r')
        entry2=$(echo "$entry2" | tr -d ' \r')
        validate_pair "$entry1" "$entry2" || continue

        local tumor_bam normal_bam
        tumor_bam=$(get_bam_path "$entry1")
        normal_bam=$(get_bam_path "$entry2")

        for bam in "$normal_bam" "$tumor_bam"; do
            if [ -f "$bam" ]; then
                manifest+=$(stat -c '%n|%s|%Y' "$bam")
            else
                manifest+="${bam}|MISSING"
            fi
            manifest+=$'\n'
        done
    done < "$FILE_PAIRS_LIST"

    if [ -z "$manifest" ]; then
        echo "EMPTY"
    else
        sha256_of_string "$manifest"
    fi
}

compute_stage_hash() {
    local stage_num="$1"
    local payload=""

    payload+="script_sha256=$(sha256_of_file_or_missing "$SCRIPT_PATH")"$'\n'
    payload+="pairs_sha256=$(hash_pairs_cleaned)"$'\n'
    payload+="bam_sha256=$(hash_paired_bam_metadata)"$'\n'
    payload+="stage=${stage_num}"$'\n'

    case "$stage_num" in
        1)
            payload+="bam_dir=${BAM_DIR}"$'\n'
            payload+="bam_suffix=${BAM_SUFFIX}"$'\n'
            ;;
        2)
            payload+="reference_sha256=$(sha256_of_file_or_missing "$GENOMEIDX1")"$'\n'
            payload+="target_bed_sha256=$(sha256_of_file_or_missing "$TARGET_BED")"$'\n'
            payload+="BRC_MAP_QUAL=${BRC_MAP_QUAL}"$'\n'
            payload+="MIN_BASE_QUAL=${MIN_BASE_QUAL}"$'\n'
            ;;
        3)
            payload+="jar_sha256=$(sha256_of_file_or_missing "$SOFTWAREDIR/VarScan.v2.3.9.jar")"$'\n'
            payload+="mpileup_manifest=$(hash_manifest "$MPILEUPDIR" '*.mpileup')"$'\n'
            payload+="VCF_SAMPLE_LIST_SHA256=$(sha256_of_file_or_missing "$VCF_SAMPLE_LIST")"$'\n'
            payload+="MIN_COVERAGE=${MIN_COVERAGE}"$'\n'
            payload+="MIN_COVERAGE_NORMAL=${MIN_COVERAGE_NORMAL}"$'\n'
            payload+="MIN_COVERAGE_TUMOR=${MIN_COVERAGE_TUMOR}"$'\n'
            payload+="MIN_VAR_FREQ=${MIN_VAR_FREQ}"$'\n'
            payload+="MIN_FREQ_FOR_HOM=${MIN_FREQ_FOR_HOM}"$'\n'
            payload+="NORMAL_PURITY=${NORMAL_PURITY}"$'\n'
            payload+="TUMOR_PURITY=${TUMOR_PURITY}"$'\n'
            payload+="P_VALUE=${P_VALUE}"$'\n'
            payload+="SOMATIC_P_VALUE=${SOMATIC_P_VALUE}"$'\n'
            payload+="STRAND_FILTER=${STRAND_FILTER}"$'\n'
            payload+="BRC_MAP_QUAL=${BRC_MAP_QUAL}"$'\n'
            payload+="MIN_BASE_QUAL=${MIN_BASE_QUAL}"$'\n'
            ;;
        4)
            payload+="jar_sha256=$(sha256_of_file_or_missing "$SOFTWAREDIR/VarScan.v2.3.9.jar")"$'\n'
            payload+="somatic_raw_manifest=$(hash_manifest "$SOMATICDIR" '*.snp.vcf')"$'\n'
            payload+="somatic_raw_indel_manifest=$(hash_manifest "$SOMATICDIR" '*.indel.vcf')"$'\n'
            payload+="MIN_TUMOR_FREQ=${MIN_TUMOR_FREQ}"$'\n'
            payload+="MAX_NORMAL_FREQ=${MAX_NORMAL_FREQ}"$'\n'
            payload+="PROCESS_P_VALUE=${PROCESS_P_VALUE}"$'\n'
            ;;
        5)
            payload+="jar_sha256=$(sha256_of_file_or_missing "$SOFTWAREDIR/VarScan.v2.3.9.jar")"$'\n'
            payload+="mpileup_manifest=$(hash_manifest "$MPILEUPDIR" '*.mpileup')"$'\n'
            payload+="flagstats_manifest=$(hash_manifest "$FLAGSTATDIR" '*.flagstats')"$'\n'
            payload+="CNV_MIN_COVERAGE=${CNV_MIN_COVERAGE}"$'\n'
            payload+="MIN_BASE_QUAL=${MIN_BASE_QUAL}"$'\n'
            payload+="BRC_MAP_QUAL=${BRC_MAP_QUAL}"$'\n'
            payload+="MIN_SEGMENT_SIZE=${MIN_SEGMENT_SIZE}"$'\n'
            payload+="MAX_SEGMENT_SIZE=${MAX_SEGMENT_SIZE}"$'\n'
            payload+="CNV_P_VALUE=${CNV_P_VALUE}"$'\n'
            ;;
        6)
            payload+="jar_sha256=$(sha256_of_file_or_missing "$SOFTWAREDIR/VarScan.v2.3.9.jar")"$'\n'
            payload+="cnv_manifest=$(hash_manifest "$COPYNUMBERDIR" '*.copynumber')"$'\n'
            payload+="CNV_AMP_THRESHOLD=${CNV_AMP_THRESHOLD}"$'\n'
            payload+="CNV_DEL_THRESHOLD=${CNV_DEL_THRESHOLD}"$'\n'
            payload+="CNV_RECENTER_UP=${CNV_RECENTER_UP}"$'\n'
            payload+="CNV_RECENTER_DOWN=${CNV_RECENTER_DOWN}"$'\n'
            ;;
        7)
            payload+="hc_manifest=$(hash_manifest "$SOMATICDIR" '*.hc.vcf')"$'\n'
            ;;
        8)
            payload+="reference_sha256=$(sha256_of_file_or_missing "$GENOMEIDX1")"$'\n'
            payload+="snp_var_manifest=$(hash_manifest "$SNPVARDIR" '*.var')"$'\n'
            payload+="indel_var_manifest=$(hash_manifest "$INDELVARDIR" '*.var')"$'\n'
            payload+="BRC_MAP_QUAL=${BRC_MAP_QUAL}"$'\n'
            payload+="BRC_BASE_QUAL=${BRC_BASE_QUAL}"$'\n'
            ;;
        9)
            payload+="fpfilter_sha256=$(sha256_of_file_or_missing "$SCRIPTSDIR/fpfilter.pl")"$'\n'
            payload+="readcount_manifest=$(hash_manifest "$READCOUNTDIR" '*.readcount')"$'\n'
            payload+="somatic_manifest=$(hash_manifest "$SOMATICDIR" '*.vcf')"$'\n'
            payload+="MIN_VAR_FREQ=${MIN_VAR_FREQ}"$'\n'
            payload+="BRC_BASE_QUAL=${BRC_BASE_QUAL}"$'\n'
            ;;
        10)
            payload+="somatic_manifest=$(hash_manifest "$SOMATICDIR" '*.vcf')"$'\n'
            payload+="copy_manifest=$(hash_manifest "$COPYNUMBERDIR" '*.copynumber*')"$'\n'
            payload+="filtered_manifest=$(hash_manifest "$FILTEREDDIR" '*.fpfilter.vcf')"$'\n'
            ;;
    esac

    sha256_of_string "$payload"
}

stage_marker_path() {
    echo "$RESUME_STATEDIR/stage_$1.sha256"
}

write_stage_marker() {
    local stage_num="$1"
    local marker
    marker=$(stage_marker_path "$stage_num")
    compute_stage_hash "$stage_num" > "$marker"
}

resume_stage_match() {
    local stage_num="$1"
    local marker
    marker=$(stage_marker_path "$stage_num")
    if [ ! -f "$marker" ]; then
        return 1
    fi
    local expected current
    expected=$(cat "$marker")
    current=$(compute_stage_hash "$stage_num")
    [ "$expected" = "$current" ] && stage_output_exists "$stage_num"
}

job_pool_init() {
    JOB_LIMIT="$1"
    JOB_PIDS=()
}

job_pool_wait_all() {
    local failed=0
    local pid
    for pid in "${JOB_PIDS[@]}"; do
        if ! wait "$pid"; then
            failed=1
        fi
    done
    JOB_PIDS=()
    [ "$failed" -eq 0 ]
}

job_pool_add() {
    local pid="$1"
    JOB_PIDS+=("$pid")
    if [ "${#JOB_PIDS[@]}" -ge "$JOB_LIMIT" ]; then
        job_pool_wait_all
    fi
}

# Execute a pipeline stage with range, resume, and dry-run awareness.
# Usage: run_stage <stage_number> <function_name> "<description>"
run_stage() {
    local stage_num="$1"
    local stage_fn="$2"
    local stage_desc="$3"

    if [ "$stage_num" -lt "$FROM_STAGE" ] || [ "$stage_num" -gt "$TO_STAGE" ]; then
        log_message "--- Stage ${stage_num} skipped (outside range ${FROM_STAGE}–${TO_STAGE})"
        return 0
    fi

    if [ "$RESUME" = true ] && resume_stage_match "$stage_num"; then
        log_message "--- Stage ${stage_num} skipped — SHA256 resume match: ${stage_desc}"
        return 0
    fi

    if [ "$DRY_RUN" = true ]; then
        log_message "--- Stage ${stage_num} [dry-run] would run: ${stage_desc}"
        return 0
    fi

    if ! "$stage_fn"; then
        log_message "ERROR: Stage ${stage_num} failed: ${stage_desc}"
        return 1
    fi

    write_stage_marker "$stage_num"
}

# Skip a pair if either entry is empty (e.g. incomplete rows in sample_pairs.csv).
# Usage: validate_pair "$entry1" "$entry2" || continue
validate_pair() {
    if [ -z "$1" ] || [ -z "$2" ]; then
        log_message "WARNING: Incomplete pair (empty field) — skipping: '${1}','${2}'"
        return 1
    fi
    return 0
}

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

setup_directories() {
    log_message "Setting up output directory structure..."

    FLAGSTATDIR="$VARSCANDIR/flagstats"
    MPILEUPDIR="$VARSCANDIR/mpileup"
    SOMATICDIR="$VARSCANDIR/somatic"
    COPYNUMBERDIR="$VARSCANDIR/copynumber"
    SNPVARDIR="$VARSCANDIR/snp-VAR"
    INDELVARDIR="$VARSCANDIR/indel-VAR"
    READCOUNTDIR="$VARSCANDIR/readcount"
    FILTEREDDIR="$VARSCANDIR/filtered"

    for dir in "$FLAGSTATDIR" "$MPILEUPDIR" "$SOMATICDIR" \
               "$COPYNUMBERDIR" "$SNPVARDIR" "$INDELVARDIR" \
               "$READCOUNTDIR" "$FILTEREDDIR" "$RESUME_STATEDIR"; do
        create_directory "$dir"
    done

    log_message "Directory setup complete"
}

# =============================================================================
# STAGE 1: FLAGSTATS
# Generates alignment statistics and mapped read counts used to calculate the
# data-ratio for CNV depth normalisation in Stage 5.
# =============================================================================

generate_flagstats() {
    log_message "=== STAGE 1: Generating BAM flagstats ==="

    if [ ! -d "$BAM_DIR" ]; then
        log_message "ERROR: BAM directory not found: $BAM_DIR"
        exit 1
    fi

    job_pool_init "$MAX_PARALLEL_JOBS"
    for bam_file in "$BAM_DIR"/*"${BAM_SUFFIX}"; do
        [ -f "$bam_file" ] || continue
        local sample
        sample=$(basename "$bam_file" "${BAM_SUFFIX}")
        log_message "Flagstats: $sample"
        (samtools flagstat "$bam_file" > "$FLAGSTATDIR/${sample}.flagstats") &
        job_pool_add "$!" || return 1
    done
    job_pool_wait_all || return 1
    log_message "Flagstats complete"
}

# =============================================================================
# STAGE 2: PAIRED MPILEUP
# Generates a two-sample pileup file required by VarScan somatic and copynumber.
#
# CRITICAL: normal BAM must be listed FIRST, tumour BAM SECOND.
# VarScan somatic assigns Somatic/Germline/LOH status based on column order.
# Reversing the order silently swaps all classifications.
#
# -B        disables BAQ recalculation — VarScan docs and FAQ both recommend this
#           as BAQ is "occasionally too stringent for variant calling".
# -q        excludes reads with mapping quality below BRC_MAP_QUAL (min 10).
# -Q        excludes bases with base quality below MIN_BASE_QUAL (20), matching
#           the quality floor applied by VarScan — keeps read sets consistent
#           between mpileup and bam-readcount.
# -F 0x400  excludes duplicate-flagged reads (PCR/optical duplicates).
#           samtools mpileup does NOT exclude duplicates by default. BAMs from
#           GATK/Picard MarkDuplicates retain flagged duplicates in the file;
#           without this flag they contribute to depth and variant counts,
#           undermining deduplication.
# -l        (WES) restricts pileup to captured regions in TARGET_BED.
#           Reduces mpileup file size by ~95%, excludes off-target noise, and
#           naturally excludes alt contigs and decoy sequences. Leave TARGET_BED
#           empty for WGS.
#
# Output file naming: <normal>_<tumour>.mpileup
# This exact naming is used by all downstream stages — do not change it.
# =============================================================================

generate_mpileup() {
    log_message "=== STAGE 2: Generating paired mpileup (normal first, tumour second) ==="

    check_file_exists "$FILE_PAIRS_LIST"
    check_file_exists "$GENOMEIDX1"
    [ -n "$TARGET_BED" ] && check_file_exists "$TARGET_BED"

    # Build optional target-regions argument
    local target_opt=()
    [ -n "$TARGET_BED" ] && target_opt=(-l "$TARGET_BED")

    job_pool_init 10
    while IFS=',' read -r entry1 entry2; do
        entry1=$(echo "$entry1" | tr -d ' \r')
        entry2=$(echo "$entry2" | tr -d ' \r')
        validate_pair "$entry1" "$entry2" || continue

        local samplen samplet normal_bam tumor_bam out_mpileup
        samplet=$(get_sample_name "$entry1")    # column 1 = tumour
        samplen=$(get_sample_name "$entry2")    # column 2 = normal
        tumor_bam=$(get_bam_path "$entry1")
        normal_bam=$(get_bam_path "$entry2")
        out_mpileup="$MPILEUPDIR/${samplen}_${samplet}.mpileup"

        check_file_exists "$normal_bam"
        check_file_exists "$tumor_bam"

        log_message "Mpileup: Normal=${samplen}  Tumour=${samplet}"

        (samtools mpileup \
            -B \
            -q "$BRC_MAP_QUAL" \
            -Q "$MIN_BASE_QUAL" \
            -F 0x400 \
            "${target_opt[@]}" \
            -f "$GENOMEIDX1" \
            "$normal_bam" \
            "$tumor_bam" > "$out_mpileup") &

        job_pool_add "$!" || return 1
    done < "$FILE_PAIRS_LIST"
    job_pool_wait_all || return 1
    log_message "Mpileup complete"
}

# =============================================================================
# STAGE 3: VARSCAN SOMATIC VARIANT CALLING
# Calls SNPs and INDELs from the paired mpileup.
# Output prefix <normal>_<tumour> — VarScan auto-appends .snp.vcf and .indel.vcf
#
# --min-map-qual and --min-base-qual are applied by VarScan internally to the
# mpileup read data; explicit values ensure reproducibility across versions.
#
# --vcf-sample-list (optional): provides correct NORMAL/TUMOR sample names in
# VCF headers. Without it, VarScan uses generic column labels.
# =============================================================================

run_varscan_somatic() {
    log_message "=== STAGE 3: VarScan somatic variant calling ==="

    check_file_exists "$FILE_PAIRS_LIST"
    check_file_exists "$SOFTWAREDIR/VarScan.v2.3.9.jar"

    # Build optional vcf-sample-list argument
    local vcf_sample_opt=()
    if [ -n "$VCF_SAMPLE_LIST" ] && [ -f "$VCF_SAMPLE_LIST" ]; then
        vcf_sample_opt=(--vcf-sample-list "$VCF_SAMPLE_LIST")
        log_message "Using VCF sample list: $VCF_SAMPLE_LIST"
    fi

    job_pool_init "$MAX_PARALLEL_JOBS"
    while IFS=',' read -r entry1 entry2; do
        entry1=$(echo "$entry1" | tr -d ' \r')
        entry2=$(echo "$entry2" | tr -d ' \r')
        validate_pair "$entry1" "$entry2" || continue

        local samplen samplet mpileup_file
        samplet=$(get_sample_name "$entry1")    # column 1 = tumour
        samplen=$(get_sample_name "$entry2")    # column 2 = normal
        mpileup_file="$MPILEUPDIR/${samplen}_${samplet}.mpileup"

        if [ ! -f "$mpileup_file" ]; then
            log_message "WARNING: Mpileup not found: $mpileup_file — skipping pair"
            continue
        fi

        log_message "VarScan somatic: Normal=${samplen}  Tumour=${samplet}"

        # Output prefix — VarScan appends .snp.vcf and .indel.vcf automatically
        (java -jar "$SOFTWAREDIR/VarScan.v2.3.9.jar" somatic \
            "$mpileup_file" \
            "$SOMATICDIR/${samplen}_${samplet}" \
            --mpileup 1 \
            --min-coverage          "$MIN_COVERAGE" \
            --min-coverage-normal   "$MIN_COVERAGE_NORMAL" \
            --min-coverage-tumor    "$MIN_COVERAGE_TUMOR" \
            --min-var-freq          "$MIN_VAR_FREQ" \
            --min-freq-for-hom      "$MIN_FREQ_FOR_HOM" \
            --normal-purity         "$NORMAL_PURITY" \
            --tumor-purity          "$TUMOR_PURITY" \
            --p-value               "$P_VALUE" \
            --somatic-p-value       "$SOMATIC_P_VALUE" \
            --strand-filter         "$STRAND_FILTER" \
            --min-map-qual          "$BRC_MAP_QUAL" \
            --min-base-qual         "$MIN_BASE_QUAL" \
            --output-vcf 1 \
            "${vcf_sample_opt[@]}") &

        job_pool_add "$!" || return 1
    done < "$FILE_PAIRS_LIST"
    job_pool_wait_all || return 1
    log_message "VarScan somatic calling complete"
}

# =============================================================================
# STAGE 4: PROCESS SOMATIC
# Classifies raw variants from Stage 3 into Somatic, Germline, and LOH
# categories and writes separate high-confidence (.hc) VCF files.
#
# Outputs (per pair, for both .snp and .indel):
#   *.Somatic.vcf / *.Somatic.hc.vcf
#   *.Germline.vcf / *.Germline.hc.vcf
#   *.LOH.vcf / *.LOH.hc.vcf
#
# Must run AFTER Stage 3. Stages 7 and 8 depend on its .hc.vcf output.
# processSomatic inherits the VCF format from its input — .hc.vcf output
# requires that Stage 3 was run with --output-vcf 1.
# =============================================================================

process_somatic_variants() {
    log_message "=== STAGE 4: processSomatic — classifying into Somatic / Germline / LOH ==="

    job_pool_init "$MAX_PARALLEL_JOBS"

    for snp_file in "$SOMATICDIR"/*.snp.vcf; do
        [ -f "$snp_file" ] || continue
        log_message "processSomatic SNP: $(basename "$snp_file")"
        (java -jar "$SOFTWAREDIR/VarScan.v2.3.9.jar" processSomatic \
            "$snp_file" \
            --min-tumor-freq  "$MIN_TUMOR_FREQ" \
            --max-normal-freq "$MAX_NORMAL_FREQ" \
            --p-value         "$PROCESS_P_VALUE") &
        job_pool_add "$!" || return 1
    done

    for indel_file in "$SOMATICDIR"/*.indel.vcf; do
        [ -f "$indel_file" ] || continue
        log_message "processSomatic INDEL: $(basename "$indel_file")"
        (java -jar "$SOFTWAREDIR/VarScan.v2.3.9.jar" processSomatic \
            "$indel_file" \
            --min-tumor-freq  "$MIN_TUMOR_FREQ" \
            --max-normal-freq "$MAX_NORMAL_FREQ" \
            --p-value         "$PROCESS_P_VALUE") &
        job_pool_add "$!" || return 1
    done
    job_pool_wait_all || return 1
    log_message "processSomatic complete"
}

# =============================================================================
# STAGE 4 OUTPUT GUARD
# Verifies that processSomatic produced .hc.vcf files before Stages 7-9
# attempt to consume them. Aborts with a diagnostic message if none found.
# Common cause of failure: Stage 3 ran without --output-vcf 1, so processSomatic
# produced tab-separated .hc files instead of .hc.vcf files.
# =============================================================================

check_stage4_output() {
    local hc_count
    hc_count=$(find "$SOMATICDIR" -name '*.hc.vcf' -type f 2>/dev/null | wc -l)

    if [ "$hc_count" -eq 0 ]; then
        log_message "ERROR: processSomatic (Stage 4) produced no .hc.vcf files."
        log_message "  Expected files matching: $SOMATICDIR/*.hc.vcf"
        log_message "  Possible causes:"
        log_message "    1. Stage 3 did not produce any variant calls (check coverage/parameters)"
        log_message "    2. Stage 3 ran without --output-vcf 1 (processSomatic outputs .hc not .hc.vcf)"
        log_message "    3. processSomatic filtered all variants (try lowering MIN_TUMOR_FREQ)"
        exit 1
    fi

    log_message "Stage 4 guard passed: $hc_count .hc.vcf file(s) found in $SOMATICDIR"
}

# =============================================================================
# STAGE 5: VARSCAN COPY NUMBER ANALYSIS
# Computes per-window log2 read-depth ratios between tumour and normal.
#
# data-ratio = normal_mapped_reads / tumour_mapped_reads
# Normalises for sequencing depth differences between the two samples.
# Derived from flagstats with scale=6 decimal places to prevent baseline bias.
#
# --min-base-qual and --min-map-qual are applied by VarScan to the mpileup
# at each window; explicit values document the quality floor and ensure
# reproducibility. VarScan2 official defaults: both 20.
# =============================================================================

run_varscan_copynumber() {
    log_message "=== STAGE 5: VarScan copy number analysis ==="

    check_file_exists "$FILE_PAIRS_LIST"

    job_pool_init "$MAX_PARALLEL_JOBS"
    while IFS=',' read -r entry1 entry2; do
        entry1=$(echo "$entry1" | tr -d ' \r')
        entry2=$(echo "$entry2" | tr -d ' \r')
        validate_pair "$entry1" "$entry2" || continue

        local samplen samplet
        samplet=$(get_sample_name "$entry1")    # column 1 = tumour
        samplen=$(get_sample_name "$entry2")    # column 2 = normal

        # Calculate data-ratio from flagstats (scale=6 for sufficient precision)
        local dataratio="1.0"
        local normal_flags="$FLAGSTATDIR/${samplen}.flagstats"
        local tumor_flags="$FLAGSTATDIR/${samplet}.flagstats"

        if [ -f "$normal_flags" ] && [ -f "$tumor_flags" ]; then
            local n_mapped t_mapped
            # Use primary mapped counts to exclude secondary and supplementary
            # alignments (split reads, chimeric pairs) which inflate total mapped
            # counts disproportionately in CIN tumours with structural rearrangements.
            n_mapped=$(grep -m1 "primary mapped" "$normal_flags" | awk '{print $1}')
            t_mapped=$(grep -m1 "primary mapped" "$tumor_flags"  | awk '{print $1}')
            if [[ "$n_mapped" =~ ^[0-9]+$ ]] && [[ "$t_mapped" =~ ^[0-9]+$ ]] && [ "$t_mapped" -ne 0 ]; then
                dataratio=$(echo "scale=6; $n_mapped / $t_mapped" | bc)
            else
                log_message "WARNING: Could not parse primary mapped read counts for ${samplen}/${samplet} — using ratio=1.0"
            fi
        else
            log_message "WARNING: Flagstats missing for ${samplen} or ${samplet} — using ratio=1.0"
        fi

        log_message "Data ratio ${samplen}/${samplet}: $dataratio"

        local mpileup_file="$MPILEUPDIR/${samplen}_${samplet}.mpileup"
        if [ ! -f "$mpileup_file" ]; then
            log_message "WARNING: Mpileup not found: $mpileup_file — skipping pair"
            continue
        fi

        (java -jar "$SOFTWAREDIR/VarScan.v2.3.9.jar" copynumber \
            "$mpileup_file" \
            "$COPYNUMBERDIR/${samplen}_${samplet}" \
            --mpileup 1 \
            --min-coverage     "$CNV_MIN_COVERAGE" \
            --min-base-qual    "$MIN_BASE_QUAL" \
            --min-map-qual     "$BRC_MAP_QUAL" \
            --min-segment-size "$MIN_SEGMENT_SIZE" \
            --max-segment-size "$MAX_SEGMENT_SIZE" \
            --p-value          "$CNV_P_VALUE" \
            --data-ratio       "$dataratio") &

        job_pool_add "$!" || return 1
    done < "$FILE_PAIRS_LIST"
    job_pool_wait_all || return 1
    log_message "VarScan copy number complete"
}

# =============================================================================
# STAGE 6: COPY CALLER
# Segments raw copy number windows into discrete CNV calls.
#
# --amp-threshold / --del-threshold: log2 ratio thresholds for classifying a
#   segment as amplification or deletion (VarScan2 defaults: 0.25 each).
#   Corresponds to roughly one copy change in a 100% pure diploid tumour.
#   Lower for impure samples; raise for high-purity samples.
#
# --recenter-up / --recenter-down: correct the log2 baseline when the copy
#   number profile is globally shifted due to widespread chromosomal instability.
#   Inspect *.copynumber files and set to the median/modal log2 ratio before
#   running this stage. Incorrect baseline = systematically wrong CNV calls.
#
# --output-homdel-file: separate output for homozygous deletions (added v2.2.12).
# =============================================================================

run_copy_caller() {
    log_message "=== STAGE 6: VarScan copyCaller ==="

    job_pool_init "$MAX_PARALLEL_JOBS"
    for cnv_file in "$COPYNUMBERDIR"/*.copynumber; do
        [ -f "$cnv_file" ] || continue
        local base
        base=$(basename "$cnv_file" .copynumber)
        log_message "copyCaller: $base"
        (java -jar "$SOFTWAREDIR/VarScan.v2.3.9.jar" copyCaller \
            "$cnv_file" \
            --output-file        "$COPYNUMBERDIR/${base}.copynumber.called" \
            --output-homdel-file "$COPYNUMBERDIR/${base}.copynumber.homdel" \
            --amp-threshold      "$CNV_AMP_THRESHOLD" \
            --del-threshold      "$CNV_DEL_THRESHOLD" \
            --recenter-up        "$CNV_RECENTER_UP" \
            --recenter-down      "$CNV_RECENTER_DOWN") &
        job_pool_add "$!" || return 1
    done
    job_pool_wait_all || return 1
    log_message "copyCaller complete"
}

# =============================================================================
# STAGE 7: PREPARE VAR POSITION FILES
# Converts .hc.vcf files from Stage 4 into the 3-column position format
# (chrom, start, end) required by bam-readcount.
#
# SNP:   position is reported as-is
# INDEL: position is shifted by 1 for deletions so bam-readcount spans
#        the deleted base rather than the anchor base
#
# SNP files → snp-VAR/    (Somatic, Germline, and LOH)
# INDEL files → indel-VAR/ (Somatic only — Germline/LOH INDELs are not
#   processed through bam-readcount or fpfilter; their hc.vcf files are the
#   final output for those classes)
#
# Must run AFTER Stage 4 and after check_stage4_output() confirms .hc.vcf exist.
# =============================================================================

prepare_filter_input() {
    log_message "=== STAGE 7: Preparing VAR position files for bam-readcount ==="

    job_pool_init "$MAX_PARALLEL_JOBS"

    for hc_file in "$SOMATICDIR"/*.snp.*.hc.vcf; do
        [ -f "$hc_file" ] || continue
        local base out_var
        base=$(basename "$hc_file" .vcf)
        out_var="$SNPVARDIR/${base}.var"
        log_message "SNP VAR: $base"
        (awk 'BEGIN {OFS="\t"} !/^#/ { print $1, $2, $2 }' \
            "$hc_file" > "$out_var") &
        job_pool_add "$!" || return 1
    done

    # Somatic INDELs only — fpfilter is most critical for somatic INDELs where
    # alignment artefacts around indel breakpoints are the dominant source of
    # false positives. Germline and LOH INDELs are not fpfiltered.
    for hc_file in "$SOMATICDIR"/*.indel.Somatic.hc.vcf; do
        [ -f "$hc_file" ] || continue
        local base out_var
        base=$(basename "$hc_file" .vcf)
        out_var="$INDELVARDIR/${base}.var"
        log_message "INDEL VAR: $base"
        # Shift position by 1 for deletions so bam-readcount spans the deleted base
        (awk 'BEGIN {OFS="\t"} !/^#/ {
            is_del = (length($4) > length($5)) ? 1 : 0
            print $1, ($2 + is_del), ($2 + is_del)
        }' "$hc_file" > "$out_var") &
        job_pool_add "$!" || return 1
    done
    job_pool_wait_all || return 1
    log_message "VAR file preparation complete"
}

# =============================================================================
# STAGE 8: BAM-READCOUNT
# Generates per-position, strand-aware allele read counts used by fpfilter
# to remove false positives caused by mapping artefacts and strand bias.
#
# BAM selection by variant class:
#   Somatic SNP / INDEL  →  tumour BAM   (mutations arose in the tumour)
#   Germline SNP         →  normal BAM   (variants present in constitutional DNA)
#   LOH SNP              →  tumour BAM   (allele loss is a tumour-specific event)
#
# -q BRC_MAP_QUAL (10): excludes ambiguously mapped reads per VarScan FAQ.
# -b BRC_BASE_QUAL (15): fpfilter documentation recommended floor.
#
# Must run AFTER Stage 7.
# =============================================================================

run_bam_readcount() {
    log_message "=== STAGE 8: bam-readcount ==="

    check_file_exists "$FILE_PAIRS_LIST"

    job_pool_init "$MAX_PARALLEL_JOBS"
    while IFS=',' read -r entry1 entry2; do
        entry1=$(echo "$entry1" | tr -d ' \r')
        entry2=$(echo "$entry2" | tr -d ' \r')
        validate_pair "$entry1" "$entry2" || continue

        local samplen samplet normal_bam tumor_bam prefix
        samplet=$(get_sample_name "$entry1")    # column 1 = tumour
        samplen=$(get_sample_name "$entry2")    # column 2 = normal
        tumor_bam=$(get_bam_path "$entry1")
        normal_bam=$(get_bam_path "$entry2")
        prefix="${samplen}_${samplet}"

        # bam-readcount requires a .bai index; verify before launching jobs
        check_bam_index "$tumor_bam"
        check_bam_index "$normal_bam"

        # Somatic SNP — tumour BAM
        local var_file="$SNPVARDIR/${prefix}.snp.Somatic.hc.var"
        if [ -f "$var_file" ]; then
            log_message "bam-readcount: Somatic SNP — tumour ${samplet}"
            (bam-readcount \
                -q "$BRC_MAP_QUAL" -b "$BRC_BASE_QUAL" \
                -f "$GENOMEIDX1" -l "$var_file" \
                "$tumor_bam" > "$READCOUNTDIR/${prefix}.snp.Somatic.hc.readcount") &
            job_pool_add "$!" || return 1
        fi

        # Somatic INDEL — tumour BAM
        local indel_var="$INDELVARDIR/${prefix}.indel.Somatic.hc.var"
        if [ -f "$indel_var" ]; then
            log_message "bam-readcount: Somatic INDEL — tumour ${samplet}"
            (bam-readcount \
                -q "$BRC_MAP_QUAL" -b "$BRC_BASE_QUAL" \
                -f "$GENOMEIDX1" -l "$indel_var" \
                "$tumor_bam" > "$READCOUNTDIR/${prefix}.indel.Somatic.hc.readcount") &
            job_pool_add "$!" || return 1
        fi

        # Germline SNP — normal BAM
        local germ_var="$SNPVARDIR/${prefix}.snp.Germline.hc.var"
        if [ -f "$germ_var" ]; then
            log_message "bam-readcount: Germline SNP — normal ${samplen}"
            (bam-readcount \
                -q "$BRC_MAP_QUAL" -b "$BRC_BASE_QUAL" \
                -f "$GENOMEIDX1" -l "$germ_var" \
                "$normal_bam" > "$READCOUNTDIR/${prefix}.snp.Germline.hc.readcount") &
            job_pool_add "$!" || return 1
        fi

        # LOH SNP — tumour BAM (LOH is a tumour-specific allele loss event)
        local loh_var="$SNPVARDIR/${prefix}.snp.LOH.hc.var"
        if [ -f "$loh_var" ]; then
            log_message "bam-readcount: LOH SNP — tumour ${samplet}"
            (bam-readcount \
                -q "$BRC_MAP_QUAL" -b "$BRC_BASE_QUAL" \
                -f "$GENOMEIDX1" -l "$loh_var" \
                "$tumor_bam" > "$READCOUNTDIR/${prefix}.snp.LOH.hc.readcount") &
            job_pool_add "$!" || return 1
        fi

    done < "$FILE_PAIRS_LIST"
    job_pool_wait_all || return 1
    log_message "bam-readcount complete"
}

# =============================================================================
# STAGE 9: FALSE POSITIVE FILTER
# Applies fpfilter.pl to each readcount + VCF pair to remove variants with
# strand bias, low base quality, excessive mismatches near the variant,
# or short distance to the 3' read end.
#
# Must run AFTER Stage 8.
# =============================================================================

run_fpfilter() {
    log_message "=== STAGE 9: False positive filtering (fpfilter.pl) ==="

    check_file_exists "$SCRIPTSDIR/fpfilter.pl"

    job_pool_init "$MAX_PARALLEL_JOBS"
    for rc_file in "$READCOUNTDIR"/*.readcount; do
        [ -f "$rc_file" ] || continue
        local base vcf_file filtered_file
        base=$(basename "$rc_file" .readcount)
        vcf_file="$SOMATICDIR/${base}.vcf"
        if [ ! -f "$vcf_file" ]; then
            log_message "WARNING: No matching VCF for ${base} — skipping fpfilter"
            continue
        fi
        filtered_file="$FILTEREDDIR/${base}.fpfilter.vcf"
        log_message "fpfilter: $base"
        (perl "$SCRIPTSDIR/fpfilter.pl" \
            --vcf-var-file    "$vcf_file" \
            --readcount-file  "$rc_file" \
            --output-file     "$filtered_file" \
            --min-var-count    3 \
            --min-var-freq     "$MIN_VAR_FREQ" \
            --min-ref-basequal "$BRC_BASE_QUAL" \
            --min-var-basequal "$BRC_BASE_QUAL") &
        job_pool_add "$!" || return 1
    done
    job_pool_wait_all || return 1
    log_message "False positive filtering complete"
}

# =============================================================================
# STAGE 10: SUMMARY REPORT
# =============================================================================

generate_summary() {
    log_message "=== STAGE 10: Generating summary report ==="

    local summary_file="$VARSCANDIR/varscan_pipeline_summary.txt"
    {
        echo "========================================================================"
        echo "VarScan2 Pipeline Summary"
        echo "Generated: $(date)"
        echo "========================================================================"
        echo ""
        echo "CONFIGURATION"
        echo "  Reference genome       : $GENOMEIDX1"
        echo "  BAM directory          : $BAM_DIR"
        echo "  BAM suffix             : $BAM_SUFFIX"
        echo "  Pairs file             : $FILE_PAIRS_LIST"
        echo "  Pairs suffix (stripped): $PAIRS_SUFFIX"
        echo "  Target regions (WES)   : ${TARGET_BED:-not set (full-genome pileup)}"
        echo "  VCF sample list        : ${VCF_SAMPLE_LIST:-not set (generic NORMAL/TUMOR headers)}"
        echo ""
        echo "CALLING PARAMETERS"
        echo "  Min coverage           : $MIN_COVERAGE (normal: $MIN_COVERAGE_NORMAL, tumour: $MIN_COVERAGE_TUMOR)"
        echo "  Min base quality       : $MIN_BASE_QUAL"
        echo "  Min map quality        : $BRC_MAP_QUAL"
        echo "  Min variant freq       : $MIN_VAR_FREQ"
        echo "  P-value (calling)      : $P_VALUE"
        echo "  Somatic p-value        : $SOMATIC_P_VALUE"
        echo "  processSomatic p-value : $PROCESS_P_VALUE"
        echo "  CNV p-value            : $CNV_P_VALUE"
        echo ""
        echo "CNV CALLER THRESHOLDS"
        echo "  Amplification threshold: $CNV_AMP_THRESHOLD log2"
        echo "  Deletion threshold     : $CNV_DEL_THRESHOLD log2"
        echo "  Recenter up            : $CNV_RECENTER_UP"
        echo "  Recenter down          : $CNV_RECENTER_DOWN"
        echo ""
        echo "RESULTS"
        if [ -d "$SOMATICDIR" ]; then
            echo "  Somatic SNP  .hc  : $(find "$SOMATICDIR"   -name '*.snp.Somatic.hc.vcf'  -type f | wc -l) files"
            echo "  Somatic INDEL .hc : $(find "$SOMATICDIR"   -name '*.indel.Somatic.hc.vcf' -type f | wc -l) files"
            echo "  Germline .hc      : $(find "$SOMATICDIR"   -name '*.Germline.hc.vcf'      -type f | wc -l) files"
            echo "  LOH .hc           : $(find "$SOMATICDIR"   -name '*.LOH.hc.vcf'           -type f | wc -l) files"
        fi
        if [ -d "$COPYNUMBERDIR" ]; then
            echo "  Raw CNV           : $(find "$COPYNUMBERDIR" -name '*.copynumber'           -type f | wc -l) files"
            echo "  Called CNV        : $(find "$COPYNUMBERDIR" -name '*.copynumber.called'    -type f | wc -l) files"
            echo "  Homdel            : $(find "$COPYNUMBERDIR" -name '*.copynumber.homdel'    -type f | wc -l) files"
        fi
        if [ -d "$FILTEREDDIR" ]; then
            echo "  FP-filtered VCFs  : $(find "$FILTEREDDIR"   -name '*.fpfilter.vcf'         -type f | wc -l) files"
        fi
        echo ""
        echo "NEXT STEPS"
        echo "  1. Inspect *.copynumber files; set CNV_RECENTER_UP/DOWN if log2 baseline is offset"
        echo "  2. Run circular binary segmentation (CBS) on *.copynumber.called"
        echo "  3. Annotate variants with VEP or ANNOVAR"
        echo "  4. Validate high-impact somatic mutations"
        echo "========================================================================"
    } > "$summary_file"

    log_message "Summary written to: $summary_file"
    cat "$summary_file"
}

# =============================================================================
# HELP
# =============================================================================

show_help() {
    cat <<'EOF'
Usage: varscan_pipeline.sh [--help]

VarScan2 paired tumour/normal somatic variant and CNV analysis pipeline.
Runs 10 sequential stages from BAM flagstats through false-positive filtering.

PIPELINE STAGES
  Stage 1   samtools flagstat        Alignment QC; read counts for CNV data-ratio
  Stage 2   samtools mpileup         Paired pileup (normal first, tumour second)
  Stage 3   VarScan somatic          Raw SNP and INDEL calls (.snp.vcf, .indel.vcf)
  Stage 4   VarScan processSomatic   Classify: Somatic / Germline / LOH (.hc.vcf)
  Stage 5   VarScan copynumber       Per-window log2 depth ratios (.copynumber)
  Stage 6   VarScan copyCaller       Segmented CNV calls (.called, .homdel)
  Stage 7   awk                      Convert .hc.vcf to position lists (.var)
  Stage 8   bam-readcount            Per-position allele read support
  Stage 9   fpfilter.pl              False positive removal (.fpfilter.vcf)
  Stage 10  Summary report           varscan_pipeline_summary.txt

QUICK START
  1. Edit the CONFIGURATION section at the top of this script.
  2. Ensure your pairs file is correct (tumour col 1, normal col 2):
       50T_CRC_final.bam,50N_CRC_final.bam
  3. Run:
       chmod +x varscan_pipeline.sh
       ./varscan_pipeline.sh

REQUIRED CONFIGURATION (set before running)
  GENOMEIDX1        Indexed reference FASTA (must match BAM chromosome naming)
  SOFTWAREDIR       Directory containing VarScan.v2.3.9.jar
  SCRIPTSDIR        Directory containing fpfilter.pl
  BAM_DIR           Directory containing all BAM files
  BAM_SUFFIX        BAM filename suffix          e.g. _final.bam
  FILE_PAIRS_LIST   Tumour/normal pairs CSV       e.g. sample_pairs.csv
  PAIRS_SUFFIX      Suffix to strip from pairs entries to get sample name

OPTIONAL CONFIGURATION
  TARGET_BED        Capture kit BED file (WES only)
                    When set: restricts mpileup to captured regions, reduces
                    file size ~95%, excludes alt/decoy contigs automatically.
                    Leave empty for WGS.
  VCF_SAMPLE_LIST   Plain-text file with sample names (normal first, tumour
                    second, one per line) for correct VCF column labels.
                    Leave empty to use VarScan generic NORMAL/TUMOR headers.

KEY PARAMETERS (defaults shown)
  MIN_COVERAGE=20           Minimum site depth (both samples)
  MIN_COVERAGE_NORMAL=10    Minimum depth in normal
  MIN_COVERAGE_TUMOR=20     Minimum depth in tumour
  MIN_BASE_QUAL=20          Minimum base quality (mpileup + VarScan)
  MIN_VAR_FREQ=0.10         Minimum variant allele frequency to call
  TUMOR_PURITY=1.0          Tumour purity estimate — SET TO ACTUAL VALUE
                            VarScan uses this in its somatic Fisher's test.
                            Leaving at 1.0 for an impure tumour will
                            under-call subclonal variants.
  SOMATIC_P_VALUE=0.05      Significance threshold for somatic classification
  CNV_AMP_THRESHOLD=0.25    log2 ratio threshold for amplification call
  CNV_DEL_THRESHOLD=0.25    log2 ratio threshold for deletion call
  CNV_RECENTER_UP/DOWN=0    Baseline correction for chromosomally unstable
                            tumours — inspect *.copynumber before setting

PREREQUISITES
  Tool              Min version   Notes
  bash              4.0           Requires (( )) and [[ ]]
  Java JRE          8             Required for VarScan.v2.3.9.jar
  samtools          1.13          Must be on PATH
  bam-readcount     0.8           Must be on PATH
  Perl              5.10          Required for fpfilter.pl
  bc                any           Floating-point data-ratio calculation

  All BAM files must be coordinate-sorted, indexed (.bai present), and
  duplicate-marked (MarkDuplicates). Duplicates are excluded from pileup
  via -F 0x400. Chromosome naming in BAMs must match the reference FASTA
  (Ensembl: 1,2,...,X,Y,MT  vs  UCSC: chr1,chr2,...,chrX,chrY,chrM).

OUTPUT DIRECTORIES
  flagstats/        samtools flagstat outputs
  mpileup/          Paired .mpileup files
  somatic/          VarScan somatic and processSomatic outputs
  copynumber/       VarScan copynumber and copyCaller outputs
  snp-VAR/          SNP position lists (.var) for bam-readcount
  indel-VAR/        Somatic INDEL position lists (.var) for bam-readcount
  readcount/        bam-readcount outputs
  filtered/         fpfilter output VCFs (.fpfilter.vcf)
  .resume_state/    SHA256 stage markers for --resume validation

POST-PIPELINE STEPS
  1. Inspect *.copynumber files; adjust CNV_RECENTER_UP/DOWN if log2
     baseline is offset from zero (common in CIN tumours).
  2. Run circular binary segmentation (CBS) on *.copynumber.called.
  3. Annotate filtered VCFs with VEP or ANNOVAR.
  4. Validate high-impact somatic mutations.

BATCH PROCESSING
  Run all stages (default):
    ./varscan_pipeline.sh

  Run a single stage:
    ./varscan_pipeline.sh --stage 5

  Run a range of stages:
    ./varscan_pipeline.sh --from 3 --to 6

  Resume after a failed or partial run (skip stages with matching SHA256 state):
    ./varscan_pipeline.sh --resume

  Resume from a specific stage:
    ./varscan_pipeline.sh --from 5 --resume

  Preview what would run without executing anything:
    ./varscan_pipeline.sh --dry-run
    ./varscan_pipeline.sh --from 3 --to 8 --dry-run

  OPTIONS
    --stage N       Run only stage N (equivalent to --from N --to N)
    --from N        Start execution at stage N (default: 1)
    --to N          Stop execution after stage N (default: 10)
    --resume        Skip a stage only when both are true:
                    1) output files are present, and
                    2) a stored SHA256 marker matches current inputs/params.
    --dry-run       Print which stages would execute without running them.
                    Prerequisite file checks are also skipped in dry-run mode.
    -h, --help      Show this help message and exit

  STAGE DEPENDENCIES (do not skip intermediate stages in a partial run)
    Stage 2 → requires Stage 1 flagstats (data-ratio for CNV)
    Stage 3 → requires Stage 2 mpileup
    Stage 4 → requires Stage 3 VCF output
    Stages 7–9 → require Stage 4 .hc.vcf output
    Stage 8 → requires Stage 7 .var files
    Stage 9 → requires Stage 8 .readcount files

SUPPORT
  See README.md for full parameter justifications and troubleshooting.
EOF
}

# =============================================================================
# MAIN
# =============================================================================

main() {
    # --------------------------------------------------------------------------
    # Argument parsing
    # --------------------------------------------------------------------------
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                show_help; exit 0 ;;
            --stage)
                if [ -z "${2:-}" ]; then
                    log_message "ERROR: --stage requires a value (1-10)"
                    exit 1
                fi
                FROM_STAGE="$2"; TO_STAGE="$2"; shift 2 ;;
            --from)
                if [ -z "${2:-}" ]; then
                    log_message "ERROR: --from requires a value (1-10)"
                    exit 1
                fi
                FROM_STAGE="$2"; shift 2 ;;
            --to)
                if [ -z "${2:-}" ]; then
                    log_message "ERROR: --to requires a value (1-10)"
                    exit 1
                fi
                TO_STAGE="$2"; shift 2 ;;
            --resume)
                RESUME=true; shift ;;
            --dry-run)
                DRY_RUN=true; shift ;;
            *)
                log_message "ERROR: Unknown argument: $1"
                echo "Run './varscan_pipeline.sh --help' for usage."
                exit 1 ;;
        esac
    done

    # --------------------------------------------------------------------------
    # Validate stage range
    # --------------------------------------------------------------------------
    if ! [[ "$FROM_STAGE" =~ ^([1-9]|10)$ ]] || ! [[ "$TO_STAGE" =~ ^([1-9]|10)$ ]]; then
        log_message "ERROR: Stage numbers must be integers between 1 and 10"
        exit 1
    fi
    if [ "$FROM_STAGE" -gt "$TO_STAGE" ]; then
        log_message "ERROR: --from ($FROM_STAGE) must not exceed --to ($TO_STAGE)"
        exit 1
    fi

    log_message "Starting VarScan2 Pipeline"
    log_message "Stages: ${FROM_STAGE}–${TO_STAGE}  |  Resume: ${RESUME}  |  Dry-run: ${DRY_RUN}"

    # --------------------------------------------------------------------------
    # Prerequisite file checks (skipped in dry-run — files may not exist yet)
    # --------------------------------------------------------------------------
    if [ "$DRY_RUN" = false ]; then
        check_command_exists samtools
        check_command_exists java
        check_command_exists bam-readcount
        check_command_exists perl
        check_command_exists bc
        check_command_exists sha256sum
        check_file_exists "$FILE_PAIRS_LIST"
        check_file_exists "$GENOMEIDX1"
        check_file_exists "$SOFTWAREDIR/VarScan.v2.3.9.jar"
        check_file_exists "$SCRIPTSDIR/fpfilter.pl"
    fi

    setup_directories

    # --------------------------------------------------------------------------
    # Stage execution
    # run_stage <number> <function> "<description>"
    # --------------------------------------------------------------------------
    run_stage  1  generate_flagstats        "samtools flagstat   — read counts for CNV data-ratio" || exit 1
    run_stage  2  generate_mpileup          "samtools mpileup    — paired pileup (normal first)" || exit 1
    run_stage  3  run_varscan_somatic       "VarScan somatic     — raw SNP and INDEL calls" || exit 1
    run_stage  4  process_somatic_variants  "processSomatic      — Somatic / Germline / LOH" || exit 1

    # Stage 4 guard: verify .hc.vcf output exists before stages 7–9 consume it.
    # Runs whenever stage 7 or later is in scope and execution is not a dry-run.
    if [ "$TO_STAGE" -ge 7 ] && [ "$DRY_RUN" = false ]; then
        check_stage4_output
    fi

    run_stage  5  run_varscan_copynumber    "VarScan copynumber  — per-window log2 depth ratios" || exit 1
    run_stage  6  run_copy_caller           "VarScan copyCaller  — segmented CNV calls" || exit 1
    run_stage  7  prepare_filter_input      "awk                 — convert .hc.vcf to .var" || exit 1
    run_stage  8  run_bam_readcount         "bam-readcount       — per-position allele counts" || exit 1
    run_stage  9  run_fpfilter              "fpfilter.pl         — false positive removal" || exit 1
    run_stage 10  generate_summary          "summary report" || exit 1

    if [ "$DRY_RUN" = false ]; then
        log_message "VarScan2 Pipeline complete"
    fi
}

main "$@"
