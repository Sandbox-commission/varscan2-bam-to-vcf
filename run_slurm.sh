#!/usr/bin/env bash
#SBATCH --job-name=varscan2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=varscan2_%j.log

set -euo pipefail

# Paths — adjust for your HPC
SIF="${SIF:-varscan2.sif}"
CONFIG="${CONFIG:-config.toml}"
BAM_DIR="${BAM_DIR:-/data/bams}"
REF_DIR="${REF_DIR:-/data/ref}"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: Apptainer image not found: $SIF" >&2
    echo "Pull with: apptainer pull varscan2.sif docker://ghcr.io/sandbox-commission/varscan2-pipeline:latest" >&2
    exit 1
fi

apptainer run \
    --bind "${BAM_DIR}:/data/bams" \
    --bind "${REF_DIR}:/data/ref" \
    --bind "$(pwd):/workspace" \
    --pwd /workspace \
    "$SIF" \
    --config "$CONFIG" \
    --resume \
    "$@"
