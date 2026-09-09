#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<USAGE
Usage: ./run_pdb_prediction.sh <input_fasta> <output_dir>

Predicts 3D PDB structures from an input FASTA using ESMFold.
Automatically uses GPU (and handles low VRAM via FP16/chunking) or falls back to CPU.

Arguments:
  <input_fasta>  Path to the input FASTA file (e.g., data/stage1/nonredundant.fasta)
  <output_dir>   Directory to save the resulting .pdb files
USAGE
    exit 1
}

if [ -z "${1:-}" ] || [ -z "${2:-}" ]; then
    usage
fi

INPUT_FASTA="$1"
OUTPUT_DIR="$2"

source "${SCRIPT_DIR}/scripts/ensure_snakemake.sh"

conda run -n "${ENV_NAME}" --no-capture-output python "${SCRIPT_DIR}/scripts/predict_structures.py" "${INPUT_FASTA}" "${OUTPUT_DIR}"
