#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<USAGE
Usage: ./run_pdb_prediction.sh <input_fasta> <output_dir>
       ./run_pdb_prediction.sh -l <list_file> [output_dir]

Predicts 3D PDB structures from an input FASTA using ESMFold.
Automatically uses GPU (and handles low VRAM via FP16/chunking) or falls back to CPU.

Arguments:
  <input_fasta>       Path to the input FASTA file (e.g., data/stage1/nonredundant.fasta)
  <output_dir>        Directory to save the resulting .pdb files

Multi-FASTA mode:
  -l, --list <file>   Text file containing file paths of FASTA files (one per line).
                      Predicts structures for each FASTA file, saving PDBs under:
                        <sample_name>/data/predicted_pdbs/ (default) or <output_dir>/<sample_name>/
                      If <sample_name>/data/stage1/nonredundant.fasta exists, it is
                      automatically used as the filtered input for prediction.
USAGE
    exit 1
}

get_sample_name() {
    local p="$1"
    local b
    b="$(basename "$p")"
    if [[ "$b" == "predicted_pdbs" || "$b" == "pdbs" ]]; then
        local parent
        parent="$(basename "$(dirname "$p")")"
        if [[ "$parent" == "data" ]]; then
            b="$(basename "$(dirname "$(dirname "$p")")")"
        elif [[ -n "$parent" && "$parent" != "." && "$parent" != "/" ]]; then
            b="$parent"
        fi
    fi
    if [[ "$b" == *.tar.gz || "$b" == *.fasta.gz || "$b" == *.fa.gz || "$b" == *.faa.gz ]]; then
        b="${b%.*}"
        b="${b%.*}"
    else
        b="${b%.*}"
    fi
    echo "$b"
}

resolve_path() {
    local p="$1"
    local list_dir="$2"
    if [ -e "$p" ]; then
        echo "$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
    elif [ -n "$list_dir" ] && [ -e "$list_dir/$p" ]; then
        echo "$(cd "$list_dir/$(dirname "$p")" && pwd)/$(basename "$p")"
    else
        echo "$p"
    fi
}

INPUT_FASTA=""
OUTPUT_DIR=""
LIST_FILE=""

while [ $# -gt 0 ]; do
    case "$1" in
        -l|--list)
            [ -z "${2:-}" ] && echo "Error: -l requires a file path" >&2 && exit 1
            LIST_FILE="$2"
            shift 2 ;;
        -h|--help)
            usage ;;
        -*)
            echo "Unknown option: $1" >&2
            usage ;;
        *)
            if [ -z "$INPUT_FASTA" ] && [ -z "$LIST_FILE" ]; then
                INPUT_FASTA="$1"
            elif [ -z "$OUTPUT_DIR" ]; then
                OUTPUT_DIR="$1"
            else
                echo "Error: Unexpected argument: $1" >&2
                usage
            fi
            shift ;;
    esac
done

if [ -z "$LIST_FILE" ] && ([ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_DIR" ]); then
    usage
fi

source "${SCRIPT_DIR}/scripts/ensure_snakemake.sh"

if [ -n "$LIST_FILE" ]; then
    if [ ! -f "$LIST_FILE" ]; then
        echo "Error: List file does not exist: $LIST_FILE" >&2
        exit 1
    fi

    LIST_DIR="$(cd "$(dirname "$LIST_FILE")" && pwd)"
    files=()
    while IFS= read -r line || [ -n "$line" ]; do
        line="$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
        [ -z "$line" ] && continue
        [[ "$line" =~ ^# ]] && continue
        files+=("$line")
    done < "$LIST_FILE"

    if [ ${#files[@]} -eq 0 ]; then
        echo "Error: No valid file paths found in $LIST_FILE" >&2
        exit 1
    fi

    echo "[INFO] Found ${#files[@]} FASTA file(s) in $LIST_FILE"

    for i in "${!files[@]}"; do
        raw_path="${files[$i]}"
        abs_fasta="$(resolve_path "$raw_path" "$LIST_DIR")"
        sample_name="$(get_sample_name "$abs_fasta")"

        # Locate appropriate input FASTA: prefer stage 1 nonredundant FASTA if available
        target_fasta=""
        if [ -f "${sample_name}/data/stage1/nonredundant.fasta" ]; then
            target_fasta="${sample_name}/data/stage1/nonredundant.fasta"
        elif [ -f "${sample_name}/nonredundant.fasta" ]; then
            target_fasta="${sample_name}/nonredundant.fasta"
        elif [ -f "$abs_fasta" ]; then
            target_fasta="$abs_fasta"
        else
            echo "Error: Could not find input FASTA for sample '${sample_name}'. Checked ${sample_name}/data/stage1/nonredundant.fasta and $abs_fasta" >&2
            exit 1
        fi

        if [ -n "$OUTPUT_DIR" ]; then
            sample_out="${OUTPUT_DIR}/${sample_name}"
        else
            sample_out="${sample_name}/data/predicted_pdbs"
        fi
        mkdir -p "${sample_out}"

        echo ""
        echo "================================================================================"
        echo "[INFO] [$((i+1))/${#files[@]}] Predicting structures for '${sample_name}': ${target_fasta} -> ${sample_out}"
        echo "================================================================================"

        conda run -n "${ENV_NAME}" --no-capture-output \
            python "${SCRIPT_DIR}/scripts/predict_structures.py" "${target_fasta}" "${sample_out}"

        # Create convenience symlinks if using standard folder structure
        if [ -z "$OUTPUT_DIR" ]; then
            (
                cd "${sample_name}"
                ln -sfn "data/predicted_pdbs" "predicted_pdbs" 2>/dev/null || true
                ln -sfn "data/predicted_pdbs" "pdbs" 2>/dev/null || true
            )
        fi
    done

    echo ""
    echo "[INFO] Successfully completed structure prediction for all ${#files[@]} sample(s)."
    exit 0
fi

conda run -n "${ENV_NAME}" --no-capture-output \
    python "${SCRIPT_DIR}/scripts/predict_structures.py" "${INPUT_FASTA}" "${OUTPUT_DIR}"
