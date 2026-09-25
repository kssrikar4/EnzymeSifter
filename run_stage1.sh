#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: ./run_stage1.sh <input> [options]
       ./run_stage1.sh -l <list_file> [options]

Stage 1: filter and cluster input amino-acid sequences, producing a
non-redundant FASTA for structure prediction.

<input> can be either:
  - a single FASTA file (.fasta/.fa/.faa), or
  - a directory containing one or more FASTA files (concatenated).

Multi-FASTA mode:
  -l, --list <file>   Text file containing file paths of FASTA files (one per line).
                      Creates a folder named after each FASTA file and runs the
                      analysis inside that folder.

Options:
  -identity <pct>     Cluster at <pct>% identity (0-100), keeping one
                      representative per cluster.
  -residues <motif>   Keep only sequences matching a motif.
                      Regex syntax: '.' = any single residue.
                      Examples: G.S.G, SHD
  -pfam <IDs>         Keep only sequences that hit any of the given
                      Pfam accessions (comma-separated, OR logic).
                      Uses Pfam gathering thresholds via hmmsearch.
                      Examples: PF07519
                                PF07519,PF00657
                      First run auto-downloads Pfam (~1.5 GB).

  -ec <IDs>           Keep only sequences whose CLEAN-predicted EC
                      number matches any of the given EC specs
                      (comma-separated, OR logic). Partial specs
                      allowed: a missing or '-' digit is a wildcard.
                      Examples: 3.13.1.8
                                3.13.-.-           (same as 3.13)
                                3.13.1.8,2.5.1.94
                      First run clones CLEAN, builds the package, and
                      downloads pretrained weights from Google Drive
                      (~few hundred MB) plus ESM-1b on first inference
                      (~7 GB to ~/.cache/torch).

  -threads <n>        Number of CPU cores to use.
                      Default: all available cores ($(nproc)).

When combined, filters apply in this order:
  residues -> pfam -> ec -> identity clustering

Examples:
  ./run_stage1.sh seqs.fasta -pfam PF07519
  ./run_stage1.sh -l data.lst -identity 90
  ./run_stage1.sh -l data.lst -residues G.S.G -ec 3.13.1.8
  ./run_stage1.sh seqs.fasta -residues G.S.G -ec 3.13.1.8,2.5.1.94
  ./run_stage1.sh seqs.fasta -pfam PF07519 -ec 2.5.1.94 -identity 90
  ./run_stage1.sh seqs.fasta -pfam PF07519 -threads 8

Outputs:
  data/stage1/nonredundant.fasta       Non-redundant sequences
  data/stage1/clustering_report.tsv    Representative-to-member mapping
  data/stage1/motif_report.tsv         motif match results
  data/stage1/pfam_report.tsv          Pfam match results
  data/stage1/ec_report.tsv            EC prediction results
  (In -l mode, outputs are placed under <sample_name>/data/stage1/ and symlinked to <sample_name>/)
EOF
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
LIST_FILE=""
IDENTITY=""
RESIDUES=""
PFAM=""
EC=""
THREADS="$(nproc)"

while [ $# -gt 0 ]; do
    case "$1" in
        -l|--list)
            [ -z "${2:-}" ] && echo "Error: -l requires a file path" >&2 && exit 1
            LIST_FILE="$2"
            shift 2 ;;
        -identity)
            [ -z "${2:-}" ] && echo "Error: -identity requires a value" >&2 && exit 1
            IDENTITY="$2"
            shift 2 ;;
        -residues)
            [ -z "${2:-}" ] && echo "Error: -residues requires a value" >&2 && exit 1
            RESIDUES="$2"
            shift 2 ;;
        -pfam)
            [ -z "${2:-}" ] && echo "Error: -pfam requires a value" >&2 && exit 1
            PFAM="$2"
            shift 2 ;;
        -ec)
            [ -z "${2:-}" ] && echo "Error: -ec requires a value" >&2 && exit 1
            EC="$2"
            shift 2 ;;
        -threads)
            [ -z "${2:-}" ] && echo "Error: -threads requires a value" >&2 && exit 1
            THREADS="$2"
            shift 2 ;;
        -h|--help)
            usage ;;
        -*)
            echo "Unknown option: $1" >&2
            usage ;;
        *)
            if [ -z "$INPUT_FASTA" ] && [ -z "$LIST_FILE" ]; then
                INPUT_FASTA="$1"
                shift
            else
                echo "Error: Unexpected argument: $1" >&2
                usage
            fi ;;
    esac
done

if [ -z "$INPUT_FASTA" ] && [ -z "$LIST_FILE" ]; then
    usage
fi

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

    # Ensure Snakemake is available. Auto-installs it into the 'enzymesifter'
    # conda environment on first run; sets $ENV_NAME for the call below.
    source "${SCRIPT_DIR}/scripts/ensure_snakemake.sh"

    for i in "${!files[@]}"; do
        raw_path="${files[$i]}"
        abs_fasta="$(resolve_path "$raw_path" "$LIST_DIR")"
        if [ ! -e "$abs_fasta" ]; then
            echo "Error: FASTA file does not exist: $raw_path (resolved to $abs_fasta)" >&2
            exit 1
        fi

        sample_name="$(get_sample_name "$abs_fasta")"
        echo ""
        echo "================================================================================"
        echo "[INFO] [$((i+1))/${#files[@]}] Processing sample '${sample_name}': ${abs_fasta}"
        echo "================================================================================"

        sample_out="${sample_name}/data/stage1"
        sample_log="${sample_name}/logs/stage1"
        mkdir -p "${sample_out}" "${sample_log}"

        SAMPLE_CONFIG="input_fasta=${abs_fasta} out_dir=${sample_out} log_dir=${sample_log}"
        [ -n "$IDENTITY" ] && SAMPLE_CONFIG+=" identity=${IDENTITY}"
        [ -n "$RESIDUES" ] && SAMPLE_CONFIG+=" residues=${RESIDUES}"
        [ -n "$PFAM" ]     && SAMPLE_CONFIG+=" pfam=${PFAM}"
        [ -n "$EC" ]       && SAMPLE_CONFIG+=" ec=${EC}"

        conda run -n "${ENV_NAME}" --no-capture-output \
            snakemake --snakefile "${SCRIPT_DIR}/Snakefile_stage1" --config ${SAMPLE_CONFIG} \
                -j "${THREADS}" \
                --quiet rules progress

        # Create convenience symlinks directly under ${sample_name}/
        (
            cd "${sample_name}"
            ln -sf "data/stage1/nonredundant.fasta" "nonredundant.fasta" 2>/dev/null || true
            ln -sf "data/stage1/clustering_report.tsv" "clustering_report.tsv" 2>/dev/null || true
            ln -sf "data/stage1/motif_report.tsv" "motif_report.tsv" 2>/dev/null || true
            ln -sf "data/stage1/pfam_report.tsv" "pfam_report.tsv" 2>/dev/null || true
            ln -sf "data/stage1/ec_report.tsv" "ec_report.tsv" 2>/dev/null || true
        )
    done

    echo ""
    echo "[INFO] Successfully completed Stage 1 for all ${#files[@]} sample(s)."
    exit 0
fi

CONFIG="input_fasta=${INPUT_FASTA}"
[ -n "$IDENTITY" ] && CONFIG+=" identity=${IDENTITY}"
[ -n "$RESIDUES" ] && CONFIG+=" residues=${RESIDUES}"
[ -n "$PFAM" ]     && CONFIG+=" pfam=${PFAM}"
[ -n "$EC" ]       && CONFIG+=" ec=${EC}"

# Ensure Snakemake is available. Auto-installs it into the 'enzymesifter'
# conda environment on first run; sets $ENV_NAME for the call below.
source "${SCRIPT_DIR}/scripts/ensure_snakemake.sh"

conda run -n "${ENV_NAME}" --no-capture-output \
    snakemake --snakefile "${SCRIPT_DIR}/Snakefile_stage1" --config ${CONFIG} \
        -j "${THREADS}" \
        --quiet rules progress
