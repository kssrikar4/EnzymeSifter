#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<EOF
Usage: ./run_stage2.sh /path/to/pdbs [filter options] [clade options]
       ./run_stage2.sh -l <list_file> [filter options] [clade options]

Stage 2: structural screening, property prediction (solubility, pH/T optimum, Tm),
phylogenetic tree generation, and clade representative selection.

Multi-sample mode:
  -l, --list <file>        Text file containing file paths of FASTA files (or PDB directories).
                           Runs Stage 2 analysis for each sample in its own directory
                           named after the sample (<sample_name>/).
                           Automatically detects PDBs under <sample_name>/data/predicted_pdbs/

Filter options (all optional):
  -solubility <min>        Minimum predicted solubility   (e.g. 0.4)
  -usability  <min>        Minimum predicted usability    (e.g. 0.35)
  -tm  <min>  | <lo:hi>    Melting temp: cutoff OR interval
                             e.g. -tm 55         → Tm >= 55, higher better
                             e.g. -tm 50:80      → 50 <= Tm <= 80, mid best
  -topt  <lo:hi>           Optimal temp interval (required form)
                             e.g. -topt 35:42
  -phopt <lo:hi>           Optimal pH interval (required form)
                             e.g. -phopt 7.0:8.5

Clade options (optional):
  -clades <n>              Partition the NJ tree into <n> clades.
                           If combined with filter options, writes
                           predictions_output/clade_representatives.tsv
                           with the best-scoring enzyme per clade.

Other options:
  -threads <n>             Number of CPU cores to use.
                           Default: all available cores ($(nproc)).

Scoring (when filters are given):
  Each passing enzyme receives a score in [0, 1]:
    cutoff properties:   normalised across the passing pool (higher=better)
    interval properties: closeness to the interval midpoint
  Combined score = mean of per-property scores over user-specified
  properties with non-NA values for that enzyme.

Examples:
  ./run_stage2.sh /data/pdbs
  ./run_stage2.sh -l data.lst -solubility 0.69 -clades 5
  ./run_stage2.sh -l data.lst -usability 0.35 -tm 50:80 -topt 35:42 -clades 10
  ./run_stage2.sh /data/pdbs -usability 0.35 -solubility 0.4
  ./run_stage2.sh /data/pdbs -tm 55
  ./run_stage2.sh /data/pdbs -tm 50:80 -topt 35:42 -phopt 7:8
  ./run_stage2.sh /data/pdbs -clades 7
  ./run_stage2.sh /data/pdbs -clades 7 -threads 8
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

find_pdb_dir() {
    local target="$1"
    local sample_name="$2"

    if [ -d "$target" ] && compgen -G "$target/*.pdb" >/dev/null; then
        echo "$target"
        return 0
    fi
    if [ -d "${sample_name}/data/predicted_pdbs" ] && compgen -G "${sample_name}/data/predicted_pdbs/*.pdb" >/dev/null; then
        echo "${sample_name}/data/predicted_pdbs"
        return 0
    fi
    if [ -d "${sample_name}/predicted_pdbs" ] && compgen -G "${sample_name}/predicted_pdbs/*.pdb" >/dev/null; then
        echo "${sample_name}/predicted_pdbs"
        return 0
    fi
    if [ -d "${sample_name}/pdbs" ] && compgen -G "${sample_name}/pdbs/*.pdb" >/dev/null; then
        echo "${sample_name}/pdbs"
        return 0
    fi
    if [ -d "${sample_name}" ] && compgen -G "${sample_name}/*.pdb" >/dev/null; then
        echo "${sample_name}"
        return 0
    fi
    return 1
}

PDB_DIR=""
LIST_FILE=""
FILTER_SOLUBILITY=""
FILTER_USABILITY=""
FILTER_PHOPT=""
FILTER_TOPT=""
FILTER_TM=""
CLADES=""
THREADS="$(nproc)"

while [ $# -gt 0 ]; do
    case "$1" in
        -l|--list)
            [ -z "${2:-}" ] && echo "Error: -l requires a file path" >&2 && exit 1
            LIST_FILE="$2"
            shift 2 ;;
        -solubility)
            [ -z "${2:-}" ] && echo "Error: -solubility requires a value" >&2 && exit 1
            FILTER_SOLUBILITY="$2"
            shift 2 ;;
        -usability)
            [ -z "${2:-}" ] && echo "Error: -usability requires a value" >&2 && exit 1
            FILTER_USABILITY="$2"
            shift 2 ;;
        -phopt)
            [ -z "${2:-}" ] && echo "Error: -phopt requires a value" >&2 && exit 1
            FILTER_PHOPT="$2"
            shift 2 ;;
        -topt)
            [ -z "${2:-}" ] && echo "Error: -topt requires a value" >&2 && exit 1
            FILTER_TOPT="$2"
            shift 2 ;;
        -tm)
            [ -z "${2:-}" ] && echo "Error: -tm requires a value" >&2 && exit 1
            FILTER_TM="$2"
            shift 2 ;;
        -clades)
            [ -z "${2:-}" ] && echo "Error: -clades requires a value" >&2 && exit 1
            CLADES="$2"
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
            if [ -z "$PDB_DIR" ] && [ -z "$LIST_FILE" ]; then
                PDB_DIR="$1"
                shift
            else
                echo "Error: Unexpected argument: $1" >&2
                usage
            fi ;;
    esac
done

if [ -z "$PDB_DIR" ] && [ -z "$LIST_FILE" ]; then
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
        echo "Error: No valid paths found in $LIST_FILE" >&2
        exit 1
    fi

    echo "[INFO] Found ${#files[@]} sample(s) in $LIST_FILE"

    for i in "${!files[@]}"; do
        raw_path="${files[$i]}"
        abs_entry="$(resolve_path "$raw_path" "$LIST_DIR")"
        sample_name="$(get_sample_name "$abs_entry")"

        pdb_path="$(find_pdb_dir "$abs_entry" "$sample_name" || true)"
        if [ -z "$pdb_path" ]; then
            echo "Error: No PDB directory containing .pdb files found for sample '${sample_name}'." >&2
            echo "Checked: '${abs_entry}', '${sample_name}/data/predicted_pdbs', '${sample_name}/predicted_pdbs', '${sample_name}/pdbs', and '${sample_name}'." >&2
            echo "Please ensure PDB prediction has been run." >&2
            exit 1
        fi

        echo ""
        echo "================================================================================"
        echo "[INFO] [$((i+1))/${#files[@]}] Running Stage 2 for '${sample_name}': PDBs from ${pdb_path} -> output in ${sample_name}/"
        echo "================================================================================"

        SAMPLE_CONFIG="pdb_dir=${pdb_path} out_dir=${sample_name}"
        [ -n "$FILTER_SOLUBILITY" ] && SAMPLE_CONFIG+=" filter_solubility=${FILTER_SOLUBILITY}"
        [ -n "$FILTER_USABILITY" ]  && SAMPLE_CONFIG+=" filter_usability=${FILTER_USABILITY}"
        [ -n "$FILTER_PHOPT" ]      && SAMPLE_CONFIG+=" filter_phopt=${FILTER_PHOPT}"
        [ -n "$FILTER_TOPT" ]       && SAMPLE_CONFIG+=" filter_topt=${FILTER_TOPT}"
        [ -n "$FILTER_TM" ]         && SAMPLE_CONFIG+=" filter_tm=${FILTER_TM}"
        [ -n "$CLADES" ]            && SAMPLE_CONFIG+=" clades=${CLADES}"

        conda run -n "${ENV_NAME}" --no-capture-output \
            snakemake --snakefile "${SCRIPT_DIR}/Snakefile_stage2" --config ${SAMPLE_CONFIG} \
                -j "${THREADS}" \
                --quiet rules progress
    done

    echo ""
    echo "[INFO] Successfully completed Stage 2 for all ${#files[@]} sample(s)."
    exit 0
fi

CONFIG="pdb_dir=${PDB_DIR}"
[ -n "$FILTER_SOLUBILITY" ] && CONFIG+=" filter_solubility=${FILTER_SOLUBILITY}"
[ -n "$FILTER_USABILITY" ]  && CONFIG+=" filter_usability=${FILTER_USABILITY}"
[ -n "$FILTER_PHOPT" ]      && CONFIG+=" filter_phopt=${FILTER_PHOPT}"
[ -n "$FILTER_TOPT" ]       && CONFIG+=" filter_topt=${FILTER_TOPT}"
[ -n "$FILTER_TM" ]         && CONFIG+=" filter_tm=${FILTER_TM}"
[ -n "$CLADES" ]            && CONFIG+=" clades=${CLADES}"

conda run -n "${ENV_NAME}" --no-capture-output \
    snakemake --snakefile "${SCRIPT_DIR}/Snakefile_stage2" --config ${CONFIG} \
        -j "${THREADS}" \
        --quiet rules progress
