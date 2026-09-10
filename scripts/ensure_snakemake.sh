#!/usr/bin/env bash

ENV_NAME="enzymesifter"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v conda &>/dev/null; then
    echo "Error: conda was not found on your PATH." >&2
    exit 1
fi

if ! conda env list | grep -qE "^[[:space:]]*${ENV_NAME}[[:space:]]"; then
    echo "[setup] First run: creating the unified '${ENV_NAME}' environment..." >&2
    conda env create -f "${SCRIPT_DIR}/envs/enzymesifter.yaml"
else
    echo "[setup] Verifying unified '${ENV_NAME}' environment packages..." >&2
    conda env update -n "${ENV_NAME}" -f "${SCRIPT_DIR}/envs/enzymesifter.yaml"
fi
