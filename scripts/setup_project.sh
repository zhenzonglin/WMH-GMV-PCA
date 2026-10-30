#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${ROOT}/config/project.yaml"
if [[ "$#" -gt 0 ]]; then
  [[ "$#" -eq 2 && "$1" == "--config" ]] || { echo 'Usage: setup_project.sh --config CONFIG'; exit 2; }
  CONFIG="$2"
fi
command -v uv >/dev/null || { echo 'Install uv or add it to PATH; see docs/REMOTE_RUN.md'; exit 2; }
uv run --no-project --with PyYAML "${ROOT}/scripts/bootstrap_project.py" --config "$CONFIG"
