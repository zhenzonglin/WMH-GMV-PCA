#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${project_dir}/src${PYTHONPATH:+:${PYTHONPATH}}"
bootstrap_python="${project_dir}/.venv/bin/python"
if [[ -r "${project_dir}/.wmh-gmv-python" ]]; then
  IFS= read -r bootstrap_python < "${project_dir}/.wmh-gmv-python"
fi
if [[ ! -x "${bootstrap_python}" ]]; then
  printf '%s\n' 'Project bootstrap environment missing. Run bash scripts/setup_project.sh first.' >&2
  exit 2
fi
exec "${bootstrap_python}" -m wmh_gmv.launcher stop "$@"
