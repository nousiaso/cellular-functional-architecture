#!/usr/bin/env bash
#SBATCH --job-name=cfa_python
#SBATCH --partition=normal_prio
#SBATCH --account=no_limit_account
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
ROOT="$(git -C "$SLURM_SUBMIT_DIR" rev-parse --show-toplevel)"
cd "$ROOT"
[[ $# -ge 1 ]] || { echo "Missing Python stage argument" >&2; exit 2; }
STAGE="$1"
shift
case "$STAGE" in scripts/diagnostics/*.py|scripts/acquisition/*.py|scripts/staging/*.py) ;; *) echo "Unsupported stage path: $STAGE" >&2; exit 2;; esac
STAGE="$(realpath "$STAGE")"
[[ "$STAGE" == "$ROOT"/scripts/* && -f "$STAGE" ]] || exit 2
git diff --quiet && git diff --cached --quiet || { echo "Tracked files changed; inspect git status." >&2; exit 2; }
PYTHON="$(command -v "${CFA_PYTHON:-python}")"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS" MKL_NUM_THREADS="$OMP_NUM_THREADS"
finish() {
    rc=$?
    printf 'CFA_JOB_END status=%s utc=%s\n' "$rc" "$(date -u +%FT%TZ)"
    exit "$rc"
}
trap finish EXIT
printf 'CFA_JOB_START job=%s node=%s utc=%s\n' "$SLURM_JOB_ID" "$(hostname)" "$(date -u +%FT%TZ)"
printf 'GIT_HEAD %s\n' "$(git rev-parse HEAD)"
printf 'PYTHON %s\n' "$PYTHON"
"$PYTHON" --version
printf 'STAGE_SHA256 '; sha256sum "$STAGE"
printf 'RESOURCES cpus=%s mem_per_node_mb=%s\n' "$SLURM_CPUS_PER_TASK" "${SLURM_MEM_PER_NODE:-unknown}"
printf 'COMMAND '; printf '%q ' "$PYTHON" -u "$STAGE" "$@"; printf '\n'
if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v "$PYTHON" -u "$STAGE" "$@"
else
    "$PYTHON" -u "$STAGE" "$@"
fi
