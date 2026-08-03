#!/usr/bin/env bash
# Attach to a running (or completed) pipeline log.
# Usage:
#   ./watch_run.sh                        # reattach to the most recent run
#   ./watch_run.sh <VM_NAME>              # attach to a specific run
#   ./watch_run.sh <VM_NAME> <PROJECT>    # ... in a project other than the last run's
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# .last_run holds the VM name on the first line and the project it was created in on the
# second, both written by `./slidr --gcp`. Older single-line files still work: the project
# then falls back to $SLIDR_PROJECT or the active gcloud config.
LAST_RUN="${SCRIPT_DIR}/.last_run"
VM_NAME="${1:-$(sed -n '1p' "$LAST_RUN" 2>/dev/null || true)}"
PROJECT="${2:-${SLIDR_PROJECT:-}}"
[[ -z "${PROJECT}" ]] && PROJECT="$(sed -n '2p' "$LAST_RUN" 2>/dev/null || true)"
[[ -z "${PROJECT}" ]] && PROJECT="$(gcloud config get-value project 2>/dev/null || true)"
[[ "${PROJECT}" == "(unset)" ]] && PROJECT=""
PROJECT_HINT="${PROJECT:-<PROJECT>}"

[[ -z "${VM_NAME}" ]] && {
  echo "[ERROR]: no active run found." >&2
  echo "Troubleshooting:" >&2
  echo " • Start a run with: ./slidr --bcl <BCL_ID> --gcp --input gs://..." >&2
  echo " • ${SCRIPT_DIR}/.last_run is written by ./slidr --gcp; it is absent because no GCP run has been launched from this checkout" >&2
  echo " • To attach to a run launched elsewhere, pass the VM name: ./watch_run.sh <VM_NAME>" >&2
  echo " • List the running VMs: gcloud compute instances list --project=${PROJECT_HINT} --filter='name~slidr-run'" >&2
  exit 1
}

# GCE instance names are DNS-label-safe; reject anything else before it is
# interpolated into the Cloud Logging filter string below.
if [[ ! "${VM_NAME}" =~ ^[a-z]([-a-z0-9]{0,61}[a-z0-9])?$ ]]; then
  echo "[ERROR]: invalid VM name: '${VM_NAME}'" >&2
  echo "Troubleshooting:" >&2
  echo " • GCE instance names are lowercase letters, digits and hyphens, starting with a letter (max 63 characters)" >&2
  echo " • slidr names its VMs slidr-run-<timestamp>-<random>; copy the name printed when the run launched" >&2
  echo " • List the running VMs: gcloud compute instances list --project=${PROJECT_HINT} --filter='name~slidr-run'" >&2
  exit 1
fi

# The project is not hardcoded: a checkout is not tied to one GCP project, and baking an ID
# into a tracked file publishes it to everyone who clones the repo.
if [[ -z "${PROJECT}" ]]; then
  echo "[ERROR]: could not determine which GCP project to read logs from" >&2
  echo "Troubleshooting:" >&2
  echo " • Pass it as the second argument: ./watch_run.sh ${VM_NAME} <PROJECT>" >&2
  echo " • Or export it: SLIDR_PROJECT=my-project-id ./watch_run.sh ${VM_NAME}" >&2
  echo " • Or set a default for gcloud: gcloud config set project my-project-id" >&2
  echo " • Runs launched by \`./slidr --gcp\` record their project on line 2 of ${LAST_RUN}; an older single-line file has no project to read" >&2
  exit 1
fi

FILTER="resource.type=gce_instance \
  AND labels.\"compute.googleapis.com/resource_name\"=\"${VM_NAME}\" \
  AND logName=\"projects/${PROJECT}/logs/pipeline_runtime\""

echo "=== ${VM_NAME} ===" >&2
echo "" >&2

# Dump all history in chronological order, then stream live entries.
# The tail command (with default --freshness=1d) also replays recent history,
# so we track the last-seen timestamp to avoid printing duplicates.

LAST_TS=""
READ_ERR="$(mktemp)"
while IFS=$'\t' read -r ts msg; do
  echo "${msg}"
  LAST_TS="${ts}"
done < <(
  gcloud logging read "${FILTER}" \
    --project="${PROJECT}" \
    --format='csv[no-heading,separator="\t"](timestamp,jsonPayload.message)' \
    --order=asc 2>"${READ_ERR}" || true
)
if [[ -s "${READ_ERR}" ]]; then
  echo "[WARNING]: failed to fetch historical logs:" >&2
  cat "${READ_ERR}" >&2
  echo " • Check you are authenticated: gcloud auth print-access-token" >&2
  echo " • Check you can read logs in the project: gcloud logging read '' --limit=1 --project=${PROJECT}" >&2
  echo " • The run itself is unaffected -- this only means the log history could not be replayed here" >&2
  echo " • The live tail below may still work; the full log is also uploaded with the run's outputs as log/runtime.log" >&2
fi
rm -f "${READ_ERR}"

echo "" >&2
echo "--- live (Ctrl+C to detach) ---" >&2
echo "" >&2

# Build a time-bounded tail filter to skip entries already shown above.
TAIL_FILTER="${FILTER}"
if [[ -n "${LAST_TS}" ]]; then
  TAIL_FILTER="${TAIL_FILTER} AND timestamp>\"${LAST_TS}\""
fi

gcloud beta logging tail "${TAIL_FILTER}" \
  --format='value(jsonPayload.message)' \
  --project="${PROJECT}" || {
  echo "" >&2
  echo "[ERROR]: could not stream live logs (gcloud's own error is above)" >&2
  echo "Troubleshooting:" >&2
  echo " • \`logging tail\` needs the beta components: gcloud components install beta" >&2
  echo " • It also needs the grpcio Python package: pip install grpcio" >&2
  echo " • Detaching does not stop the run; re-attach at any time with ./watch_run.sh ${VM_NAME}" >&2
  echo " • As a fallback, poll the history instead: gcloud logging read \"${FILTER}\" --project=${PROJECT} --order=asc" >&2
  exit 1
}
