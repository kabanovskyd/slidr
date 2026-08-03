#!/usr/bin/env bash
# Slurm batch payload: stage a run's inputs out of GCS onto cluster storage, then run the pipeline.
#
# This is the Slurm counterpart of workflow/bash/slidr_gcp.sh. It does the same staging work so
# slidr can run on a cluster with no access to the lab filesystem, minus everything GCE-specific:
# there is no VM to create, no instance metadata to read, no privilege drop (Slurm already runs the
# job as the submitting user) and no self-delete. Unlike the GCP script it also does not clone the
# repo -- the job runs inside the checkout it was submitted from.
#
# Submitted by ./slidr --slurm --stage-gcs ; not intended to be run by hand. Inputs arrive as
# environment variables (sbatch exports the submitting environment by default):
#
#   INPUT           GCS prefix holding config.yaml, auth_key.json and the <BCL_ID> data folder
#   BCL_ID          BCL run ID
#   SLIDR_WORKDIR   directory on cluster storage to stage inputs into and write outputs to
#   SLIDR_REPO      path to the slidr checkout to run
#   RESTAGE         optional; set to 1 to re-download data that is already staged
#
# Anything else the run needs is set up outside slidr and inherited through the exported environment:
# gcloud credentials via `gcloud auth login` / `activate-service-account`, and any of config.py's
# PATH_ENV_OVERRIDES (notably SLIDR_SOFTWARE_PATH) exported before submitting.
#
# Pipeline flags are passed as positional arguments.

set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# die MESSAGE [HINT...] -- report a fatal error, then list each remaining argument as a
# troubleshooting bullet. Keeping the hints in the call keeps every failure in this script
# self-explanatory in the batch job's output, which is often the only thing the user ever sees.
die() {
    echo "[ERROR]: $1" >&2
    shift
    if [[ $# -gt 0 ]]; then
        echo "Troubleshooting:" >&2
        for hint in "$@"; do echo " • ${hint}" >&2; done
    fi
    exit 1
}

: "${INPUT:?INPUT (GCS input prefix) must be set}"
: "${BCL_ID:?BCL_ID must be set}"
: "${SLIDR_WORKDIR:?SLIDR_WORKDIR must be set}"
: "${SLIDR_REPO:?SLIDR_REPO must be set}"

log "slidr Slurm run starting"
log "  job id:    ${SLURM_JOB_ID:-<not under slurm>}"
log "  node:      $(hostname)"
log "  repo:      ${SLIDR_REPO}"
log "  workdir:   ${SLIDR_WORKDIR}"
log "  input:     ${INPUT}"
log "  bcl id:    ${BCL_ID}"
log "  args:      $*"

# uv installs to a leaf bindir; re-add it here because a batch job may start from a login shell
# profile that never saw it
export PATH="${HOME}/.local/slidr:${PATH}"

# --------------------------------------------------------------------------------------- #
#                                    check the toolchain                                  #
# --------------------------------------------------------------------------------------- #

command -v gcloud >/dev/null 2>&1 || die \
    "\`gcloud\` is not on PATH, so inputs cannot be staged from GCS" \
    "Many clusters ship the CLI as a module -- load it before submitting (e.g. \`module load google-cloud-sdk\`)" \
    "Otherwise install it: https://cloud.google.com/sdk/docs/install" \
    "A module loaded interactively does not reach the compute node unless it is loaded before \`sbatch\`, or in your shell profile" \
    "PATH seen by this job: ${PATH}"

command -v uv >/dev/null 2>&1 || die \
    "\`uv\` is not on PATH, so the pipeline's Python dependencies cannot be installed" \
    "Run ./slidr once on the submit node to install it into \${HOME}/.local/slidr" \
    "Or install it yourself: https://docs.astral.sh/uv/getting-started/installation/" \
    "This job looked in \${HOME}/.local/slidr, which must be on a filesystem the compute node can see" \
    "PATH seen by this job: ${PATH}"

# Credentials are established outside slidr. A GCE VM gets them from its attached service account; an
# unrelated cluster has none, so the user authenticates once on the submit node and gcloud persists it
# to ~/.config/gcloud -- which compute nodes see on a shared home. slidr deliberately does not take a
# key file itself: `gcloud auth activate-service-account` already does exactly that, persistently, and
# threading a second path for it through the submit script bought nothing.
if ! gcloud auth print-access-token >/dev/null 2>&1; then
    die "no active gcloud credentials on this node" \
        "Authenticate once on the submit node with \`gcloud auth login\` -- it persists to ~/.config/gcloud, which compute nodes see on a shared home" \
        "For a headless cluster, use a service-account key instead: \`gcloud auth activate-service-account --key-file=/path/to/key.json\`" \
        "Check which account gcloud thinks it has: \`gcloud auth list\`" \
        "If your home directory is NOT shared with compute nodes, set CLOUDSDK_CONFIG to a shared location before submitting"
fi
log "gcloud credentials OK ($(gcloud config get-value account 2>/dev/null || echo unknown))"

# --------------------------------------------------------------------------------------- #
#                                     stage the inputs                                    #
# --------------------------------------------------------------------------------------- #

CONFIG_FILE="${SLIDR_WORKDIR}/config.yaml"
AUTH_KEY="${SLIDR_WORKDIR}/auth_key.json"
DATA_DIR="${SLIDR_WORKDIR}/data"
OUT_DIR="${SLIDR_WORKDIR}/outs"

# NB: create the DATA_DIR parent but NOT the ${DATA_DIR}/${BCL_ID} leaf. `gcloud storage cp -r SRC
# DEST` nests SRC under DEST when DEST already exists, which would yield data/<BCL_ID>/<BCL_ID>/...;
# leaving the leaf absent makes the copy land the contents directly where the pipeline expects them.
mkdir -p "$SLIDR_WORKDIR" "$DATA_DIR" "$OUT_DIR" || die \
    "could not create the working directories under ${SLIDR_WORKDIR}" \
    "Check the parent directory exists and is writable from this node" \
    "Point --workdir at scratch you own, e.g. /scratch/\${USER}/slidr/${BCL_ID}" \
    "Check your quota and the filesystem's free space: \`df -h ${SLIDR_WORKDIR}\`"

log "Staging config.yaml and auth_key.json from ${INPUT}"
gcloud storage cp "${INPUT}/config.yaml" "$CONFIG_FILE" \
    || die "could not download ${INPUT}/config.yaml" \
           "Check the object exists: \`gcloud storage ls ${INPUT}/\`" \
           "config.yaml must be uploaded to the \`settings.input_bucket\` prefix by hand before submitting" \
           "Check the active account can read the bucket: \`gcloud storage ls ${INPUT}/\`" \
           "Check \`settings.input_bucket\` names the prefix holding config.yaml, not the bucket root"
gcloud storage cp "${INPUT}/auth_key.json" "$AUTH_KEY" \
    || die "could not download ${INPUT}/auth_key.json" \
           "Check the object exists: \`gcloud storage ls ${INPUT}/\`" \
           "auth_key.json must be uploaded to the \`settings.input_bucket\` prefix by hand before submitting" \
           "This is the Google service-account key the pipeline reads its metadata sheet with" \
           "Check the active account can read the bucket: \`gcloud storage ls ${INPUT}/\`"
chmod 600 "$AUTH_KEY"

# Sequencing data is the expensive part of staging (routinely hundreds of GB), so skip it when a
# previous job in this workdir already brought it down. RESTAGE=1 forces a re-download.
if [[ -d "${DATA_DIR}/${BCL_ID}" ]] && [[ -n "$(ls -A "${DATA_DIR}/${BCL_ID}" 2>/dev/null)" ]] \
   && [[ "${RESTAGE:-}" != "1" ]]; then
    log "Sequencing data already staged at ${DATA_DIR}/${BCL_ID} (set RESTAGE=1 to re-download)"
else
    log "Staging sequencing data from ${INPUT}/${BCL_ID} (this can take a while)"
    rm -rf "${DATA_DIR}/${BCL_ID}"
    gcloud storage cp -r "${INPUT}/${BCL_ID}" "${DATA_DIR}/${BCL_ID}" \
        || die "could not download the sequencing data from ${INPUT}/${BCL_ID}" \
               "Check the folder exists: \`gcloud storage ls ${INPUT}/\`" \
               "Check the BCL ID matches the folder name in the bucket exactly, including case" \
               "A BCL run is routinely hundreds of GB -- check the free space and your quota: \`df -h ${DATA_DIR}\`" \
               "A partial download was left behind; it is removed and retried automatically on the next submission"
    log "Sequencing data staged"
fi

# --------------------------------------------------------------------------------------- #
#                                      run the pipeline                                   #
# --------------------------------------------------------------------------------------- #

# Point the bucket's config.yaml at this cluster's locations. The config was written for whatever
# machine produced it and cannot know this workdir, so these four host-specific paths are overridden
# in the environment (config.py's PATH_ENV_OVERRIDES) rather than by rewriting the YAML.
export SLIDR_INPUT_PATH="$DATA_DIR"
export SLIDR_OUTPUT_PATH="$OUT_DIR"
export SLIDR_AUTH_KEY_PATH="$AUTH_KEY"
# SLIDR_SOFTWARE_PATH is deliberately not set here. config.py already honours it from the environment
# (PATH_ENV_OVERRIDES), and sbatch exports the submitting environment by default, so
# `export SLIDR_SOFTWARE_PATH=/opt/cellranger` before submitting reaches the job untouched -- which is
# why the old --software flag was removed rather than reimplemented here.

cd "$SLIDR_REPO" || die \
    "could not enter the slidr checkout at ${SLIDR_REPO}" \
    "The job runs from the checkout it was submitted from, so that path must be visible to the compute node" \
    "Check the directory still exists and has not been moved or deleted since submission"

log "Syncing Python dependencies"
uv sync || die \
    "\`uv sync\` failed, so the pipeline's Python dependencies are not installed" \
    "uv's own error is above -- it names the real cause" \
    "This usually means no network access from the compute node; run \`uv sync\` on the submit node first so the cache is warm" \
    "Check there is free space for the virtualenv: \`df -h ${SLIDR_REPO}\`" \
    "If the lockfile is stale, refresh it on the submit node with \`uv lock\`"

log "Launching the pipeline"
# --stage-gcs makes the pipeline pull the reference genome, puck maps and raw barcodes from the
# gs:// locations in the staged config; --config points at the staged copy so concurrent jobs never
# fight over the repo's own config/config.yaml.
uv run python workflow/main.py \
    --bcl "$BCL_ID" \
    --config "$CONFIG_FILE" \
    --stage-gcs \
    "$@" \
    || die "the pipeline exited with an error" \
           "The pipeline's own error and troubleshooting notes are above, and in ${OUT_DIR}/${BCL_ID}/log/runtime.log" \
           "Per-stage tool output is in ${OUT_DIR}/${BCL_ID}/log/ (mkfastq.log, count.log, cellbender.log, ...)" \
           "If the job was killed rather than failing, check the scheduler's accounting for an OOM or time-limit kill: \`sacct -j ${SLURM_JOB_ID:-<job id>} --format=JobID,State,ExitCode,MaxRSS,Elapsed\`" \
           "Staged inputs are kept in ${SLIDR_WORKDIR}, so a resubmission does not re-download them"

log "slidr Slurm run finished; outputs in ${OUT_DIR}"
