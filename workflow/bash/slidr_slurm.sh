#!/usr/bin/env bash
# Slurm batch payload: stage a run's inputs out of GCS onto cluster storage, then run the pipeline.
#
# This is the Slurm counterpart of workflow/bash/slidr_gcp.sh. It does the same staging work so
# slidr can run on a cluster with no access to the lab filesystem, minus everything GCE-specific:
# there is no VM to create, no instance metadata to read, no privilege drop (Slurm already runs the
# job as the submitting user) and no self-delete. Unlike the GCP script it also does not clone the
# repo -- the job runs inside the checkout it was submitted from.
#
# It also does NOT stage config.yaml, which is the other place it diverges from the GCP script. A
# Slurm job runs from a real checkout on a cluster the submitter administers, so config/config.yaml is
# already right there and is the file they just edited; downloading a copy written for a different
# machine on top of it only invited the two to disagree. The GCP script still needs its own copy
# because a fresh VM has no checkout to read. Consequences: the bucket needs no config.yaml for a
# --slurm --stage-gcs run, and settings the submitter changes locally take effect on the next
# submission with nothing to re-upload.
#
# Submitted by ./slidr --slurm --stage-gcs ; not intended to be run by hand. Inputs arrive as
# environment variables (sbatch exports the submitting environment by default):
#
#   INPUT           GCS prefix holding the <BCL_ID> data folder (and optionally auth_key.json)
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

# The service-account key is staged when the bucket has one, but its absence is not fatal: the local
# config this job runs with may already point `auth_key_path` at a key on cluster storage or at a gs://
# object, and a run whose metadata is a local .tsv/.csv needs no key at all. Requiring one here would
# force those users to upload a key nothing reads. When the download does succeed it wins, on the
# grounds that a key placed in the run's own input prefix was put there for this job.
log "Staging auth_key.json from ${INPUT} (optional)"
if gcloud storage cp "${INPUT}/auth_key.json" "$AUTH_KEY" 2>/dev/null; then
    chmod 600 "$AUTH_KEY"
    export SLIDR_AUTH_KEY_PATH="$AUTH_KEY"
    log "  Staged the service-account key to ${AUTH_KEY}"
else
    rm -f "$AUTH_KEY"
    log "  No auth_key.json at ${INPUT} -- using \`paths.auth_key_path\` from the local config"
    log "  (a Google Sheet metadata source needs a key there; a local .tsv/.csv needs none)"
fi

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

# Point the run at the locations this job staged into. The checkout's config.yaml cannot know the
# workdir sbatch was given, so these host-specific paths are overridden in the environment (config.py's
# PATH_ENV_OVERRIDES) rather than by rewriting the YAML. Everything else -- threads, memory, buckets,
# reference/puck/barcode locations, metadata source -- comes from config/config.yaml as committed.
# SLIDR_AUTH_KEY_PATH is exported above, but only if a key was actually staged.
export SLIDR_INPUT_PATH="$DATA_DIR"
export SLIDR_OUTPUT_PATH="$OUT_DIR"
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
# --stage-gcs makes the pipeline pull the reference genome, puck maps, raw barcodes and -- when
# `software_path` is a gs:// location -- the software, from the locations named in the checkout's own
# config/config.yaml. No --config: this job runs with the repo's config, not a staged copy of one.
uv run python workflow/main.py \
    --bcl "$BCL_ID" \
    --stage-gcs \
    "$@" \
    || die "the pipeline exited with an error" \
           "The pipeline's own error and troubleshooting notes are above, and in ${OUT_DIR}/${BCL_ID}/log/runtime.log" \
           "Per-stage tool output is in ${OUT_DIR}/${BCL_ID}/log/ (mkfastq.log, count.log, cellbender.log, ...)" \
           "If the job was killed rather than failing, check the scheduler's accounting for an OOM or time-limit kill: \`sacct -j ${SLURM_JOB_ID:-<job id>} --format=JobID,State,ExitCode,MaxRSS,Elapsed\`" \
           "Staged inputs are kept in ${SLIDR_WORKDIR}, so a resubmission does not re-download them"

log "slidr Slurm run finished; outputs in ${OUT_DIR}"
