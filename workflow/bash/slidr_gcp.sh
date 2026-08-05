#!/usr/bin/env bash
set -euo pipefail

mkdir -p /var/log /etc/google-cloud-ops-agent
exec > >(tee /var/log/runtime.log /dev/ttyS0) 2>&1
set -x

META="http://metadata.google.internal/computeMetadata/v1/instance/attributes"
meta() { curl -sf -H "Metadata-Flavor: Google" "${META}/$1"; }

# Self-delete this VM on exit (success or failure) so the pipeline is never
# dependent on the caller's process staying alive.
_self_delete() {
  local name zone project
  name=$(curl -sf "http://metadata.google.internal/computeMetadata/v1/instance/name" \
    -H "Metadata-Flavor: Google")
  zone=$(curl -sf "http://metadata.google.internal/computeMetadata/v1/instance/zone" \
    -H "Metadata-Flavor: Google" | sed 's|.*/||')
  project=$(curl -sf "http://metadata.google.internal/computeMetadata/v1/project/project-id" \
    -H "Metadata-Flavor: Google")
  gcloud compute instances delete "${name}" --zone="${zone}" --project="${project}" --quiet || true
}
trap _self_delete EXIT

# die MESSAGE [HINT...] -- report a fatal error, then list each remaining argument as a
# troubleshooting bullet. Everything this script prints goes to /var/log/runtime.log and the serial
# console, which (because the VM self-deletes) is often the only record of why a run failed -- so a
# failure has to explain itself here rather than leaving the user with a vanished instance.
die() {
    echo "[ERROR]: $1" >&2
    shift
    if [[ $# -gt 0 ]]; then
        echo "Troubleshooting:" >&2
        for hint in "$@"; do echo " • ${hint}" >&2; done
    fi
    exit 1
}

INPUT=$(meta input-gcs) || die \
    "could not read the \`input-gcs\` instance metadata value" \
    "This VM was not created by ./slidr, or was created without --metadata input-gcs=..." \
    "Launch the run with \`./slidr --bcl <BCL_ID> --gcp\` rather than creating the VM by hand"
BCL_ID=$(meta bcl-id) || die \
    "could not read the \`bcl-id\` instance metadata value" \
    "This VM was not created by ./slidr, or was created without --metadata bcl-id=..." \
    "Launch the run with \`./slidr --bcl <BCL_ID> --gcp\` rather than creating the VM by hand"
# `./slidr` (the producer) builds this metadata value with `printf '%q'` per
# argument before embedding it, so every element is already shell-escaped -
# eval here just reverses that encoding to rebuild the array, it does not
# execute arbitrary content from the args themselves. If the producer is ever
# changed to serialize slidr-args differently, this must change too.
SLIDR_ARGS_RAW=$(meta slidr-args) || die \
    "could not read the \`slidr-args\` instance metadata value" \
    "This VM was not created by ./slidr, or was created without --metadata slidr-args=..." \
    "Launch the run with \`./slidr --bcl <BCL_ID> --gcp\` rather than creating the VM by hand"
eval "SLIDR_ARGS=(${SLIDR_ARGS_RAW})" || die \
    "could not decode the pipeline flags from instance metadata: ${SLIDR_ARGS_RAW}" \
    "./slidr shell-quotes each flag with printf '%q' before embedding it; this value is not in that form" \
    "A ',' or '=' in a flag value truncates the metadata -- re-launch without those characters in --fastqs"

THRESHOLD=10      # CPU % below which VM is considered idle
WAIT_MINUTES=30   # How long it must stay idle before shutting down
INTERVAL=60       # Check every 60 seconds

mkdir -p "/pipeline/out/${BCL_ID}/log"
touch "/pipeline/out/${BCL_ID}/log/runtime.log"
chmod 644 "/pipeline/out/${BCL_ID}/log/runtime.log"

sudo tee /etc/google-cloud-ops-agent/config.yaml > /dev/null <<EOF
logging:
  receivers:
    pipeline_runtime:
      type: files
      include_paths:
        - /pipeline/out/${BCL_ID}/log/runtime.log
      record_log_file_path: true
  service:
    pipelines:
      pipeline_pipeline:
        receivers: [pipeline_runtime]

EOF

sudo systemctl restart google-cloud-ops-agent

RUNNER_USER="runner"

# --- hand off the actual pipeline run to the unprivileged runner user ---
# Everything above this point (ops-agent setup, logging) legitimately needs
# root. Everything below (cloning the repo, staging data, running the
# pipeline) does not, and previously ran as root only because `su - $USER`
# with no `-c`/stdin is a no-op under a non-interactive startup script.
# Pre-create and hand over ownership of the directories the runner-context
# work needs to write to, then drop privileges for the rest of the run.
# The sequencing data no longer gets its own directory here: the pipeline stages it inside its own run
# directory (/pipeline/out/${BCL_ID}/data), the way it already stages the reference, pucks and barcodes.
mkdir -p /slidr /pipeline/out /pipeline/software /pipeline/pucks /pipeline/reference /pipeline/barcodes
chown -R "${RUNNER_USER}:${RUNNER_USER}" /slidr /pipeline

echo "Launching slidr as ${RUNNER_USER}..."
runuser -u "${RUNNER_USER}" -- env \
  HOME="/home/${RUNNER_USER}" \
  INPUT="${INPUT}" \
  BCL_ID="${BCL_ID}" \
  PATH="/home/${RUNNER_USER}/.local/bin:/opt/miniforge3/bin:/home/${RUNNER_USER}/.juliaup/bin:${PATH}" \
  bash -s -- "${SLIDR_ARGS[@]}" <<'RUNNER_SCRIPT'
set -euo pipefail

# Same helper as the privileged half of this script. Repeated rather than shared because this block
# is a separate shell fed on stdin, and it matters more here: every step below can fail, the VM
# self-deletes afterwards, and this log is all the user is left with.
die() {
    echo "[ERROR]: $1" >&2
    shift
    if [[ $# -gt 0 ]]; then
        echo "Troubleshooting:" >&2
        for hint in "$@"; do echo " • ${hint}" >&2; done
    fi
    exit 1
}

echo "Cloning the slidr git repository..."
git clone -b stable https://github.com/kabanovskyd/slidr.git /slidr || die \
    "could not clone the slidr repository" \
    "The VM pulls the \`stable\` branch at boot, so it needs outbound access to github.com" \
    "Check the VM's network (vpc1 / my-subnet-central) allows egress, via Cloud NAT or an external IP" \
    "Check the \`stable\` branch still exists in https://github.com/kabanovskyd/slidr"
cd /slidr

# The sequencing data is deliberately NOT copied here. Which run folders this run needs is stated by
# the metadata -- a split-BCL sample merges reads from further BCLs via `Merge RNA/Spatial From BCL` --
# and nothing outside the pipeline reads that metadata, so this script could only ever fetch the one
# BCL ID it was handed. pipeline.stage_input_data does the whole job instead, out of the same ${INPUT}
# prefix (it is `paths.input_path` in the config downloaded just below), into the run's own directory,
# and only if a stage actually needs the reads.
#
# What is left here is the bootstrap that must happen before the pipeline exists: its own config, and
# the key it reads the metadata sheet with.
echo "Staging config and auth key..."
mkdir -p /slidr/config
gcloud storage cp "${INPUT}/config.yaml" /slidr/config/config.yaml || die \
    "could not download ${INPUT}/config.yaml" \
    "config.yaml must be uploaded to the \`paths.input_path\` prefix by hand before launching" \
    "Check the object exists: \`gcloud storage ls ${INPUT}/\`" \
    "Check \`paths.input_path\` names the prefix holding config.yaml, not the bucket root" \
    "Check the slidr-runner service account has Storage Object Viewer on that bucket"
gcloud storage cp "${INPUT}/auth_key.json" /slidr/auth_key.json || die \
    "could not download ${INPUT}/auth_key.json" \
    "auth_key.json must be uploaded to the \`paths.input_path\` prefix by hand before launching" \
    "This is the Google service-account key the pipeline reads its metadata sheet with" \
    "Check the object exists: \`gcloud storage ls ${INPUT}/\`"
chmod 600 /slidr/auth_key.json

echo "Launching slidr..."
uv sync || die \
    "\`uv sync\` failed, so the pipeline's Python dependencies are not installed" \
    "uv's own error is above -- it names the real cause" \
    "This usually means the VM has no outbound access to PyPI; check Cloud NAT or the external IP" \
    "Check the boot disk has free space: relaunch with a larger --disk"
uv run python workflow/main.py --bcl "${BCL_ID}" "$@" || die \
    "the pipeline exited with an error" \
    "The pipeline's own error and troubleshooting notes are above, and in /pipeline/out/${BCL_ID}/log/runtime.log" \
    "Per-stage tool output is in /pipeline/out/${BCL_ID}/log/ (mkfastq.log, count.log, cellbender.log, ...)" \
    "This VM self-deletes shortly, so copy anything you need now -- the log is also streamed to Cloud Logging and readable with ./watch_run.sh" \
    "Outputs are only uploaded to \`settings.output_bucket\` on success, so nothing was written there"
RUNNER_SCRIPT

# signal completion
echo "done"

idle_count=0
required_count=$(( WAIT_MINUTES * 60 / INTERVAL ))

while true; do
  sleep "$INTERVAL"
  # `|| true` keeps a grep miss (pipeline exit non-zero under pipefail) from aborting the script
  cpu=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d. -f1 || true)
  # default to "busy" (100) if the parse yields nothing OR a non-numeric value, so an unexpected
  # top output format can't make `[ "$cpu" -lt ... ]` error out and (under set -e) kill the script
  cpu=${cpu:-100}
  [[ "$cpu" =~ ^[0-9]+$ ]] || cpu=100
  if [ "$cpu" -lt "$THRESHOLD" ]; then
    (( idle_count++ )) || true
    echo "Idle count: $idle_count / $required_count (CPU: ${cpu}%)"
  else
    idle_count=0
  fi
  if [ "$idle_count" -ge "$required_count" ]; then
    echo "Idle threshold reached; deleting VM"
    # break out and let the EXIT trap's _self_delete run to completion (a real
    # `gcloud compute instances delete`), rather than `sudo shutdown -h now`, which merely halts
    # the VM and races the trap — potentially leaving a stopped-but-not-deleted (still billable)
    # instance
    break
  fi
done
