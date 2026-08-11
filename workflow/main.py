# written by Daniel Kabanovsky at the Broad Institute with ❤️ for 🧬

import sys
import time
import atexit
import pandas as pd

from datetime import datetime
from pipeline.helpers import (
    need_run_module,
    log_write,
    log_ts,
    log_detail,
    log_run_started,
    notify_complete,
    log_skipped,
    run_banner,
    run_relative,
    format_duration
)
from pipeline.pipeline import (
    stage_input_data,
    upload_outputs,
    upload_diagnostics,
    run_mkfastq,
    write_metadata_to_file,
    create_samplesheet,
    run_count,
    run_cellbender,
    run_spatial_analysis,
    run_spatial_positioning,
    run_takara_spatial_profiling
)
from config import args, cfg


# unpack globals from config dictionary
START_TIME = cfg['start_time']
ROOT_PATH = cfg['root_path']
OUTPUT_PATH = cfg['output_path']
MKFASTQ_OUTS = cfg['mkfastq_outs']
COUNT_OUTS = cfg['count_outs']
CELLBENDER_OUTS = cfg['cellbender_outs']
SPATIAL_COUNT_OUTS = cfg['spatial_count_outs']
SPATIAL_ANALYSIS_OUTS = cfg['spatial_analysis_outs']
SUMMARY_PATH = cfg['summary_path']
SAMPLESHEET_PATH = cfg['samplesheet_path']
SUMMARY_LOG = cfg['summary_log']
INPUT_BUCKET = cfg['input_bucket']
FASTQ_INPUT = cfg['fastq_input']
STAGE_FASTQ_INPUT = cfg['stage_fastq_input']
BCL_ID = cfg['bcl_id']


# --------------------------------------------------------------------------------------- #
#                               load experimental metadata                                #
# --------------------------------------------------------------------------------------- #

# Preserve this run's logs however it ends. Registered here, before the first stage can fail, because
# a --gcp VM self-deletes as soon as this process exits: without it a crash takes runtime.log and the
# per-stage tool logs down with the disk, and the "see log/<stage>.log for the tool's own output"
# that job_failure prints points at a file nobody can ever read. atexit rather than a call at the
# bottom of this script, since the failing paths leave through sys.exit and never reach it.
atexit.register(upload_diagnostics)

# save sample metadata locally
write_metadata_to_file()

# One Slack "started" alert per run, rather than the one-per-stage stream the old "<stage> started"
# messages produced.
#
# Ordered after write_metadata_to_file, not before it: slack_message resolves the recipient by reading
# the `Email` column out of metadata_summary.csv at send time, and that file is what
# write_metadata_to_file writes. Announcing the run first meant the alert either found no file at all
# -- a first run into a fresh output directory, where it gave up with "User email could not be
# extracted, so Slack alerts have been disabled" and no start alert was ever delivered to anyone -- or,
# re-running into a directory that already had one, found the *previous* run's metadata and delivered
# this run's start alert to whoever that run belonged to. The second case is the reason for the
# ordering: a corrected `Email` could not take effect until after the one message it was meant to fix.
log_run_started(f"run {BCL_ID} started")

# generate samplesheet from metadata
samplesheet_paths = create_samplesheet()

# first samplesheet is the one that contains all sample names
# the rest are auxiliary ones for merging extra RNA/spatial barcodes
samplesheet_path = samplesheet_paths[0] if samplesheet_paths else ""

# read metadata file
metadata_df = pd.read_csv(SUMMARY_PATH)

# Flex chemistry uses cellranger multi (its own cell calling) and the Takara/Trekker spatial path,
# so the standard CellBender and spatial-barcode-counting stages do not apply to Flex runs
is_flex = 'Flex' in metadata_df['Chemistry'].tolist()

samples = [str(row['Sample Name']) for _, row in metadata_df.iterrows()]
chemistry = str(metadata_df['Chemistry'].values[0])
log_ts(f"metadata loaded · {len(samples)} sample(s), {chemistry} chemistry",
       [f"samples → {', '.join(samples)}",
        f"summary → {run_relative(SUMMARY_PATH)}"])
if samplesheet_paths:
    log_detail([f"sheets  → {run_relative(samplesheet_paths[0])}"]
               + [f"{'':<10}{run_relative(extra)}" for extra in samplesheet_paths[1:]])

log_write("  Samples:", SUMMARY_LOG, terminal=False)
for sample in samples:
    log_write(f"   • {sample}", SUMMARY_LOG, terminal=False)


# --------------------------------------------------------------------------------------- #
#                                  run cellranger mkfastq                                 #
# --------------------------------------------------------------------------------------- #


# Whether the input is BCLs to demultiplex or already-demultiplexed FASTQs is stated explicitly by
# the presence of --fastqs, never inferred: --fastqs means skip mkfastq, no --fastqs means the
# input is BCLs. That makes run_mkfastq's BCL run-folder check a real validation with one meaning --
# "the BCLs you asked me to demultiplex are not usable" -- so it can hard-fail instead of being
# reinterpreted as "these must have been FASTQs after all".
#
# A staged bare --fastqs skips mkfastq but still has to fetch its reads, which config.py could not do
# (they are a several-hundred-GB download, and it only resolves paths). Do it here, before the first
# stage that reads them, and apply the checks config.py defers for exactly this case.
if STAGE_FASTQ_INPUT:
    stage_input_data([BCL_ID])
    if not FASTQ_INPUT.is_dir() or not any(FASTQ_INPUT.rglob('*.fastq.gz')):
        log_write(f"[ERROR]: no .fastq.gz files were staged to the FASTQ directory for --fastqs: {FASTQ_INPUT}")
        log_write("Troubleshooting:")
        log_write(f" • --fastqs was given without a directory, so the reads are expected in {INPUT_BUCKET}/{BCL_ID}")
        log_write(f" • Check what that prefix holds: `gcloud storage ls {INPUT_BUCKET}/{BCL_ID}`")
        log_write(" • Check the FASTQs are gzipped (cellranger requires .fastq.gz, not plain .fastq)")
        log_write(" • Drop --fastqs if that folder holds BCLs to be demultiplexed rather than FASTQs")
        log_write(" • Or name a directory already on this machine instead: `--fastqs /path/to/fastqs`")
        sys.exit(1)

# check if cellranger mkfastq needs to be run
if not args.fastqs and (args.run_all or args.mkfastq):
    if need_run_module("mkfastq", metadata_df) or args.force:
        run_mkfastq()
    else:
        log_skipped("mkfastq", f"outputs already present in {run_relative(MKFASTQ_OUTS)} (use --force to regenerate)")


# --------------------------------------------------------------------------------------- #
#                                  run cellranger count                                   #
# --------------------------------------------------------------------------------------- #


# check if cellranger count needs to be run
if args.count or args.run_all:
    if need_run_module("count", metadata_df) or args.force:
        run_count()
    else:
        log_skipped("count", f"outputs already present in {run_relative(COUNT_OUTS)} (use --force to regenerate)")


# --------------------------------------------------------------------------------------- #
#                                     run cellbender                                      #
# --------------------------------------------------------------------------------------- #


# check if cellbender needs to be run
if args.cellbender or (args.run_all and not args.no_cellbender):
    if is_flex:
        log_skipped("cellbender", "not applicable to Flex chemistry (cellranger multi handles cell calling)")
    elif need_run_module("cellbender", metadata_df) or args.force:
        run_cellbender()
    else:
        log_skipped("cellbender", f"outputs already present in {run_relative(CELLBENDER_OUTS)} (use --force to regenerate)")


# --------------------------------------------------------------------------------------- #
#                             run spatial positioning script                              #
# --------------------------------------------------------------------------------------- #

# check if spatial barcodes count needs to be run
if args.spatial_count or args.run_all:
    if is_flex:
        log_skipped("spatial-count", "not applicable to Flex chemistry (spatial barcodes go through the Trekker path in spatial analysis)")
    elif need_run_module("spatial_positioning", metadata_df) or args.force:
        run_spatial_positioning()
    else:
        log_skipped("spatial-count", f"outputs already present in {run_relative(SPATIAL_COUNT_OUTS)} (use --force to regenerate)")


# --------------------------------------------------------------------------------------- #
#                               run spatial analysis script                               #
# --------------------------------------------------------------------------------------- #


# check if spatial analysis / Flex pipeline needs to be run
if args.spatial_analysis or args.run_all:
    if need_run_module("spatial_analysis", metadata_df) or args.force:
        if is_flex:
            run_takara_spatial_profiling()
        else:
            run_spatial_analysis()
    else:
        log_skipped("spatial-analysis", f"outputs already present in {run_relative(SPATIAL_ANALYSIS_OUTS)} (use --force to regenerate)")


# --------------------------------------------------------------------------------------- #
#                                  notify user and exit                                   #
# --------------------------------------------------------------------------------------- #

# Closing banner. The Slack completion alert is sent explicitly rather than by log_write spotting a
# '[DONE]' substring, so the marker does not have to appear inside the banner text itself.
elapsed = format_duration(time.time() - START_TIME)
log_write([
    "",
    run_banner(f"✨ Analysis complete in {elapsed}! ✨"),
    f"  Outputs       {OUTPUT_PATH}",
])
notify_complete(f"Analysis complete in {elapsed}! ✨")

# upload outputs to GCS, into a folder named for this run rather than over an existing one
upload_outputs()

# end timer
end_time = time.time()
elapsed_time = end_time - START_TIME
hours, remainder = divmod(elapsed_time, 3600)
minutes, seconds = divmod(remainder, 60)
log_write(f"  Total runtime:        {int(hours):02d}:{int(minutes):02d}:{seconds:02.0f}", SUMMARY_LOG, terminal=False)
log_write(f"  Result:               SUCCESS", SUMMARY_LOG, terminal=False)


