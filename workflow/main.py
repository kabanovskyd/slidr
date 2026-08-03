# written by Daniel Kabanovsky at the Broad Institute with ❤️ for 🧬

import subprocess
import sys
import time
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
    format_duration,
    resolve_bcl_dir,
    validate_bcl_dir
)
from pipeline.pipeline import (
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
INPUT_PATH = cfg['input_path']
MKFASTQ_OUTS = cfg['mkfastq_outs']
COUNT_OUTS = cfg['count_outs']
CELLBENDER_OUTS = cfg['cellbender_outs']
SPATIAL_COUNT_OUTS = cfg['spatial_count_outs']
SPATIAL_ANALYSIS_OUTS = cfg['spatial_analysis_outs']
SUMMARY_PATH = cfg['summary_path']
SAMPLESHEET_PATH = cfg['samplesheet_path']
SUMMARY_LOG = cfg['summary_log']
OUTPUT_BUCKET = cfg['output_bucket']
BCL_ID = cfg['bcl_id']


# --------------------------------------------------------------------------------------- #
#                               load experimental metadata                                #
# --------------------------------------------------------------------------------------- #

# one Slack "started" alert per run, rather than the one-per-stage stream the old
# "<stage> started" messages produced
log_run_started(f"run {BCL_ID} started")

# save sample metadata locally
write_metadata_to_file()

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


# Check if cellranger mkfastq needs to be run.
#
# Whether the input is BCLs to demultiplex or already-demultiplexed FASTQs is stated explicitly by
# the presence of --fastqs, never inferred: --fastqs means skip mkfastq, no --fastqs means the
# input is BCLs. That makes the BCL run-folder check below a real validation with one meaning --
# "the BCLs you asked me to demultiplex are not usable" -- so it can hard-fail instead of being
# reinterpreted as "these must have been FASTQs after all".
if not args.fastqs and (args.run_all or args.mkfastq):
    bcl_dir = resolve_bcl_dir(BCL_ID)
    if not validate_bcl_dir(bcl_dir):
        log_write(f"[ERROR]: {bcl_dir} is not a usable Illumina BCL run directory (see the warnings above for what is missing)")
        log_write("Troubleshooting:")
        log_write(f" • Check that `input_path` points at the directory containing the {BCL_ID} run folder")
        log_write(" • Check that the BCL run finished transferring (RunInfo.xml, RunParameters.xml and Data/Intensities/BaseCalls must all be present)")
        log_write(" • If these reads are already demultiplexed, pass the FASTQ directory with --fastqs instead")
        sys.exit(1)
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

# upload outputs to GCS
if OUTPUT_BUCKET is not None:
    if not OUTPUT_BUCKET.startswith('gs://'):
        OUTPUT_BUCKET = 'gs://' + OUTPUT_BUCKET
    log_write(f'  Uploading the results to GCP: {OUTPUT_BUCKET}... ', terminator='')
    try:
        subprocess.run(
            ['gcloud', 'storage', 'cp', '-r', f'{OUTPUT_PATH}', OUTPUT_BUCKET],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        log_write(' Done.')
    except Exception as exc:
        log_write(f'\n[WARNING]: Could not upload outputs to Google Cloud: {exc}')
        log_write("Troubleshooting:")
        log_write(" • Make sure gcloud is installed: `gcloud --version`")
        log_write(f" • Make sure you are authenticated: `gcloud auth print-access-token`")
        log_write(f" • If not, authenticate with `gcloud auth login`")
        log_write(f" • Make sure the output bucket {OUTPUT_BUCKET.replace('gs://', '').split('/')[0]} exists: `gcloud storage buckets describe {OUTPUT_BUCKET.replace('gs://', '').split('/')[0]}`")
        log_write(f" • Run the upload manually: `gcloud storage cp -r {OUTPUT_BUCKET}`")

# end timer
end_time = time.time()
elapsed_time = end_time - START_TIME
hours, remainder = divmod(elapsed_time, 3600)
minutes, seconds = divmod(remainder, 60)
log_write(f"  Total runtime:        {int(hours):02d}:{int(minutes):02d}:{seconds:02.0f}", SUMMARY_LOG, terminal=False)
log_write(f"  Result:               SUCCESS", SUMMARY_LOG, terminal=False)


