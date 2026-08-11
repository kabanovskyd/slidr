import os
import sys
import re
import csv
import subprocess
import signal
import types
import shutil
import pynvml
import zipfile

import pandas as pd
import xml.etree.ElementTree as ET

from datetime import datetime
from pathlib import Path
from rich.console import Console

sys.path.append(str(Path(__file__).parent))
from config import (
    args,
    cfg,
    gcs_uri,
    is_int,
    is_number,
    bool_value_hint,
    GCS_METADATA_TIMEOUT,
    GCS_TRANSFER_TIMEOUT,
)

from helpers import (
    log_write,
    log_ts,
    log_detail,
    log_stage_start,
    format_duration,
    run_relative,
    stage_from_gcs,
    run_gcloud_transfer,
    job_crash,
    job_success,
    format_multi_samplesheet,
    ensure_conda_env,
    load_metadata,
    merge_bcls,
    row_declares,
    resolve_bcl_dir,
    validate_bcl_dir,
    test_and_install_software,
    create_tmp_dir,
    parse_cellranger_html,
    missing_fastqs,
    find_sample_fastqs,
    stage_spatial_fastqs,
    declared_lanes,
    FASTQ_NAME_RE,
    retrieve_takara_bead_barcode_file,
    sanitize_path_component,
    split_probe_barcodes,
    resolve_cellranger_version,
    downsample_spatial
)


# unpack globals from config dictionary
ROOT_PATH = cfg['root_path']
OUTPUT_PATH = cfg['output_path']
SCRIPT_PATH = cfg['script_path']
LOG_PATH = cfg['log_path']
INPUT_PATH = cfg['input_path']
FASTQ_INPUT = cfg['fastq_input']
STAGE_GCS = cfg['stage_gcs']
METADATA_PATH = cfg['metadata_path']
TMP_PATH = cfg['tmp_path']
METADATA_SRC = cfg['metadata_src']
RAW_BARCODES_PATH = cfg['raw_barcodes_path']
PUCK_PATH = cfg['puck_path']
REF_PATH = cfg['ref_path']
INPUT_BUCKET = cfg['input_bucket']
OUTPUT_BUCKET = cfg['output_bucket']
OUTPUT_DEST = cfg['output_dest']
SUMMARY_PATH = cfg['summary_path']
SUMMARY_LOG = cfg['summary_log']
REF_GENOME = cfg['ref_genome']
NUM_THREADS = cfg['num_threads']
MEM_SIZE = cfg['mem_size']
MKFASTQ_OUTS = cfg['mkfastq_outs']
COUNT_OUTS = cfg['count_outs']
CELLBENDER_OUTS = cfg['cellbender_outs']
CELLBENDER_CELLS = cfg['cellbender_cells']
CELLBENDER_DROPLETS = cfg['cellbender_droplets']
CELLBENDER_EPOCHS = cfg['cellbender_epochs']
CELLBENDER_RATE = cfg['cellbender_rate']
SPATIAL_DOWNSAMPLING = cfg['spatial_downsampling']
PERCENT_UMI_FILTERING = cfg['percent_umi_filtering']
FLEX_OUTS = cfg['flex_outs']
SPATIAL_COUNT_OUTS = cfg['spatial_count_outs']
SPATIAL_ANALYSIS_OUTS = cfg['spatial_analysis_outs']
SAMPLESHEET_PATH = cfg['samplesheet_path']
GENERATE_BAM = cfg['generate_bam']
EMPTYDROPS_MIN_UMI = cfg['emptydrops_min_umis']
FLEX_PROBE_SET = cfg['flex_probe_set']
FLEX_R1_PATH = cfg['flex_r1_path']
FLEX_R2_PATH = cfg['flex_r2_path']
FLEX_GEX_FASTQS = cfg['flex_gex_fastqs']
BCL_ID = cfg['bcl_id']

# Local directories GCS-hosted resources are staged into when --stage-gcs/--gcp is in effect. They
# live under OUTPUT_PATH so a run is self-contained: on a cluster the pipeline may have no writable
# location other than its own output tree, and keeping them there means nothing is left behind
# outside it.
#   REF_DEST       reference genome, assigned by create_samplesheet (via `global REF_DEST`) and read
#                  by run_count/run_spatial_analysis. Declared here so a reference before assignment
#                  is a clear guarded error rather than a NameError.
#   PUCK_DEST      puck CSVs -- both those downloaded from puck_path and those generated locally
#                  from raw barcodes, which is why it must be a writable local directory and not the
#                  gs:// puck_path itself
#   BARCODES_DEST  raw barcode directories downloaded from raw_barcodes_path
REF_DEST = None
PUCK_DEST = OUTPUT_PATH / "pucks"
BARCODES_DEST = OUTPUT_PATH / "barcodes"

# Directories under output/ that upload_outputs does not send to `settings.output_bucket`.
#
# Both are staging destinations for a --stage-gcs/--gcp run's own *inputs*, not results: `reference`
# is where the `settings.reference_genome` subdirectory of `paths.reference_path` is downloaded
# (assigned to REF_DEST below), and `software` is where a gs:// `paths.software_path` is unpacked
# (helpers.test_and_install_software). They live under output/ only because a staged run needs a
# local, writable place to put them, and uploading them again returns bytes to the bucket they were
# just read out of -- 14 GiB of reference genome and 2.5 GiB of Cellranger on a mouse run, against
# roughly 5 GiB of actual analysis output.
#
# `pucks` and `barcodes` are deliberately NOT here: puck CSVs can be *generated* during a run from
# raw barcodes, so that directory holds real output for the runs that produce it, and both are small.
UPLOAD_SKIP = frozenset({'reference', 'software'})


def check_barcode_validity(
    sample_name: str,
    count_dir: Path,
    chemistry: str,
    threshold: float = 50.0
) -> None:
    """
    Warn when a sample's cell-barcode validity is too low to be explained by anything but the wrong
    `Chemistry`.

    Cellranger reports `Valid Barcodes` -- the share of reads whose cell barcode is in the whitelist
    for the chemistry it was told to assume. A healthy library sits at 85-95%. A chemistry mismatch
    puts it near the rate at which random 16-mers happen to collide with the *other* whitelist, which
    is a few percent, because 3' and 5' v3 use different whitelists.

    Cellranger treats this as a metric, not an error: it exits 0, writes a full output directory, and
    the pipeline carries on. The damage surfaces two stages later and in a form that names nothing
    useful -- with almost no cells called, CellBender finds no empty droplets to learn an ambient
    profile from and dies inside scipy with "ZeroDivisionError: float division by zero". Diagnosing
    that from the traceback means testing barcodes against candidate whitelists by hand.

    So this is a warning, not a failure: the run is almost certainly wrong, but cellranger's own
    output is intact and the operator may have a reason to continue. It says the one thing the
    traceback two stages later will not -- look at the `Chemistry` column.

    Inputs:
     - sample_name: sample being checked, for the message
     - count_dir:   this sample's count output directory, holding metrics_summary.csv
     - chemistry:   chemistry string handed to cellranger, quoted back in the warning
     - threshold:   percent below which the rate is treated as implausible
    Output:
     - none; a missing or unparseable metrics file is silently ignored, since this is a diagnostic
       and must never be the reason a successful run fails
    """

    metrics = Path(count_dir) / 'metrics_summary.csv'
    if not metrics.is_file():
        return

    try:
        with open(metrics, newline='') as handle:
            row = next(csv.DictReader(handle), None)
        if not row or 'Valid Barcodes' not in row:
            return
        # cellranger writes this as a percentage string, e.g. "85.7%"
        rate = float(str(row['Valid Barcodes']).strip().rstrip('%'))
    except (OSError, StopIteration, TypeError, ValueError):
        return

    if rate >= threshold:
        return

    log_write(f"[WARNING]: only {rate:.1f}% of {sample_name}'s reads carry a valid cell barcode, which normally means `Chemistry` is wrong")
    log_write("Troubleshooting:")
    log_write(f" • cellranger was told this library is {chemistry}; a correct chemistry gives 85-95% here")
    log_write(" • 3' and 5' 10x kits use different cell-barcode whitelists, so the wrong one leaves only a few percent matching by chance")
    log_write(f" • Check the `Chemistry` metadata column for {sample_name} against how the library was actually prepared")
    log_write(f" • Full metrics are in {run_relative(metrics)}, and the barcode-rank plot in {run_relative(Path(count_dir) / 'web_summary.html')}")
    log_write(" • Left uncorrected, cellbender fails later with an unrelated-looking error, having found no empty droplets")


def conda_subprocess_env() -> dict:
    """
    This process's environment, with the per-user site directory switched off, for running a tool out
    of a conda environment `ensure_conda_env` just resolved.

    Python puts `~/.local/lib/pythonX.Y/site-packages` on sys.path ahead of the interpreter's own
    site-packages whenever the two version tags match. A conda environment is therefore only as
    isolated as the invoking user's home directory is empty: `<env>/bin/cellbender` is the entry
    point, but `import cellbender` (and `torch`, and `pyro`) can resolve to a completely different
    pip install under ~/.local. The environment slidr carefully built or found is then not the one
    the tool actually runs against, and nothing says so.

    That is not hypothetical. The `--gcp` VM image carries a stale user-site install, so cellbender
    trained all 150 epochs against ~/.local's torch and then failed to save its checkpoint with
    "TypeError: cannot pickle 'weakref.ReferenceType' object" -- a torch/pyro version
    incompatibility that does not exist in the r_env the image also ships. The stage died on the
    assertion that follows, since the posterior step requires the checkpoint that was never written.

    PYTHONNOUSERSITE drops the user directory from sys.path, so the environment resolved is the
    environment used. Applied to every tool slidr launches from a conda env, rather than only the one
    that broke, because the failure mode is silent everywhere else too.

    Output:
     - an environment mapping suitable for passing as subprocess `env=`
    """

    env = os.environ.copy()
    env['PYTHONNOUSERSITE'] = '1'
    return env


def staged_ref_dir() -> Path:
    """
    Directory holding the reference genome for this run: the local copy when staging from GCS, the
    configured reference_path otherwise.
    """

    if not STAGE_GCS:
        return Path(REF_PATH) / REF_GENOME
    if REF_DEST is None:
        log_write("[ERROR]: the reference genome has not been staged yet (REF_DEST unset)")
        log_write("Troubleshooting:")
        log_write(" • This stage needs the reference; run it together with --count/--run-all, which stages it first")
        sys.exit(1)
    return REF_DEST / REF_GENOME


# set up a rich console
console = Console()

# define a global process state and a signal handler for capturing SIGINT/SIGTERM signals
proc = None

def handle_signal(
    signum: int,
    frame: types.FrameType | None
) -> None:
    """
    Gracefully terminate the active subprocess on SIGINT or SIGTERM

    Inputs:
     - signum: signal number received from the OS
     - frame:  current stack frame (required by signal.signal interface, unused)
    """

    if proc is not None:
        log_write(f"[ERROR]: signal {signum} received, terminating subprocess...")
        log_write("Troubleshooting:")
        log_write(" • Partial outputs from the interrupted stage may be left behind; re-run with --force to overwrite them")
        log_write(" • If you did not interrupt this yourself, the job was killed externally -- on a cluster check the scheduler's accounting for an OOM or time-limit kill (e.g. `sacct -j $SLURM_JOB_ID`)")
        # send SIGTERM to subprocess
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            # force kill if process doesn't stop
            proc.kill()
    sys.exit(1)

# set signal handlers to monitor for signals
signal.signal(signal.SIGTERM, handle_signal)
signal.signal(signal.SIGINT, handle_signal)


def write_metadata_to_file() -> None:
    """
    Load sample metadata from Google Sheet or CSV, filter and verify the data, and save the metadata to a CSV
    """

    # log execution
    log_ts("loading sample metadata...")

    # load sample metadata from Google Sheet or CSV
    metadata_df = load_metadata(METADATA_SRC)

    # check if any rows with matching BCL have the "Run" column set to "YES" and exit if not
    mask = (metadata_df['BCL'] == BCL_ID) & (metadata_df['Run'].str.upper().isin(['YES', 'Y']))
    matching_rows = metadata_df[mask]
    if matching_rows.empty:
        if len(metadata_df[metadata_df['BCL'] == BCL_ID].values) > 0:
            log_write(f"[ERROR]: All the runs with the specified BCL have their run status set to `NO` in the Google Sheet: {BCL_ID}")
            log_write(f"Troubleshooting:")
            log_write(f" • Find your runs in the samplesheet and set their value in the `Run` column to `YES`")
            log_write(f" • Metadata source: {METADATA_SRC}")
            sys.exit(1)
        else:
            log_write(f"[ERROR]: No runs matching the provided BCL found in the input table/Google Sheets: {BCL_ID}")
            log_write("Troubleshooting:")
            log_write(f" • Check that the `BCL` column contains a row whose value is exactly `{BCL_ID}` (no trailing spaces, and matching case)")
            log_write(" • Check that --bcl was given the BCL ID as written in the metadata, not the path to the run folder")
            log_write(f" • Metadata source: {METADATA_SRC}")
            sys.exit(1)

    # validate that species and chemistry are the same across all samples
    if len(set(matching_rows['Species'].dropna().tolist())) > 1:
        log_write(f"[ERROR]: More than a single species is specified: ")
        for species in set(matching_rows['Species'].dropna().tolist()):
            log_write(f" • {species}")
        log_write("Troubleshooting:")
        log_write(" • One run uses one reference genome, so all samples in a BCL must share a `Species` value")
        log_write(f" • Correct the `Species` column for the {BCL_ID} rows, or split the mismatched samples into a separate BCL and run them one at a time")
        log_write(f" • Metadata source: {METADATA_SRC}")
        sys.exit(1)

    if len(set(matching_rows['Chemistry'].dropna().tolist())) > 1:
        log_write(f"[ERROR]: More than a single sequencing chemistry is specified: ")
        for chemistry in set(matching_rows['Chemistry'].dropna().tolist()):
            log_write(f" • {chemistry}")
        log_write("Troubleshooting:")
        log_write(" • Chemistry selects which cellranger invocation the whole run uses, so all samples in a BCL must share a `Chemistry` value")
        log_write(f" • Correct the `Chemistry` column for the {BCL_ID} rows, or set `Run` to `NO` on the mismatched samples and run them separately")
        log_write(f" • Metadata source: {METADATA_SRC}")
        sys.exit(1)

    # write the metadata sheet to file
    matching_rows.to_csv(SUMMARY_PATH, index=False)


def create_samplesheet() -> list[Path]:
    """
    Load metadata from a local CSV file and format according to Cellranger specifications

    Output:
     - A list of paths to generated samplesheet files in CSV format
    """

    # log execution
    log_ts("generating samplesheets...")

    # load the metadata DataFrame
    try:
        metadata_df = pd.read_csv(SUMMARY_PATH)
    except Exception as exc:
        log_write(f"[ERROR]: Could not read {SUMMARY_PATH}: {exc}")
        log_write("Troubleshooting:")
        log_write(" • This file is written by the metadata-loading step at the start of every run; a failure here usually means the run's output directory is not writable or ran out of space")
        log_write(f" • Check the permissions and free space on {SUMMARY_PATH.parent}")
        log_write(" • Delete the file and re-run to have it regenerated from the metadata source")
        sys.exit(1)

    if args.run_all or args.count or args.mkfastq or args.spatial_analysis:
        # match species name with a reference genome folder
        global REF_GENOME
        species = metadata_df['Species'].values[0]

        match species:
            case 'Mouse':
                REF_GENOME = 'refdata-gex-GRCm39-2024-A'
            case 'Human':
                REF_GENOME = 'refdata-gex-GRCh38-2024-A'
            case _:
                if REF_GENOME is None:
                    log_write(f'[ERROR]: Analyzing {species} samples requiring a custom reference genome, but `reference_genome` is not set in the configuration file')
                    log_write("Troubleshooting:")
                    log_write(f" • Set `settings.reference_genome` to the name of the reference directory to use for {species}")
                    log_write(f" • That directory must exist under `paths.reference_path` ({REF_PATH})")
                    log_write(" • `Human` and `Mouse` select a reference automatically; check the `Species` column for a typo if you meant one of those")
                    sys.exit(1)
        
        # stage the reference genome out of GCS when the local filesystem does not hold it
        if STAGE_GCS:
            global REF_DEST
            REF_DEST = OUTPUT_PATH / "reference"
            REF_DEST.mkdir(exist_ok=True)
            # skip the download if a previous run in this output tree already staged it
            if (REF_DEST / REF_GENOME).is_dir():
                log_write(f"  Reference genome {REF_GENOME} already staged at {REF_DEST}")
            else:
                log_write(f"  Staging reference genome {REF_GENOME} from {REF_PATH}... ", terminator="")
                with console.status(f"  Staging reference genome {REF_GENOME}..."):
                    stage_from_gcs(
                        f'{gcs_uri(REF_PATH)}/{REF_GENOME}',
                        REF_DEST,
                        recursive=True,
                        description=f'reference genome {REF_GENOME}'
                    )
                log_write("Done.")

        # set up data for cellranger multi samplesheet
        if metadata_df['Chemistry'].values[0] == 'Flex':
            # when staging, the reference lives in the local copy under REF_DEST rather than at the
            # gs:// reference_path, which cannot be path-joined
            ref_dir = staged_ref_dir()

            # validate flex parameters
            if not Path(FLEX_PROBE_SET).is_file():
                log_write(f'[ERROR]: the Flex probe set {FLEX_PROBE_SET} is not found or is not a file')
                log_write('Troubleshooting:')
                log_write(" • Set `workflow.flex_probe_set` to the probe-set CSV that ships with your Flex kit's cellranger release")
                log_write(" • The probe set is shipped inside the cellranger install, e.g. <cellranger>/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv")
                log_write(" • Check the path for typos and that the file is readable")
                sys.exit(1)
            if not ref_dir.is_dir():
                log_write(f'[ERROR]: the reference genome directory {ref_dir} is not found or is not a directory')
                log_write('Troubleshooting:')
                log_write(f" • Check that `settings.reference_genome` ({REF_GENOME}) names a directory that exists under `paths.reference_path` ({REF_PATH})")
                if STAGE_GCS:
                    log_write(f" • Check the reference exists in the bucket: `gcloud storage ls {gcs_uri(REF_PATH)}/{REF_GENOME}`")
                else:
                    log_write(" • If the reference lives in a bucket rather than on this filesystem, pass --stage-gcs to download it")
                sys.exit(1)
            # Resolve the emptydrops cutoff into a local rather than reassigning the module global.
            # Assigning to EMPTYDROPS_MIN_UMI here would make the name local to this whole function,
            # so the isinstance() read below would raise UnboundLocalError before the default could
            # ever be applied. A local is also the honest scope: nothing outside this function reads
            # the value (unlike REF_GENOME/REF_DEST above, which staged_ref_dir()/run_count() do read
            # and which therefore genuinely need `global`).
            emptydrops_min_umi = EMPTYDROPS_MIN_UMI
            if emptydrops_min_umi is None:
                emptydrops_min_umi = 500
            elif not is_int(emptydrops_min_umi):
                log_write(f'[ERROR]: unrecognized value for the `flex_emptydrops_minimum_umis` field in the configuration file (should be an integer): {emptydrops_min_umi!r}')
                log_write('Troubleshooting:')
                log_write(" • Set `workflow.flex_emptydrops_minimum_umis` to a plain integer (no quotes, no decimal point), e.g. 100")
                log_write(" • Remove the field entirely to use the default of 500")
                hint = bool_value_hint(emptydrops_min_umi, 'workflow.flex_emptydrops_minimum_umis')
                if hint:
                    log_write(hint)
                sys.exit(1)
            # cellranger multi consumes FASTQs, never BCLs, so a Flex run's GEX reads come from
            # --fastqs when it is given and from `input_path` otherwise
            flex_fastq_path = FASTQ_INPUT if FASTQ_INPUT is not None else Path(INPUT_PATH)

            if not isinstance(FLEX_GEX_FASTQS, list):
                log_write(f'[ERROR]: the `flex_gex_fastqs` field of the configuration file should be a list: {FLEX_GEX_FASTQS}')
                log_write('Troubleshooting:')
                log_write(" • Write the field as a YAML list of FASTQ name prefixes, one per line:")
                log_write("     flex_gex_fastqs:")
                log_write("       - fastq_prefix_1")
                log_write("       - fastq_prefix_2")
                log_write(" • A single prefix still needs the leading `- ` so YAML parses it as a list rather than a bare string")
                sys.exit(1)
            for fastq in FLEX_GEX_FASTQS:
                if not isinstance(fastq, str) or not fastq.strip():
                    log_write(f'[ERROR]: each entry in the `flex_gex_fastqs` field of the configuration file should be a non-empty string, but one is: {fastq!r}')
                    log_write('Troubleshooting:')
                    log_write(" • Remove the blank/malformed list entry -- a stray `- ` on its own line parses as an empty value")
                    log_write(" • Each entry is the FASTQ filename prefix cellranger matches on, e.g. `mysample_GEX` for mysample_GEX_S1_L001_R1_001.fastq.gz")
                    sys.exit(1)
                fastq_glob = list(flex_fastq_path.rglob(f'{fastq}*'))
                if len(fastq_glob) < 2:
                    log_write(f'[ERROR]: fewer than 2 GEX FASTQs found for the `flex_gex_fastqs` prefix `{fastq}` under {flex_fastq_path}')
                    log_write(f'  Found: {[str(f.name) for f in fastq_glob] or "nothing"}')
                    log_write('Troubleshooting:')
                    log_write(" • cellranger multi needs at least an R1 and an R2 per library; check both reads are present and gzipped (.fastq.gz)")
                    log_write(f" • Check the prefix matches the start of the real filenames: `ls {flex_fastq_path}`")
                    log_write(" • Flex GEX reads are taken from --fastqs when it is given, and from `paths.input_path` otherwise -- make sure the one you are using points at the FASTQs")
                    sys.exit(1)

            # assemble the samplesheet
            gene_expression = {
                "reference": ref_dir,
                "probe-set": FLEX_PROBE_SET,
                "create-bam": GENERATE_BAM
            }
            samples = []
            libraries = []
            for _, sample in metadata_df.iterrows():
                probe_barcodes = split_probe_barcodes(sample['Flex Probe Barcode IDs'])
                if not probe_barcodes:
                    log_write(f"[ERROR]: sample {sample['Sample Name']} uses Flex chemistry but its `Flex Probe Barcode IDs` metadata field is empty")
                    log_write('Troubleshooting:')
                    log_write(" • Fill in the `Flex Probe Barcode IDs` column for every Flex sample -- cellranger multi needs it to split the pooled libraries")
                    log_write(" • Use the probe barcode IDs from your Flex kit (e.g. BC001), separated by ',' or '|' when a sample carries more than one")
                    log_write(" • If this sample is not actually Flex, correct its `Chemistry` column instead")
                    log_write(f" • Metadata source: {METADATA_SRC}")
                    sys.exit(1)

                # emit one cellranger-multi sample per probe barcode, with sample_id
                # "<Sample Name>_<BC>" and a single probe barcode. cellranger names each
                # per_sample_outs directory after its sample_id, producing
                # count/flex/<Sample Name>_<BC>/... -- exactly the per-barcode layout the
                # downstream demux / Trekker pipeline (sc_outdir, sub-sample naming) expects. A
                # single sample_id carrying all barcodes would instead pool them into one
                # count/flex/<Sample Name>/ output that nothing downstream can find.
                sample_name = str(sample['Sample Name'])
                for bc in probe_barcodes:
                    samples.append({
                        "sample_id": f"{sample_name}_{bc}",
                        "probe_barcode_ids": bc,
                        "emptydrops_minimum_umis": emptydrops_min_umi
                    })
        
            for fastq in FLEX_GEX_FASTQS:
                libraries.append({
                    "fastq_id": fastq,
                    "fastqs": flex_fastq_path,
                    "feature_types": "Gene Expression"
                })

            # generate a Flex-specific input samplesheet. Write it to METADATA_PATH (alongside the
            # other samplesheets), which is exactly where run_count reads it back as
            # METADATA_PATH/multi_samplesheet.csv -- the FASTQ paths inside are absolute, so the
            # samplesheet's own location is independent of the input FASTQ dir.
            flex_samplesheet_path = format_multi_samplesheet(
                gene_expression,
                libraries,
                samples,
                outdir=METADATA_PATH
            )

            # return a Flex samplesheet
            return [flex_samplesheet_path]

    if metadata_df['Chemistry'].values[0] == 'Flex':
        return []

    # Format a samplesheet for cellranger mkfastq.
    #
    # A blank `RNA Index` is legitimate in exactly one case: a library sequenced entirely on a merged-in
    # run declares no primary index and carries `Merge RNA From BCL` + `Add RNA Index` instead, so its
    # one library belongs to that BCL's samplesheet further down rather than to this one. Writing a row
    # here anyway would hand mkfastq an indexless line for the primary BCL, which bcl2fastq rejects. The
    # same blank is why `_library_complete` does not require a primary read pair for such a sample, so
    # `row_declares` is shared with it: the writer and the completeness check must agree on which
    # libraries this samplesheet is expected to produce.
    #
    # A blank index with no merge column to explain it is a metadata mistake, and is worth stopping for
    # rather than filtering: the row names no index in any run folder, so there is no gene-expression
    # library to demultiplex at all. Skipping it quietly is the worse failure of the two -- the
    # completeness check has no read pair to require either, so the sample counts as done, mkfastq is
    # skipped, and cellranger fails hours later with nothing pointing back at the metadata. A --fastqs
    # run is exempt: it never demultiplexes, and reads its library names off the FASTQ filenames.
    samplesheet_rows = []
    indexless = []
    for _, row in metadata_df.iterrows():
        if not row_declares(row, 'RNA Index'):
            if FASTQ_INPUT is None and not row_declares(row, 'Merge RNA From BCL'):
                indexless.append(str(row['Sample Name']))
            # merge-only: written to the merged-in BCL's samplesheet below, not this one
            continue
        lanes = str(row['Lane']).split(',')
        for lane in lanes:
            if lane == '':
                continue
            samplesheet_rows.append(row[['Lane', 'Sample Name', 'RNA Index']])
            samplesheet_rows[len(samplesheet_rows) - 1]['Lane'] = lane

    if indexless:
        log_write("[ERROR]: these samples declare no gene-expression library to demultiplex — `RNA Index` is empty and no `Merge RNA From BCL` names another run to take the reads from:")
        for name in indexless:
            log_write(f"  • {name}")
        log_write("Troubleshooting:")
        log_write(" • Set `RNA Index` to the cellranger index the sample's RNA library carries in this run")
        log_write(" • If that library was sequenced on a different run, leave `RNA Index` empty and set `Merge RNA From BCL` to that run's BCL ID and `Add RNA Index` to the index it carries there")
        log_write(f" • Check the row belongs to this run at all — `Run` should be YES only for samples sequenced in {BCL_ID}")
        log_write(" • If these reads are already demultiplexed, pass them with --fastqs, which needs no index")
        sys.exit(1)

    if samplesheet_rows:
        samplesheet_df = pd.DataFrame(samplesheet_rows)
        samplesheet_df.columns=['Lane', 'Sample', 'Index']
    else:
        # every gene-expression library is merge-only, so the primary BCL has no RNA row. Built with
        # explicit columns because an empty DataFrame has none to rename, and mkfastq still needs the
        # header -- the spatial libraries below, or the split samplesheets, carry the actual work.
        samplesheet_df = pd.DataFrame(columns=['Lane', 'Sample', 'Index'])

    # The spatial-barcode libraries, selected with the same predicate for the same reason. A blank
    # `SB Index` is filtered rather than fatal: a merge-only spatial library legitimately has none, and
    # the spatial stages are skippable in a way the gene-expression side is not. It is still named,
    # since a typo looks identical here and would otherwise surface only when spatial-count found
    # nothing to count.
    declares_sb = [row_declares(row, 'SB Index') for _, row in metadata_df.iterrows()]
    spatial_barcodes = metadata_df.loc[declares_sb, ['Sample Name', 'SB Index', 'SB Lane']]

    no_spatial = [str(row['Sample Name']) for _, row in metadata_df.iterrows()
                  if not row_declares(row, 'SB Index') and not row_declares(row, 'Merge Spatial From BCL')]
    if no_spatial:
        log_write("[WARNING]: these samples declare no spatial-barcode library — `SB Index` and `Merge Spatial From BCL` are both empty, so none will be demultiplexed for them:")
        for name in no_spatial:
            log_write(f"  • {name}")
        log_write(" • The spatial-count and spatial-analysis stages will have nothing to run on for these samples")
        log_write(" • Set `SB Index`, or `Merge Spatial From BCL` + `Add SB Index` if the library was sequenced on another run")

    # add spatial barcode samples to the samplesheet
    for _, row in spatial_barcodes.iterrows():
        for lane in str(row['SB Lane']).split(','):
            new_row = {'Lane': lane, 'Sample': (str(row['Sample Name']) + '_sb'), 'Index': str(row['SB Index'])}
            samplesheet_df = pd.concat([samplesheet_df, pd.DataFrame([new_row])], ignore_index=True)

    # identify samples with split BCLs. merge_bcls() is the one place the metadata's merge columns
    # are read into a list of BCL IDs, so the samplesheet written per BCL here is named for exactly
    # the ID run_mkfastq later stages, demultiplexes and reads that samplesheet back for.
    split_samplesheets = {
        bcl: pd.DataFrame(columns=['Lane', 'Sample', 'Index'])
        for bcl in merge_bcls(metadata_df)
    }

    # create individual samplesheets for samples with split BCLs
    missing_add_index = []
    for _, sample in metadata_df.iterrows():
        for column, suffix, index_column in (
            ('Merge RNA From BCL', '_split_rna', 'Add RNA Index'),
            ('Merge Spatial From BCL', '_split_sb', 'Add SB Index'),
        ):
            if column not in metadata_df:
                continue
            bcl = sample[column]
            if bcl is None or pd.isna(bcl) or not str(bcl).strip():
                continue
            bcl = str(bcl).strip()
            if bcl not in split_samplesheets:
                # merge_bcls() drops a merge column that names this run's own BCL, since there is no
                # second run folder to demultiplex; its rows are already covered by the main
                # samplesheet, and writing them to a <BCL_ID>_samplesheet.csv would only be
                # overwritten by the main one below
                log_write(f"[WARNING]: sample {sample['Sample Name']} sets `{column}` to {BCL_ID}, this run's own BCL; ignoring it")
                log_write(f" • Reads from {BCL_ID} are already demultiplexed by the main samplesheet")
                log_write(f" • To merge reads from a second sequencing run, set `{column}` to that run's BCL ID")
                continue
            # A merge column names the run to demultiplex; `Add ... Index` names the index that library
            # carries there. Without the second, the row below writes a blank index -- the same line
            # bcl2fastq rejects that the primary-index check above exists to prevent, and one the
            # completeness check would then wait on forever. Checked after the own-BCL guard, since a
            # merge column pointing at this run is ignored anyway and needs no index. Collected rather
            # than raised on the spot so a single pass names every offending sample.
            if not row_declares(sample, index_column):
                missing_add_index.append((str(sample['Sample Name']), column, bcl, index_column))
                continue
            new_row = {'Lane': '*', 'Sample': (str(sample['Sample Name']) + suffix), 'Index': str(sample[index_column])}
            split_samplesheets[bcl] = pd.concat([split_samplesheets[bcl], pd.DataFrame([new_row])], ignore_index=True)

    if missing_add_index:
        log_write("[ERROR]: these samples merge reads from another run but do not say which index that library carries there:")
        for name, column, bcl, index_column in missing_add_index:
            log_write(f"  • {name}: `{column}` is {bcl}, but `{index_column}` is empty")
        log_write("Troubleshooting:")
        log_write(" • Set the `Add ... Index` column to the cellranger index the library was sequenced with in that run")
        log_write(" • It is usually a different index than the primary one; the same index in both runs produces identically named FASTQs")
        log_write(" • Clear the `Merge ... From BCL` column instead if the sample does not in fact merge reads from another run")
        sys.exit(1)

    # set up the output samplesheets list
    output_samplesheets = []

    # write samplesheets to CSV
    for bcl in split_samplesheets:
        split_samplesheets[bcl] = split_samplesheets[bcl][['Lane', 'Sample', 'Index']]
        samplesheet_name = METADATA_PATH / f"{bcl}_samplesheet.csv"
        split_samplesheets[bcl].to_csv(samplesheet_name, index=False)
        output_samplesheets.append(samplesheet_name)

    # save the main samplesheet to CSV
    global SAMPLESHEET_PATH
    if split_samplesheets:
        main_samplesheet = METADATA_PATH / f"{BCL_ID}_samplesheet.csv"
        samplesheet_df.to_csv(main_samplesheet, index=False)

        # insert the main samplesheet at the beginning of the list
        output_samplesheets.insert(0, main_samplesheet)
    else:
        main_samplesheet = METADATA_PATH / f"samplesheet.csv"
        samplesheet_df.to_csv(main_samplesheet, index=False)
        output_samplesheets.append(main_samplesheet)

    SAMPLESHEET_PATH = main_samplesheet

    return output_samplesheets


def stage_input_data(bcl_ids: list[str]) -> None:
    """
    Download this run's input folders from the `gs://` `paths.input_path` prefix into the local
    directory `settings.gcs_download_dest` names (config.py re-points `input_path` at it when staging).

    Every folder of reads a staged run needs comes through here: the primary `<BCL_ID>`, the extra run
    folders a split-BCL run merges from (`Merge RNA From BCL` / `Merge Spatial From BCL`), and the
    already-demultiplexed FASTQ folder of a `--fastqs` run. The launcher scripts deliberately no longer
    copy any of it -- `slidr_gcp.sh`/`slidr_slurm.sh` know only the primary BCL ID, so leaving the whole
    job here is what lets a split-BCL run work at all, and keeps one implementation with one
    already-staged rule and one destination.

    A no-op unless this run stages from GCS: without --stage-gcs/--gcp `input_path` is a local
    directory of run folders, and a missing one is reported by the caller's own check rather than
    silently downloaded.

    Inputs:
     - bcl_ids: BCL run IDs to stage
    Output:
     - none; each folder is left at resolve_bcl_dir(<bcl_id>)
    """

    if not STAGE_GCS or not bcl_ids:
        return

    for bcl in bcl_ids:
        dest = resolve_bcl_dir(bcl)

        # Sequencing data is the expensive part of staging (routinely hundreds of GB per run), so a
        # folder an earlier run staged into this same destination is left alone. RESTAGE=1 forces a
        # re-download, which is also the way out of a download that was interrupted so early that the
        # folder looks populated.
        if dest.is_dir() and any(dest.iterdir()) and os.environ.get('RESTAGE') != '1':
            log_write(f"  Input data for {bcl} already staged at {dest} (set RESTAGE=1 to re-download)")
            continue

        # `gcloud storage cp -r SRC DEST` requires DEST to be an existing directory and nests SRC
        # under it as DEST/<leaf>. The source leaf here *is* <bcl>, so copying into dest.parent lands
        # the run folder exactly at dest -- the same shape the reference genome and raw-barcode call
        # sites rely on. Passing `dest` itself is what will not work: absent it is rejected outright
        # ("Destination URL must name an existing directory"), and present it nests the run folder at
        # <dest>/<bcl>/ where nothing looks for it.
        if dest.exists():
            shutil.rmtree(dest)  # clear out a partial earlier download
        dest.parent.mkdir(parents=True, exist_ok=True)

        log_write(f"  Staging input data for {bcl} from {INPUT_BUCKET}... ", terminator="")
        with console.status(f"  Staging input data for {bcl} (this can take a while)..."):
            stage_from_gcs(
                f'{INPUT_BUCKET}/{bcl}',
                dest.parent,
                recursive=True,
                description=f'the {bcl} input data'
            )
        log_write("Done.")


def gcs_dest_taken(uri: str) -> bool:
    """
    True if `uri` already names something in the bucket -- an object or a prefix with contents.

    Used to keep an upload from landing on top of an earlier run's results. `gcloud storage ls` exits
    non-zero with "matched no objects" for a free name, which is the signal wanted here; any *other*
    failure (no such bucket, no permission, no credentials) says nothing about whether the name is
    free, so it is reported and treated as free -- the upload that follows fails for the same reason
    moments later, with gcloud's own diagnosis, rather than being silently renamed around a problem
    that has nothing to do with collisions.

    Inputs:
     - uri: gs:// URI to test
    Output:
     - True if the name is already in use
    """

    cmd = ['gcloud', 'storage', 'ls', uri]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=GCS_METADATA_TIMEOUT)
    except FileNotFoundError:
        log_write("[WARNING]: `gcloud` was not found on PATH, so the output bucket could not be checked for existing results")
        return False
    except subprocess.TimeoutExpired:
        # Non-fatal, like every other failure here: the results exist on local disk either way, and
        # refusing to upload them because a name check timed out would be the worse outcome. Reported
        # as "not in use", so the caller uploads to `uri` -- which is why the overwrite risk is named.
        log_write(f"[WARNING]: checking whether {uri} already exists timed out after {GCS_METADATA_TIMEOUT}s")
        log_write("  Uploading to that location anyway; anything already there may be overwritten")
        return False

    if proc.returncode == 0:
        return bool((proc.stdout or '').strip())

    # gcloud's wording for "nothing there", which is the one non-zero exit that is not a problem
    stderr = (proc.stderr or '').lower()
    if 'matched no objects' in stderr or 'not found' in stderr:
        return False

    detail = (proc.stderr or proc.stdout or '').strip().splitlines()
    log_write(f"[WARNING]: could not check whether {uri} already exists (`{' '.join(cmd)}` exited {proc.returncode})")
    for line in detail[-3:] or ['(no output)']:
        log_write(f"  {line}")
    log_write("  Uploading to that location anyway; anything already there may be overwritten")
    return False


def unique_gcs_dest(base: str, attempts: int = 100) -> str:
    """
    First unused variant of `base`: `base` itself when the bucket has nothing there, else `base_2`,
    `base_3`, and so on.

    This is what keeps a re-run of a BCL from overwriting the results of the earlier one. Numbering
    rather than replacing is deliberate: an upload is the last thing a run does, hours of compute after
    the point where a mistake could still be corrected, and the pipeline cannot tell a deliberate
    re-run from an accidental second launch of the same BCL. Keeping both and letting the operator
    delete one is recoverable; overwriting is not.

    Two things it is honestly not: an exclusive claim -- two runs finishing at the same moment can
    resolve to the same free name, since GCS offers no way to reserve a prefix -- and a merge, so a
    re-run of a single stage uploads to a new folder holding only what that run produced rather than
    updating the complete set next to it.

    Inputs:
     - base:     gs:// URI the results would go to if nothing were in the way
     - attempts: how many numbered variants to try before falling back to a timestamp
    Output:
     - a gs:// URI that was unused when checked
    """

    if not gcs_dest_taken(base):
        return base

    for suffix in range(2, attempts + 1):
        candidate = f'{base}_{suffix}'
        if not gcs_dest_taken(candidate):
            return candidate

    # Reaching here means 100 numbered folders are already in the bucket, which is far likelier to be
    # a misconfigured output_bucket than a hundredth re-run -- but the results still have to go
    # somewhere, and a timestamp needs no further probing to be distinct.
    return f"{base}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


def upload_outputs() -> None:
    """
    Copy this run's `output/` tree to `settings.output_bucket` once the run has succeeded.

    The destination is `<output_bucket>/<BCL_ID>`, so results are identifiable in a bucket that
    collects more than one run, and `unique_gcs_dest` moves it to `<BCL_ID>_2` (then `_3`, ...) rather
    than writing over a folder that is already there.

    A `--gcp` run is the one case where the folder is not chosen here. `./slidr` has to resolve it
    before creating the VM, since that is where it uploads the config.yaml the VM boots from, and it
    passes the answer down as `SLIDR_OUTPUT_DEST`. Re-deciding here would ignore that and put the
    results somewhere other than beside their own config — and would in fact always pick the *next*
    name, because the config sitting there makes the folder look taken.

    The copy names the contents of `output/` as its sources rather than `output/` itself, with a
    trailing slash on the destination. Both of those tell gcloud the destination is a container to place
    the sources in; `cp -r output <dest>` would instead depend on whether `<dest>` happens to exist
    already — writing output/'s contents as `<dest>` when it does not, and nesting them at
    `<dest>/output/` when it does. Under `SLIDR_OUTPUT_DEST` it always does exist, since the launcher
    just put a config.yaml there, so leaving that to chance would put the results one level deeper than
    a locally-resolved run's.

    Only `output/` is uploaded, as before: `data/` holds staged sequencing reads that came out of a
    bucket in the first place, and `log/`, `metadata/` and `tmp/` stay on the machine that ran the job.

    Within `output/`, the directories a staged run downloads its *inputs* into are skipped as well --
    see UPLOAD_SKIP. They sit there only because a staged run needs somewhere local to put them, and
    sending them back is pure cost: on a mouse run that is the 14 GiB reference genome and 2.5 GiB of
    Cellranger, re-uploaded into the results folder they were copied out of moments earlier.

    Output:
     - none; a failure is a warning, since the outputs themselves are already on local disk
    """

    if OUTPUT_DEST is None and (OUTPUT_BUCKET is None or not str(OUTPUT_BUCKET).strip()):
        return

    entries = sorted(
        entry for entry in (OUTPUT_PATH.iterdir() if OUTPUT_PATH.is_dir() else [])
        if entry.name not in UPLOAD_SKIP
    )
    if not entries:
        log_write(f"[WARNING]: nothing to upload to {OUTPUT_DEST or OUTPUT_BUCKET} — {run_relative(OUTPUT_PATH)} is empty")
        return

    bucket = gcs_uri(OUTPUT_DEST or OUTPUT_BUCKET)
    if OUTPUT_DEST is not None:
        dest, base = bucket, bucket
    else:
        base = f'{bucket}/{BCL_ID}'
        dest = unique_gcs_dest(base)
    if dest != base:
        log_write(f"  {base} already holds results, so this run is uploaded alongside them rather than over them")

    log_write(f'  Uploading the results to GCP: {dest}... ', terminator='')

    # sources are output/'s contents and the destination carries a trailing slash: see the docstring
    cmd = ['gcloud', 'storage', 'cp', '-r'] + [str(entry) for entry in entries] + [f'{dest}/']
    try:
        # run_gcloud_transfer, not subprocess.run(capture_output=True): this exact call is the one
        # that deadlocked for 3h45m when a gcloud worker outlived its parent still holding the
        # stdout pipe, and it is the worst place in the pipeline to hang -- the analysis is already
        # finished by the time it runs. See run_gcloud_transfer's docstring for the mechanism.
        with console.status("  Uploading the results (this can take a while)..."):
            proc = run_gcloud_transfer(cmd)
        if proc.returncode != 0:
            # surface gcloud's own last words rather than a bare exit code, which is undiagnosable
            detail = (proc.stderr or proc.stdout or '').strip().splitlines()
            raise RuntimeError(f"`gcloud storage cp` exited {proc.returncode}: {detail[-1] if detail else '(no output)'}")
        log_write(' Done.')
        log_write(f"  Outputs uploaded to:  {dest}", SUMMARY_LOG, terminal=False)
    except Exception as exc:
        bucket_name = bucket.replace('gs://', '').split('/')[0]
        log_write(f'\n[WARNING]: Could not upload outputs to Google Cloud: {exc}')
        log_write("Troubleshooting:")
        log_write(" • Make sure gcloud is installed: `gcloud --version`")
        log_write(" • Make sure you are authenticated: `gcloud auth print-access-token`")
        log_write(" • If not, authenticate with `gcloud auth login`")
        log_write(f" • Make sure the output bucket {bucket_name} exists: `gcloud storage buckets describe {bucket_name}`")
        if isinstance(exc, subprocess.TimeoutExpired):
            log_write(f" • The transfer was killed after {GCS_TRANSFER_TIMEOUT}s; raise the limit with SLIDR_GCS_TRANSFER_TIMEOUT=<seconds>")
            log_write(" • `gcloud storage rsync -r` resumes a partial upload, where `cp -r` restarts it")
        log_write(f" • The outputs are still on this machine at {OUTPUT_PATH}; run the upload manually:")
        log_write(f"    `{' '.join(cmd)}`")


def upload_diagnostics() -> None:
    """
    Copy this run's `log/` and `metadata/` to the results folder, however the run ended.

    `upload_outputs` uploads `output/`, and only once the run has succeeded, which leaves a `--gcp`
    run with no record of itself in either outcome: the VM self-deletes when the startup script exits,
    taking runtime.log, the run summary and every per-stage tool log (mkfastq.log, count.log,
    cellbender.log, ...) with it. A failure is then undiagnosable in precisely the case where the
    machine cannot be inspected afterwards -- `job_failure` names `log/<stage>.log` as the place the
    real cause is written, and by the time anyone reads that sentence the file is gone with the disk --
    and a success loses its own provenance.

    Registered with `atexit` by main.py rather than called at the end of it, so that it also runs on
    the paths that leave through `sys.exit` instead of reaching the bottom of the script: every tool
    crash (`helpers.job_failure`) and every validation error. It cannot cover a hard kill -- SIGKILL
    from the OOM killer, a preempted VM -- because nothing runs then.

    Deliberately not `output/`: that tree holds the FASTQs mkfastq wrote, routinely tens of GB, and
    pushing it to a bucket on every crash would cost far more than it saves. Preserving a failed run's
    partial results is a separate problem, and one that would need the pipeline to fetch them back at
    the start of the next attempt to be worth anything.

    A no-op unless `SLIDR_OUTPUT_DEST` names the folder to upload to, which is the `--gcp` case: local
    and Slurm runs keep their logs on a filesystem that outlives the run, and re-resolving a
    destination here could pick a different `<BCL_ID>_N` folder than `upload_outputs` settled on,
    scattering one run's artifacts across two.

    Output:
     - none; every failure is swallowed, because this runs while the process is already on its way out
       and must not replace the error the user actually needs to read
    """

    try:
        if OUTPUT_DEST is None:
            return

        entries = [
            Path(p) for p in (LOG_PATH, METADATA_PATH)
            if Path(p).is_dir() and any(Path(p).iterdir())
        ]
        if not entries:
            return

        dest = gcs_uri(OUTPUT_DEST)
        cmd = ['gcloud', 'storage', 'cp', '-r'] + [str(entry) for entry in entries] + [f'{dest}/']

        # same reasoning as upload_outputs: a bulk `gcloud storage cp` must not be read through a
        # pipe, or an orphaned worker holding the write end can block the call forever. Here that
        # would hang a process that is already exiting, which is worse -- the VM's self-delete trap
        # never fires and the instance bills until someone notices.
        proc = run_gcloud_transfer(cmd, timeout=GCS_METADATA_TIMEOUT)
        if proc.returncode == 0:
            log_write(f"  Logs uploaded to: {dest}")
        else:
            detail = (proc.stderr or proc.stdout or '').strip().splitlines()
            log_write(f"[WARNING]: could not upload this run's logs to {dest}: {detail[-1] if detail else '(no output)'}")
    except Exception:
        # bookkeeping must never mask the failure that brought us here
        pass


def run_mkfastq() -> None:
    """
    Run cellranger mkfastq on samples specified in the samplesheet to generate FASTQ files from BCLs in DATA_PATH
    """

    # read metadata file and set up the logfile and output directory
    global proc
    metadata_df = pd.read_csv(SUMMARY_PATH)
    mkfastq_log = LOG_PATH / "mkfastq.log"
    MKFASTQ_OUTS.mkdir(exist_ok=True)

    # Retrieve executable paths. mkfastq demultiplexes every sample in one invocation, so unlike
    # count it cannot honour a per-sample `Cellranger` override; take the single declared version, and
    # say so plainly if the samples disagree rather than silently picking one.
    declared_versions = sorted({
        resolve_cellranger_version(v)
        for v in (metadata_df['Cellranger'] if 'Cellranger' in metadata_df.columns else [None])
    })
    mkfastq_version = declared_versions[0]
    if len(declared_versions) > 1:
        log_write(f"[WARNING]: the `Cellranger` metadata column names more than one version "
                  f"({', '.join(declared_versions)}); demultiplexing the whole run with {mkfastq_version}")
        log_write(" • mkfastq processes every sample in a single invocation, so only one version can be used")
        log_write(" • The per-sample versions are still honoured by the count stage")

    log_detail(f"cellranger {mkfastq_version} for demultiplexing", terminal=False)
    cellranger_path = test_and_install_software('cellranger', version=mkfastq_version)
    bcl2fastq_path = test_and_install_software('bcl2fastq')

    # enable cellranger to call bcl2fastq by modifying PATH
    env = os.environ.copy()
    env["PATH"] = str(Path(bcl2fastq_path).parent) + ':' + os.environ["PATH"]

    stage_started = log_stage_start("cellranger mkfastq")

    # every BCL ID this run demultiplexes: the primary one, plus the run folders a split-BCL sample
    # merges extra RNA/spatial reads from
    extra_bcls = merge_bcls(metadata_df)
    bcls = [BCL_ID] + extra_bcls
    tmp_dirs = {}

    # Bring the run folders down when this run stages from GCS. Which ones are needed depends on the
    # metadata, so this is the earliest point it is knowable -- and doing it here rather than at
    # startup means a re-run whose FASTQs already exist never spends the download at all.
    stage_input_data(bcls)

    # Refuse to hand cellranger a run folder that is not one. Checked after staging, since a staged
    # folder does not exist until then, and for the primary BCL as well as the merged-in ones: this is
    # the only place all of them are known together.
    for bcl in bcls:
        bcl_dir = resolve_bcl_dir(bcl)
        if not validate_bcl_dir(bcl_dir):
            log_write(f"[ERROR]: {bcl_dir} is not a usable Illumina BCL run directory (see the warnings above for what is missing)")
            log_write("Troubleshooting:")
            if bcl == BCL_ID:
                log_write(f" • Check that `input_path` ({INPUT_PATH}) holds the {BCL_ID} run folder")
                log_write(" • Check the BCL run finished transferring (RunInfo.xml, RunParameters.xml and Data/Intensities/BaseCalls must all be present)")
                log_write(" • If these reads are already demultiplexed, pass them with --fastqs instead -- the pipeline will not guess this for you")
            else:
                log_write(f" • {bcl} is named by a `Merge RNA From BCL` / `Merge Spatial From BCL` column, so its run folder has to be demultiplexed alongside {BCL_ID}")
                log_write(f" • Check that run folder sits beside the {BCL_ID} one, under `input_path` ({INPUT_PATH})")
                log_write(" • Clear the column for that sample if its reads should not be merged after all")
            if STAGE_GCS:
                log_write(f" • Check it exists in the bucket: `gcloud storage ls {INPUT_BUCKET}/{bcl}`")
                log_write(" • Re-download a partially staged folder by re-running with RESTAGE=1 in the environment")
            sys.exit(1)

    # create output directories
    for _, sample in metadata_df.iterrows():
        sample_path = MKFASTQ_OUTS / sample['Sample Name']
        sample_path.mkdir(parents=True, exist_ok=True)

        spatial_barcode_path = MKFASTQ_OUTS / f"{sample['Sample Name']}_sb"
        spatial_barcode_path.mkdir(parents=True, exist_ok=True)

    # build an exact library-name -> destination-directory map for routing demultiplexed FASTQs.
    # Each samplesheet library (primary GEX <sample>, primary spatial <sample>_sb, and the merged
    # <sample>_split_rna / <sample>_split_sb libraries created by create_samplesheet) maps to the
    # folder its reads belong in. Matching the cellranger sample name exactly -- rather than a
    # substring of the file path -- avoids misrouting when one sample name is a prefix of another
    # (e.g. "tumor" vs "tumor_ln").
    route_map = {}
    for _, sample in metadata_df.iterrows():
        sname = str(sample['Sample Name'])
        route_map[sname] = MKFASTQ_OUTS / sname                        # primary GEX library
        route_map[f"{sname}_split_rna"] = MKFASTQ_OUTS / sname         # merged-in GEX library
        route_map[f"{sname}_sb"] = MKFASTQ_OUTS / f"{sname}_sb"        # primary spatial library
        route_map[f"{sname}_split_sb"] = MKFASTQ_OUTS / f"{sname}_sb"  # merged-in spatial library

    # cellranger/bcl2fastq FASTQ names are "<library>_S<n>_L<lane>_<R1|R2|I1|I2>_001.fastq.gz";
    # capture the library prefix so it can be matched exactly against route_map
    fastq_name_re = re.compile(r'^(?P<lib>.+)_S\d+_L\d+_[RI][12]_\d+\.fastq\.gz$')

    # track every destination written during this run so that two source FASTQs from different
    # BCLs resolving to the same destination are caught (and refused) rather than silently
    # overwriting each other via shutil.move
    moved_dests = set()

    # run mkfastq with appropriate inputs
    for bcl in bcls:
        # define samplesheet paths
        if len(bcls) > 1:
            samplesheet_path = METADATA_PATH / f"{bcl}_samplesheet.csv"
        else:
            samplesheet_path = METADATA_PATH / 'samplesheet.csv'

        # create temporary output directory to avoid cellranger output conflicts
        tmp_dir = create_tmp_dir(bcl)
        tmp_dirs[bcl] = tmp_dir

        # specify data path
        data_path = resolve_bcl_dir(bcl)

        # launch mkfastq
        log_write(f"  Generating FASTQs from BCLs in {bcl}...", terminal=False, terminator="")
        with open(mkfastq_log, "a") as logfile:
            with console.status(f"  Generating FASTQs from BCLs in {bcl}..."):
                proc = subprocess.Popen(
                    ['time', 'stdbuf', '-oL', '-eL', cellranger_path, 'mkfastq', 
                    f'--run={data_path}',
                    f'--localcores={NUM_THREADS}',
                    f'--localmem={MEM_SIZE}',
                    f'--id={bcl}',
                    f'--output-dir={tmp_dir}',
                    f'--csv={samplesheet_path}',
                    f'--delete-undetermined'],
                    stdout=logfile,
                    stderr=subprocess.STDOUT,
                    env=env,
                    cwd=tmp_dir
                )

                proc.wait()

        # check that the job finished successfully
        if proc.returncode != 0:
            job_crash("cellranger mkfastq", proc.returncode, mkfastq_log)
        else:
            log_ts(f"demultiplexed {bcl}")

        # extract the flowcell ID for each BCL run
        run_info = data_path / 'RunInfo.xml'
        try:
            tree = ET.parse(run_info)
        except Exception as exc:
            log_write(f"[ERROR]: Could not parse {run_info}: {exc}")
            log_write("Troubleshooting:")
            log_write(" • RunInfo.xml is written by the sequencer; a parse failure usually means the run folder was copied while the run was still in progress")
            log_write(f" • Check the file is complete and well-formed: `xmllint --noout {run_info}`")
            log_write(f" • Re-copy the {bcl} run folder from the sequencer and re-run with --force")
            sys.exit(1)
        try:
            flowcell_id = tree.getroot().find('.//Flowcell').text
        except Exception as exc:
            log_write(f"[ERROR]: could not extract flowcell ID from {run_info}: {exc}")
            log_write("Troubleshooting:")
            log_write(f" • {run_info} parsed but has no <Flowcell> element; check it with `grep -i flowcell {run_info}`")
            log_write(" • This is how the pipeline locates cellranger's output directory, so the run folder cannot be used without it")
            log_write(f" • Re-copy the {bcl} run folder from the sequencer and re-run with --force")
            sys.exit(1)
        flowcell_dir = tmp_dirs[bcl] / flowcell_id
        if not flowcell_dir.is_dir():
            log_write(f"[ERROR]: cellranger mkfastq reported success but wrote no output directory for flowcell {flowcell_id}: {flowcell_dir}")
            log_write("Troubleshooting:")
            log_write(f" • Check the tail of {mkfastq_log} -- mkfastq can exit 0 after demultiplexing zero reads")
            log_write(f" • Check the flowcell ID in {run_info} matches the one in the sequencing data (a mismatched run folder demultiplexes into a differently named directory)")
            log_write(f" • Check the samplesheet indices for {bcl} match the run: `cat {samplesheet_path}`")
            sys.exit(1)

        # move FASTQ files from cellranger's flowcell output into per-sample/library folders,
        # routing by the exact parsed library name
        for fastq_file in flowcell_dir.rglob("*.fastq.gz"):
            match = fastq_name_re.match(fastq_file.name)
            if match is None:
                # not a standard per-sample FASTQ (e.g. Undetermined_*); leave it in place
                continue

            dest_dir = route_map.get(match.group('lib'))
            if dest_dir is None:
                # a library that isn't in this run's metadata; don't guess a destination
                log_write(f"  [WARNING]: FASTQ {fastq_file.name} from {bcl} does not match any sample library; leaving it in place", terminal=False)
                continue

            dest = dest_dir / fastq_file.name
            if dest in moved_dests:
                # two source FASTQs (from different BCLs) collapse to the same destination file;
                # moving the second would silently destroy the first
                log_write(f"[ERROR]: FASTQ filename collision: '{fastq_file.name}' was produced for the same library by more than one BCL and both resolve to {dest}.")
                log_write("Distinct BCLs must not yield identically named FASTQs for the same sample library; moving the second would silently destroy the first.")
                log_write("Troubleshooting:")
                log_write(" • Check the `Merge RNA From BCL` / `Merge Spatial From BCL` columns -- a BCL listed there must not be the same as the `BCL` column value")
                log_write(" • Check `Add RNA Index` / `Add SB Index` differ from the primary `RNA Index` / `SB Index`; the same index in both BCLs produces identically named FASTQs")
                log_write(f" • The colliding libraries came from these BCLs: {', '.join(str(b) for b in bcls)}")
                log_write(f" • Metadata source: {METADATA_SRC}")
                sys.exit(1)

            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(fastq_file), str(dest))
            moved_dests.add(dest)

    # logging the job results
    job_success("cellranger mkfastq", mkfastq_log, MKFASTQ_OUTS, started=stage_started)


def run_cellranger_count(
    cellranger_path: str | Path,
    sample_index: str,
    sample_name: str | list,
    fastq_dir: str | Path,
    output_path: str | Path,
    ref_genome: str,
    chemistry: str,
    threads: str,
    memory: str,
    generate_bam: str,
    logpath
) -> subprocess.Popen:
    """
    Run cellranger count on a single sample and block until completion

    Inputs:
     - cellranger_path: path to the cellranger executable
     - sample_index:    sample name used as the cellranger --id and --sample arguments
     - fastq_dir:       directory containing input FASTQ files for this sample
     - output_path:     directory where cellranger count outputs will be written
     - ref_genome:      path to the cellranger-compatible reference transcriptome
     - chemistry:       sequencing chemistry string passed to cellranger (e.g. SC3Pv3)
     - threads:         number of local cores to allocate
     - memory:          local memory limit in GB
     - generate_bam:    whether to generate a BAM file ("true" or "false")
     - logpath:         path to the file where stdout/stderr will be written
    Output:
     - a completed subprocess.Popen object; check .returncode for exit status
    """

    # define the global process variable and create output directory
    global proc
    output_path = Path(output_path)
    output_path.mkdir(exist_ok=True, parents=True)

    # if processing multiple samples, concatenate them with commas
    if isinstance(sample_name, list):
        sample_name = ','.join(sample_name)

    # launch cellranger count as a child process
    with open(logpath, "a") as logfile:
        proc = subprocess.Popen(
            ['time', 'stdbuf', '-oL', '-eL', cellranger_path, 'count',
            f'--id={sample_index}',
            f'--sample={sample_name}',
            f'--fastqs={fastq_dir}',
            f'--transcriptome={ref_genome}',
            f'--chemistry={chemistry}',
            f'--localcores={threads}',
            f'--localmem={memory}',
            f'--jobmode=local',
            f'--project=fastqs',
            f'--disable-ui',
            f'--nosecondary',
            f'--include-introns=true',
            f'--create-bam={generate_bam}'],
            cwd=output_path.parent,
            stdout=logfile,
            stderr=subprocess.STDOUT
        )

        proc.wait()

    return proc


def run_cellranger_multi(
    cellranger_path: str | Path,
    sample_index: str,
    samplesheet: str | Path,
    output_path: str | Path,
    threads: str,
    memory: str,
    logpath: str | Path
) -> subprocess.Popen:
    """
    Run cellranger multi on a single sample and block until completion

    Inputs:
     - cellranger_path: path to the cellranger executable
     - sample_index:    sample name used as the cellranger --id and --sample arguments
     - samplesheet:     path to formatted cellranger multi samplesheet (CSV)
     - output_path:     directory where cellranger count outputs will be written
     - threads:         number of local cores to allocate
     - memory:          local memory limit in GB
     - logpath:         path to the file where stdout/stderr will be written
    Output:
     - a completed subprocess.Popen object; check .returncode for exit status
    """

    # define the global process variable and create output directory
    global proc
    output_path = Path(output_path)
    output_path.mkdir(exist_ok=True, parents=True)

    # launch cellranger count as a child process
    with open(logpath, "a") as logfile:
        proc = subprocess.Popen(
            ['time', 'stdbuf', '-oL', '-eL', cellranger_path, 'multi',
            f'--id={sample_index}',
            f'--csv={samplesheet}',
            f'--localcores={threads}',
            f'--localmem={memory}',
            f'--jobmode=local',
            f'--disable-ui'],
            cwd=output_path,
            stdout=logfile,
            stderr=subprocess.STDOUT
        )

        proc.wait()

    return proc


def run_count() -> None:
    """
    Run cellranger count on samples specified in the samplesheet to generate readcount files
    """

    # read metadata file, set up output directory, and extract the species
    count_log = LOG_PATH / "count.log"
    COUNT_OUTS.mkdir(exist_ok=True)
    metadata_df = pd.read_csv(SUMMARY_PATH)
    is_flex = 'Flex' in metadata_df['Chemistry'].values

    # set up run parameters
    generate_bam = 'true' if GENERATE_BAM else 'false'

    # verify reference genome folder exists
    genome = staged_ref_dir()
    if not genome.is_dir():
        log_write(f'[ERROR]: reference genome directory {genome} not found')
        log_write("Troubleshooting:")
        log_write(f" • Check that `settings.reference_genome` ({REF_GENOME}) names a directory that exists under `paths.reference_path` ({REF_PATH})")
        if STAGE_GCS:
            log_write(f" • Check the reference exists in the bucket: `gcloud storage ls {gcs_uri(REF_PATH)}/{REF_GENOME}`")
            log_write(f" • Delete the partial local copy at {genome} and re-run to stage it again")
        else:
            log_write(" • If the reference lives in a bucket rather than on this filesystem, pass --stage-gcs to download it")
            log_write(" • Cellranger references are downloaded from https://www.10xgenomics.com/support/software/cell-ranger/downloads")
        sys.exit(1)

    # The FASTQ-input validation below belongs to the standard (non-Flex) count path only. Flex runs
    # from a cellranger-multi samplesheet built (and whose GEX FASTQs are validated) in
    # create_samplesheet, and takes its FASTQ input directory from that samplesheet.
    if not is_flex:
        # Where the FASTQs come from is decided by --fastqs alone, not by probing the input.
        if FASTQ_INPUT is not None:
            # The directory itself was already validated in config.py; check here that it actually
            # contains every library this run needs, which requires the metadata and so can only
            # happen once the run's samples are known. The spatial-barcode library is required only
            # when spatial counting is part of this invocation -- a count-only run must not be
            # blocked by FASTQs it never reads -- but when it is, checking now fails fast instead of
            # after hours of cellranger count.
            fastq_path = FASTQ_INPUT
            required_modes = ('GEX', 'SB') if (args.run_all or args.spatial_count) else ('GEX',)
            missing_samples = missing_fastqs(fastq_path, metadata_df, required_modes)
            if missing_samples:
                log_write(f"[ERROR]: the following FASTQ files are missing from {fastq_path}:")
                for sample, sample_modes in missing_samples.items():
                    log_write(f" • {sample}:")
                    for mode, reads in sample_modes.items():
                        log_write(f"  - {mode}: {', '.join(reads)}")
                log_write("Troubleshooting:")
                log_write(" • FASTQ filenames must contain the sample name plus a `GEX` (gene expression) or "
                          "`SB` (spatial barcode) token, e.g. mysample_GEX_S1_L001_R1_001.fastq.gz")
                log_write(" • Check that the lanes in the `Lane`/`SB Lane` metadata columns match the lanes in the filenames")
                log_write(" • Omit --fastqs and pass --mkfastq/--run-all to demultiplex the BCLs instead")
                sys.exit(1)
        else:
            # no --fastqs: the FASTQs are the ones mkfastq wrote under output/mkfastq/
            fastq_path = MKFASTQ_OUTS
            if not fastq_path.is_dir() or not any(fastq_path.rglob('*.fastq.gz')):
                log_write(f"[ERROR]: no demultiplexed FASTQs found at {fastq_path}")
                log_write("Troubleshooting:")
                log_write(" • Run the demultiplexing stage first with --mkfastq (or the whole pipeline with --run-all)")
                log_write(" • If the reads were demultiplexed elsewhere, pass that directory with --fastqs")
                sys.exit(1)

    # ensure count logfile exists
    count_log.touch(exist_ok=True)

    # print launch message and retrieve path to cellranger executable
    stage_started = log_stage_start("cellranger count")
    log_detail(f"reference → {genome}")

    # run cellranger on single-cell Flex data
    if 'Flex' in metadata_df['Chemistry'].values:
        # create temporary output directory, read samplesheet, and print out subsamples
        # Flex goes through `cellranger multi`, which needs a v9 release, so that -- not the standard
        # 8.0.1 -- is the default here. The `Cellranger` metadata column still overrides it; a Flex
        # run is a single multi invocation, so the value is taken from the first sample rather than
        # per-sample.
        flex_version = resolve_cellranger_version(
            metadata_df['Cellranger'].iloc[0] if 'Cellranger' in metadata_df.columns else None,
            default="9.0.1"
        )
        log_detail(f"cellranger {flex_version} (Flex / cellranger multi)", terminal=False)
        cellranger_path = test_and_install_software('cellranger', version=flex_version)
        tmp_dir = create_tmp_dir('flex')
        flex_samplesheet = METADATA_PATH / "multi_samplesheet.csv"
        log_write(f"  Processing partitioned samples (10x Flex):", terminal=False)
        for sample_name in metadata_df['Sample Name'].astype(str).tolist():
            log_write(f"   - {sample_name}\n", terminal=False)

        # launch cellranger multi as a child process
        with console.status(f"  Processing partitioned samples (10x Flex)..."):
            proc = run_cellranger_multi(
                cellranger_path,
                BCL_ID,
                flex_samplesheet,
                tmp_dir,
                NUM_THREADS,
                MEM_SIZE,
                count_log
            )

        # check exit status
        if proc.returncode != 0:
            job_crash("cellranger count", proc.returncode, count_log)

        # move output files to the expect directory
        src_dir = tmp_dir / BCL_ID / 'outs' / 'per_sample_outs'
        dst_dir = COUNT_OUTS / 'flex'
        dst_dir.mkdir(exist_ok=True, parents=True)
        if not src_dir.is_dir():
            log_write(f"\n[ERROR]: cellranger multi exited successfully but its expected per-sample output directory does not exist: {src_dir}")
            log_write("Troubleshooting:")
            log_write(f" • Check the tail of {count_log} for what cellranger multi actually did")
            log_write(f" • Check the multi samplesheet is correct: `cat {flex_samplesheet}`")
            log_write(" • A samplesheet whose `[samples]` section has no probe barcodes produces no per_sample_outs; check the `Flex Probe Barcode IDs` metadata column")
            sys.exit(1)
        if len(list(src_dir.iterdir())) == 0:
            log_write(f"\n[ERROR]: cellranger multi produced an empty per-sample output directory: {src_dir}")
            log_write("Troubleshooting:")
            log_write(f" • Check the tail of {count_log} -- cellranger multi can finish with no per-sample outputs when no reads matched any probe barcode")
            log_write(" • Check the probe barcode IDs in the `Flex Probe Barcode IDs` metadata column match the ones actually used in the library prep")
            log_write(f" • Check `workflow.flex_probe_set` ({FLEX_PROBE_SET}) is the probe set for this kit")
            sys.exit(1)
        for item in src_dir.iterdir():
            shutil.move(item, dst_dir / item.name)

        log_ts("counted all Flex probe-barcode partitions")

    else:
        # iterate over samples and call cellranger. The executable is resolved per sample rather than
        # once up front, because the `Cellranger` metadata column is a per-sample override and
        # cellranger count runs one invocation per sample -- so two samples in the same run can
        # legitimately be counted with different releases.
        samples = metadata_df['Sample Name'].astype(str).tolist()

        # Map each sample to the cellranger --sample value(s) for its gene-expression library.
        # mkfastq names its output FASTQs after the samplesheet library, which is the sample name
        # itself, so no discovery is needed there. A directory supplied with --fastqs was named by
        # somebody else's demultiplexer, so the library prefixes have to be read back off the
        # filenames -- keyed on --fastqs being set, not on guessing which case we're in.
        true_gex_names = {}
        if FASTQ_INPUT is not None:
            for sample in samples:
                sample_row = metadata_df.loc[metadata_df['Sample Name'].astype(str) == sample].iloc[0]
                gex_fastqs = find_sample_fastqs(
                    FASTQ_INPUT,
                    sample,
                    'GEX',
                    declared_lanes(sample_row, 'Lane')
                )
                # missing_fastqs already guaranteed at least one R1/R2 pair per sample above; the
                # guard is here so a logic slip can never hand cellranger an empty --sample
                libraries = sorted({FASTQ_NAME_RE.match(f.name).group('lib') for f in gex_fastqs})
                if not libraries:
                    log_write(f"[ERROR]: no gene-expression FASTQs found for sample {sample} in {FASTQ_INPUT}")
                    log_write("Troubleshooting:")
                    log_write(f" • FASTQs passed with --fastqs must be named after the sample plus an optional `GEX` token, "
                              f"e.g. {sample}_GEX_S1_L001_R1_001.fastq.gz or {sample}_S1_L001_R1_001.fastq.gz")
                    log_write(" • Check that the lanes in the `Lane` metadata column match the `_L00N_` tokens in the filenames")
                    log_write(f" • List what is actually there: `ls {FASTQ_INPUT}`")
                    sys.exit(1)
                # keep a single library as a plain name so the multi-library log below stays truthful
                true_gex_names[sample] = libraries if len(libraries) > 1 else libraries[0]
        else:
            # FASTQs came from mkfastq, whose library names are the sample names
            for sample in samples:
                true_gex_names[sample] = sample

        for sample_name in true_gex_names:
            # create output directory
            tmp_dir = create_tmp_dir('count')
            tmp_path = tmp_dir / sample_name

            # define sequencing type
            chemistry = metadata_df['Chemistry'].values[0]
            if chemistry.startswith(('3P', '5P')):
                chemistry = 'SC' + chemistry

            # resolve this sample's Cellranger release from the `Cellranger` metadata column
            sample_row = metadata_df.loc[metadata_df['Sample Name'].astype(str) == sample_name].iloc[0]
            declared_version = sample_row['Cellranger'] if 'Cellranger' in sample_row.index else None
            cellranger_version = resolve_cellranger_version(declared_version)
            log_detail(f"cellranger {cellranger_version} for sample {sample_name}", terminal=False)
            cellranger_path = test_and_install_software('cellranger', version=cellranger_version)

            # run cellranger on single-cell RNA data
            with console.status(f"  Processing sample {sample_name} (scRNA-seq)... "):
                if isinstance(true_gex_names[sample_name], list):
                    print('   Analyzing subsamples:')
                    for name in true_gex_names[sample_name]:
                        print(f'    - {name}')

                proc = run_cellranger_count(
                    cellranger_path,
                    sample_name,
                    true_gex_names[sample_name],
                    fastq_path,
                    tmp_path,
                    genome,
                    chemistry,
                    NUM_THREADS,
                    MEM_SIZE,
                    generate_bam,
                    count_log
                )

            # check exit status
            if proc.returncode != 0:
                job_crash("cellranger count", proc.returncode, count_log)

            # move output files to the expect directory
            src_dir = tmp_path / 'outs'
            dst_dir = COUNT_OUTS / sample_name
            dst_dir.mkdir(exist_ok=True)
            if not src_dir.is_dir():
                log_write(f"[ERROR]: cellranger count exited successfully but its expected output directory does not exist: {src_dir}")
                log_write("Troubleshooting:")
                log_write(f" • Check the tail of {count_log} for what cellranger count actually did")
                log_write(f" • Check for a cellranger failure report at {tmp_path / '_log'}")
                log_write(" • A cellranger run interrupted part-way leaves no outs/ directory; re-run this stage with --force")
                sys.exit(1)
            if len(list(src_dir.iterdir())) == 0:
                log_write(f"[ERROR]: cellranger count produced an empty output directory: {src_dir}")
                log_write("Troubleshooting:")
                log_write(f" • Check the tail of {count_log} -- this normally means no reads were assigned to sample {sample_name}")
                log_write(f" • Check the FASTQs for this sample exist and are non-empty: `ls -l {fastq_path}`")
                log_write(f" • Check the `Chemistry` metadata column ({chemistry}) matches how the library was prepared")
                sys.exit(1)
            for item in src_dir.iterdir():
                shutil.move(item, dst_dir / item.name)

            log_ts(f"counted {sample_name}")
            check_barcode_validity(sample_name, dst_dir, chemistry)
    
    # logging the job results
    job_success("cellranger count", count_log, COUNT_OUTS, started=stage_started)


def run_cellbender() -> None:
    """
    Run cellbender on samples specified in the metadata sheet and store the outputs
    """

    # read metadata file and create cellbender logfile
    global proc
    metadata_df = pd.read_csv(SUMMARY_PATH)
    CELLBENDER_OUTS.mkdir(exist_ok=True)
    cellbender_log = LOG_PATH / "cellbender.log"
    CUDA_BIN_PATH = None

    stage_started = log_stage_start("cellbender")

    # retrieve the cellbender executable installed in a conda environment
    conda_lib = ensure_conda_env("r_env")
    conda_bender = conda_lib / "cellbender"
        
    # search for GPU support to speed up cellbender training
    env = conda_subprocess_env()
    try:
        pynvml.nvmlInit()
            
        # get driver version
        driver_version = pynvml.nvmlSystemGetDriverVersion()
        
        # get CUDA version (supported by the driver)
        cuda_version = pynvml.nvmlSystemGetCudaDriverVersion()
        
        # CUDA version is returned as an integer (e.g., 12010 for 12.1)
        major = cuda_version // 1000
        minor = (cuda_version % 1000) // 10
        
        log_write("  CUDA is available!")
        log_write(f"  - NVIDIA driver version: {driver_version}")
        log_write(f"  - Maximum supported CUDA: {major}.{minor}")
        
        pynvml.nvmlShutdown()
        CUDA_BIN_PATH = shutil.which('nvidia-smi')
        if not CUDA_BIN_PATH:
            log_write("  [WARNING]: NVIDIA drivers for CUDA (nvidia-smi) not found on the system. Running on CPU (this may take a while)...")
            log_write("   • The GPU was detected but `nvidia-smi` is not on PATH, so CellBender cannot be pointed at it")
            log_write("   • Install the NVIDIA driver utilities, or add the directory containing nvidia-smi to PATH before running\n")
        else:
            # export location of CUDA drivers to PATH
            env["PATH"] = f"{CUDA_BIN_PATH}:{env['PATH']}"
    except Exception as e:
        log_write(f"  [WARNING]: CUDA is not available ({e}). Running cellbender on CPU (this may take a while)...")
        log_write("   • On GCP, pass --gpu to ./slidr to attach a GPU to the VM")
        log_write("   • On Slurm, pass --gpu to ./slidr so the job requests one (--gres=gpu:1)")
        log_write("   • If the machine does have a GPU, check the driver is loaded: `nvidia-smi`\n")
        env["CUDA_VISIBLE_DEVICES"] = ""

    # run cellbender for each sample in the dataset
    for _, sample in metadata_df.iterrows():
        sample_name = str(sample['Sample Name'])
        log_write(f"  Processing sample {sample_name}... ", terminal=False, terminator="")

        # define inputs/outputs and parameters
        cellranger_html_path = COUNT_OUTS / sample_name / 'web_summary.html'
        input_path = COUNT_OUTS / sample_name / 'raw_feature_bc_matrix.h5'
        sample_path = CELLBENDER_OUTS / sample_name
        sample_path.mkdir(exist_ok=True)
        output_path = sample_path / 'cellbender_output'
        expected_cells = CELLBENDER_CELLS
        total_droplets = CELLBENDER_DROPLETS
        if total_droplets is not None and not is_int(total_droplets):
            log_write(f'[ERROR]: incorrect configuration file value for the `cellbender_total_droplets` field (should be an integer): {total_droplets!r}')
            log_write("Troubleshooting:")
            log_write(" • Set `workflow.cellbender_total_droplets` to a plain integer (no quotes, no decimal point), e.g. 25000")
            log_write(" • Set it to `null` to have the value detected automatically from the cellranger barcode-rank plot")
            hint = bool_value_hint(total_droplets, 'workflow.cellbender_total_droplets')
            if hint:
                log_write(hint)
            sys.exit(1)

        if expected_cells is not None and not is_int(expected_cells):
            log_write(f'[ERROR]: incorrect configuration file value for the `cellbender_estimated_cells` field (should be an integer): {expected_cells!r}')
            log_write("Troubleshooting:")
            log_write(" • Set `workflow.cellbender_estimated_cells` to a plain integer (no quotes, no decimal point), e.g. 5000")
            log_write(" • Set it to `null` to let CellBender estimate the cell count itself")
            hint = bool_value_hint(expected_cells, 'workflow.cellbender_estimated_cells')
            if hint:
                log_write(hint)
            sys.exit(1)

        if CELLBENDER_EPOCHS is not None and not is_int(CELLBENDER_EPOCHS):
            log_write(f'[ERROR]: incorrect configuration file value for the `cellbender_epochs` field (should be an integer): {CELLBENDER_EPOCHS!r}')
            log_write("Troubleshooting:")
            log_write(" • Set `workflow.cellbender_epochs` to a plain integer (no quotes, no decimal point), e.g. 160")
            log_write(" • Remove the field to use CellBender's default of 150")
            hint = bool_value_hint(CELLBENDER_EPOCHS, 'workflow.cellbender_epochs')
            if hint:
                log_write(hint)
            sys.exit(1)

        # no bool guard needed here (nor on spatial_downsampling below): bool subclasses int, not
        # float, so isinstance(True, float) is already False and a boolean is rejected as it stands
        if CELLBENDER_RATE is not None and not isinstance(CELLBENDER_RATE, float):
            log_write(f'[ERROR]: incorrect configuration file value for the `cellbender_learn_rate` field (should be a float): {CELLBENDER_RATE!r}')
            log_write("Troubleshooting:")
            log_write(" • Set `workflow.cellbender_learn_rate` to a decimal number, e.g. 0.5 -- note that `1` parses as an integer, so write `1.0`")
            log_write(" • Remove the field to use CellBender's default of 0.0001")
            sys.exit(1)

        # calculate the total number of droplets to include
        if total_droplets is None or expected_cells is None:
            msg = []
            if total_droplets is None:
                msg.append("total-droplets-included")
            if expected_cells is None:
                msg.append("expected-cells")
            if not cellranger_html_path.is_file():
                if len(msg) == 2:
                    msg = f"{msg[0]} and {msg[1]} Cellbender parameters"
                else:
                    msg = f"{msg[0]} Cellbender parameter"
                log_write(f"[ERROR]: {msg} not set, and {cellranger_html_path} is not present to extract these parameters automatically.")
                log_write("Troubleshooting:")
                log_write(f" • Run the count stage for {sample_name} first (--count, or the whole pipeline with --run-all); web_summary.html is one of its outputs")
                log_write(" • Or set `workflow.cellbender_total_droplets` and `workflow.cellbender_estimated_cells` explicitly in the config file to skip auto-detection")
                log_write(" • Cellbender's own guidance on choosing these values: https://cellbender.readthedocs.io/en/latest/usage/index.html")
                sys.exit(1)
            cellbender_params = parse_cellranger_html(cellranger_html_path)
            if total_droplets is None:
                total_droplets = cellbender_params['total_droplets_included']
            if expected_cells is None:
                expected_cells = cellbender_params['expected_cells']

        # begin assembling cellbender command
        cmd = [
            conda_bender,
            "remove-background"
        ]

        if CUDA_BIN_PATH:
            cmd.append("--cuda")

        msg = []
        if total_droplets is not None:
            msg.append(f"   - Total droplets included: {total_droplets}")
            cmd.extend(["--total-droplets-included", str(total_droplets)])
        if expected_cells is not None:
            msg.append(f"   - Expected cells: {expected_cells}")
            cmd.extend(["--expected-cells", str(expected_cells)])
        if CELLBENDER_RATE is not None:
            msg.append(f"   - Learning rate: {CELLBENDER_RATE}")
            cmd.extend(["--learning-rate", str(CELLBENDER_RATE)])     
        else:
            msg.append(f"   - Learning rate: 0.0001 (default)")       
        if CELLBENDER_EPOCHS is not None:
            msg.append(f"   - Epochs: {CELLBENDER_EPOCHS} \n")
            cmd.extend(["--epochs", str(CELLBENDER_EPOCHS)])     
        else:
            msg.append(f"   - Epochs: 150 (default) \n")   

        with open(SUMMARY_LOG, "a") as summary:
            summary.write(f"  Cellbender parameters:\n")
            summary.write(f"   • Droplets included: {total_droplets}\n")
            summary.write(f"   • Expected cells:    {expected_cells}\n")
            summary.write(f"   • Learning rate:     {CELLBENDER_RATE}\n")
            summary.write(f"   • Training epochs:   {CELLBENDER_EPOCHS}\n")

        if len(msg) > 0:
            log_write(msg)
            
        cmd.extend(["--input", input_path])
        cmd.extend(["--output", output_path])

        # run cellbender with a modified environment
        with open(cellbender_log, "a") as logfile:
            with console.status(f"  Processing sample {sample_name}... "):
                proc = subprocess.Popen(
                    cmd,
                    stdout=logfile,
                    stderr=subprocess.STDOUT,
                    env=env,
                    cwd=sample_path
                )

                proc.wait()
        
        # check exit status
        if proc.returncode == 0:
            log_ts(f"denoised {sample_name}")
        else:
            job_crash("cellbender", proc.returncode, cellbender_log)
    
    # log job execution result
    job_success("cellbender", cellbender_log, CELLBENDER_OUTS, started=stage_started)


def run_spatial_positioning() -> None:
    """
    For each sample, generate a puck CSV if not already present and run spatial_count.jl
    to produce a SBcounts.h5 file mapping cell barcodes to spatial coordinates
    """

    # retrieve path to Julia interpreter and load metadata CSV and spatial barcode scripts
    global proc
    julia_interp_path = test_and_install_software("julia")
    script_base_path = SCRIPT_PATH / "spatial_barcodes"
    spatial_barcodes_log = LOG_PATH / "spatial_barcodes.log"
    SPATIAL_COUNT_OUTS.mkdir(exist_ok=True)
    metadata_df = pd.read_csv(SUMMARY_PATH)

    stage_started = log_stage_start("spatial barcode counting")

    # verify that environment files are present and create required subdirectories
    main_env_path = ROOT_PATH / "envs"
    if not (main_env_path / "julia" / "Manifest.toml").is_file():
        log_write(f"[ERROR]: envs/julia/Manifest.toml is not present in the project folder ({main_env_path / 'julia'})")
        log_write("Troubleshooting:")
        log_write(" • This file is checked into the repository; restore it with `git checkout envs/julia/Manifest.toml`")
        log_write(" • Or run `git pull` to restore the whole envs/julia directory")
        sys.exit(1)

    if not (main_env_path / "julia" / "Project.toml").is_file():
        log_write(f"[ERROR]: envs/julia/Project.toml is not present in the project folder ({main_env_path / 'julia'})")
        log_write("Troubleshooting:")
        log_write(" • This file is checked into the repository; restore it with `git checkout envs/julia/Project.toml`")
        log_write(" • Or run `git pull` to restore the whole envs/julia directory")
        sys.exit(1)

    (main_env_path / 'julia' / 'packages').mkdir(exist_ok=True)
    (main_env_path / 'julia' / 'tmp').mkdir(exist_ok=True)


    # export environment variables
    env = os.environ.copy()
    env["JULIA_PACKAGES_PATH"] = str(main_env_path / 'julia' / 'packages')
    env["JULIA_PROJECT_PATH"] = str(main_env_path / 'julia')
    env["JULIA_DEPOT_PATH"] = str(main_env_path / 'julia' / 'tmp')
    env["JULIA_NUM_THREADS"] = str(NUM_THREADS)

    # puck_path is required by this stage but validated as optional at config load, because only here
    # is it known that the run actually needs it (a Flex run never reaches this function). Report it
    # now, before any puck is looked up.
    if PUCK_PATH is None:
        log_write("[ERROR]: `puck_path` is not set in the configuration file, but the spatial barcode counting stage needs it")
        log_write("Troubleshooting:")
        if STAGE_GCS:
            log_write(" • Set `paths.puck_path` to the GCS location holding the puck CSVs (e.g. gs://my-bucket/pucks)")
        else:
            log_write(" • Set `paths.puck_path` to the local directory holding the puck CSVs")
        log_write(" • The puck CSV for each sample's `Puck ID` is read from there, or generated into it from `raw_barcodes_path`")
        log_write(" • Skip this stage instead by omitting --spatial-count/--run-all")
        sys.exit(1)

    # generate_puck_csv.jl reads raw barcode directories from LIB_PUCK_IN and *writes* the puck CSVs
    # it builds to LIB_PUCK_PATH, so when staging from GCS both must point at writable local
    # directories -- passing the gs:// puck_path as the output directory makes Julia try to create a
    # literal 'gs:' directory instead.
    if STAGE_GCS:
        puck_dir = PUCK_DEST
        barcodes_dir = BARCODES_DEST
        puck_dir.mkdir(exist_ok=True, parents=True)
        barcodes_dir.mkdir(exist_ok=True, parents=True)
    else:
        puck_dir = Path(PUCK_PATH)
        barcodes_dir = Path(RAW_BARCODES_PATH) if RAW_BARCODES_PATH is not None else None

    env["LIB_PUCK_PATH"] = str(puck_dir)
    env["LIB_PUCK_IN"] = str(barcodes_dir)

    # for each sample, ensure the correct puck is present and run the spatial barcode script
    for _, sample in metadata_df.iterrows():
        sample_name = str(sample['Sample Name'])
        if 'Add Puck ID' in sample and not pd.isna(sample['Add Puck ID']):
            puck_id = sanitize_path_component(sample['Add Puck ID'], "Add Puck ID")
        else:
            puck_id = sanitize_path_component(sample['Puck ID'], "Puck ID")
        # Pick the spatial-barcode FASTQ directory from --fastqs, not from which stage flags were
        # passed: `--spatial-count` on its own is a perfectly normal way to re-run this stage
        # against FASTQs mkfastq already produced, and the old flag-based test sent that case to
        # the BCL input directory instead.
        if FASTQ_INPUT is not None:
            # spatial_count.jl selects files by substring-matching the sample ID inside a directory,
            # so a mixed --fastqs directory has to be narrowed to this sample's spatial library
            # first or its GEX reads would be counted as spatial reads
            fastq_dir = stage_spatial_fastqs(FASTQ_INPUT, sample_name, declared_lanes(sample, 'SB Lane'))
        else:
            fastq_dir = MKFASTQ_OUTS / f"{sample_name}_sb"
            if not fastq_dir.is_dir() or not any(fastq_dir.glob('*.fastq.gz')):
                log_write(f"[ERROR]: no spatial-barcode FASTQs found for sample {sample_name} at {fastq_dir}")
                log_write("Troubleshooting:")
                log_write(" • Run the demultiplexing stage first with --mkfastq (or the whole pipeline with --run-all)")
                log_write(" • If the reads were demultiplexed elsewhere, pass that directory with --fastqs")
                sys.exit(1)
        log_write(f"  Processing sample {sample_name}...", terminal=False)

        # create output sample directory
        sample_output = SPATIAL_COUNT_OUTS / sample_name
        sample_output.mkdir(parents=True, exist_ok=True)

        # if a puck file already exists, copy it into the working directory
        input_puck_file = puck_dir / f'{puck_id}.csv'
        if STAGE_GCS and not input_puck_file.is_file():
            # A puck map absent from the bucket is NOT fatal: fall through to generating it from raw
            # barcodes, exactly as a local run does when the CSV is missing. Exiting here instead
            # would make the documented fallback unreachable.
            stage_from_gcs(
                f'{gcs_uri(PUCK_PATH)}/{puck_id}.csv',
                input_puck_file,
                required=False,
                description=f'puck map {puck_id}.csv'
            )
        if input_puck_file.is_file():
            log_write(f"  Using existing puck file: {input_puck_file}\n")
            shutil.copy(input_puck_file, sample_output / f'{puck_id}.csv')
        else:
            # run the corresponding script to generate the puck file
            log_write(f"  Puck file {input_puck_file} not found. Generating from raw barcodes... ", terminal=False, terminator="")
            if barcodes_dir is None:
                log_write(f"\n[ERROR]: puck map {puck_id}.csv is unavailable and `raw_barcodes_path` is not set, "
                          f"so it cannot be generated")
                log_write("Troubleshooting:")
                log_write(f" • Set `raw_barcodes_path` to the location of {puck_id}'s BeadBarcodes.txt/BeadLocations.txt")
                log_write(f" • Or place a pre-built {puck_id}.csv in `puck_path`")
                sys.exit(1)
            if STAGE_GCS:
                stage_from_gcs(
                    f'{gcs_uri(RAW_BARCODES_PATH)}/{puck_id}',
                    barcodes_dir,
                    recursive=True,
                    description=f'raw barcodes for puck {puck_id}'
                )

            with open(spatial_barcodes_log, "a") as barcodes_logfile:
                with console.status(f"  Puck file {input_puck_file} not found. Generating from raw barcodes... "):
                    proc = subprocess.Popen(
                        [
                            julia_interp_path,
                            script_base_path / 'generate_puck_csv.jl'
                        ],
                        stdout=barcodes_logfile,
                        stderr=subprocess.STDOUT,
                        env=env
                    )

                    proc.wait()

            if proc.returncode == 0:
                log_ts(f"generated puck {puck_id}")
            else:
                job_crash("generate_puck_csv.jl", proc.returncode, spatial_barcodes_log)

        # run the spatial barcode count script with modified environment
        log_write(f"  Processing sample {sample_name}... ", terminal=False, terminator="")
        with open(spatial_barcodes_log, "a") as barcodes_logfile:
            with console.status(f"  Processing sample {sample_name}... "):
                proc = subprocess.Popen(
                    [
                        julia_interp_path,
                        script_base_path / 'spatial_count.jl',
                        fastq_dir,
                        sample_output,
                        sample_name
                    ],
                    stdout=barcodes_logfile,
                    stderr=subprocess.STDOUT,
                    env=env
                )

                proc.wait()

        # check exit status
        if proc.returncode == 0:
            log_ts(f"counted spatial barcodes for {sample_name}")
        else:
            job_crash("spatial_count.jl", proc.returncode, spatial_barcodes_log)
    
    # log job execution result
    job_success("spatial barcode counting", spatial_barcodes_log, SPATIAL_COUNT_OUTS, started=stage_started)


def run_spatial_analysis() -> None:
    """
    Run the spatial analysis scripts on samples in spatial_barcodes
    """

    # load metadata CSV and spatial analysis scripts and create 
    global proc
    script_base_path = SCRIPT_PATH / 'spatial_analysis'
    spatial_analysis_log = LOG_PATH / "spatial_analysis.log"
    SPATIAL_ANALYSIS_OUTS.mkdir(exist_ok=True)
    metadata_df = pd.read_csv(SUMMARY_PATH)

    # top N percent of beads (by total UMI count) to filter out as bead aggregates/contamination
    percent_umi_filtering = PERCENT_UMI_FILTERING
    if percent_umi_filtering is None:
        percent_umi_filtering = 1
    elif not is_number(percent_umi_filtering):
        log_write(f"[WARNING]: top_n_percent_umi_filter must be numeric, but is currently set to {percent_umi_filtering!r}")
        log_write("Defaulting to the top 1% of beads by UMI count for this run")
        log_write(" • Set `workflow.top_n_percent_umi_filter` to a number between 0 and 100 (e.g. 1), unquoted")
        log_write(" • Remove the field entirely to accept the default of 1 without this warning")
        hint = bool_value_hint(percent_umi_filtering, 'workflow.top_n_percent_umi_filter')
        if hint:
            log_write(hint)
        percent_umi_filtering = 1

    if SPATIAL_DOWNSAMPLING is not None:
        if not isinstance(SPATIAL_DOWNSAMPLING, float):
            log_write(f"[WARNING]: spatial downsampling rate must be a float, but is currently set to {SPATIAL_DOWNSAMPLING!r}")
            log_write(f"Spatial downsampling will be ignored for this run (this may cause an OOM crash if your samples are sequenced deeply)")
            log_write(" • Set `workflow.spatial_downsampling` to a decimal fraction of reads to keep, e.g. 0.5 -- note that `1` parses as an integer, so write `1.0`")
            log_write(" • Remove the field entirely to run on all reads without this warning")
        else:
            for sample in metadata_df['Sample Name'].astype(str).values:
                downsample_spatial(
                    SPATIAL_COUNT_OUTS / sample / "SBcounts.h5",
                    SPATIAL_COUNT_OUTS / sample / f"SBcounts_downsampled_{str(SPATIAL_DOWNSAMPLING).replace('.', '')}.h5"
                )

    stage_started = log_stage_start("spatial analysis")

    # activate conda environment and install R libraries
    conda_lib = ensure_conda_env("r_env")

    # export environment variables
    env = conda_subprocess_env()
    env["R_FUNC"] = str(SCRIPT_PATH / 'spatial_analysis' / 'functions')
    env["R_LIBS"] = str(conda_lib)
    env["DATA_PATH"] = str(OUTPUT_PATH.parent.parent)
    # the R scripts take the directory *containing* reference genomes, not the genome itself
    env["REF_PATH"] = str(staged_ref_dir().parent)

    # define samples list
    samples = metadata_df['Sample Name'].astype(str).tolist()
    sample_list = "[" + ", ".join(samples) + "]"
    log_write([
        f"  Using R environment at {conda_lib}",
        f"  Samples: {sample_list}\n"
    ])

    # check if spatial analysis script should use Cellbender or regular .h5 files
    if not CELLBENDER_OUTS.is_dir() or len(list(CELLBENDER_OUTS.iterdir())) == 0:
        log_write(f"  Cellbender outputs not found. Using direct cellranger count outputs in {COUNT_OUTS}\n")
        use_cellbender = "FALSE"
    elif not args.no_cellbender:
        log_write(f"  Using denoised cellbender outputs at {CELLBENDER_OUTS} (provide the --no-cellbender flag to force the spatial analysis to use direct cellranger count .h5 files)\n")
        use_cellbender = "TRUE"
    else:
        log_write(f"  Using direct cellranger count .h5 files at {COUNT_OUTS} (--no-cellbender flag provided)\n")
        use_cellbender = "FALSE"

    log_write("  Running spatial analysis... ", terminal=False, terminator="")

    # run the spatial analysis R scripts with the modified environment
    with open(LOG_PATH / 'spatial_analysis.log', "w") as logfile:
        with console.status("  Running spatial analysis... "):
            proc = subprocess.Popen(
                [
                    'stdbuf', '-oL', '-eL', 'mamba', 'run', '-n', 'r_env', 'Rscript',
                    script_base_path / 'run_spatial.R',
                    BCL_ID,
                    sample_list,
                    use_cellbender,
                    str(SPATIAL_DOWNSAMPLING).replace('.', '') if SPATIAL_DOWNSAMPLING is not None else "",
                    str(NUM_THREADS),
                    str(percent_umi_filtering)
                ],
                stdout=logfile,
                stderr=subprocess.STDOUT,
                env=env
            )

            proc.wait()
            log_write('Done.')

    # check exit status
    if proc.returncode == 0:
        job_success("spatial analysis", spatial_analysis_log, SPATIAL_ANALYSIS_OUTS, started=stage_started)
    else:
        job_crash("run_spatial.R", proc.returncode, spatial_analysis_log)


def run_takara_spatial_profiling() -> None:
    """
    Run custom scripts from Takara for 10x Flex sample processing
    """

    stage_started = log_stage_start("spatial analysis (Flex / Trekker)")

    # activate Trekker conda environment, load scripts, and create required subdirectories
    # No definition file: the Trekker environment is not built by slidr (this repository carries no
    # trekker env spec), so it must already exist. Passing None makes a missing environment report
    # that -- and point at ~/.condarc's envs_dirs, which is how a shared prebuilt env is picked up --
    # rather than complaining about a definition file that was never shipped.
    trekker_env = ensure_conda_env("trekker", environment_yml=None)
    takara_pipeline_log = LOG_PATH / "takara_pipeline.log"
    takara_path = SCRIPT_PATH / "takara"
    flex_outputs_path = OUTPUT_PATH / "flex"
    flex_puck_path = flex_outputs_path / "pucks"
    flex_samplesheets_path = flex_outputs_path / "samplesheets"
    flex_outputs_path.mkdir(exist_ok=True)
    flex_puck_path.mkdir(exist_ok=True)
    flex_samplesheets_path.mkdir(exist_ok=True)
    metadata_df = pd.read_csv(SUMMARY_PATH)

    # add the Trekker conda environment to PATH
    env = conda_subprocess_env()
    env['PATH'] = f"{trekker_env}:{env.get('PATH', '')}"
    env.pop('PYTHONPATH', None)

    # generate demux samplesheet
    # Build parallel lists (one row per (sample, spatial barcode) pair) rather than a dict keyed on
    # the spatial barcode. Keying on the barcode silently dropped a sample whenever the same probe
    # barcode ID appeared under two samples; emitting every pair means a genuinely reused barcode
    # instead surfaces as trekker_demux.py's "Duplicate barcode label" error rather than a silent
    # drop. `demux_samples` is the ordered list of partition names reused by the loops below.
    demux_samples = []   # column 1: TrekkerFX_<sample>_<AB> partition names
    demux_barcodes = []  # column 2: corresponding <AB> spatial barcodes
    for _, sample in metadata_df.iterrows():
        barcode_group = sample['Flex Probe Barcode IDs']
        sample_barcodes = split_probe_barcodes(barcode_group)
        if not sample_barcodes:
            log_write(f"[ERROR]: no value provided for the `Flex Probe Barcode IDs` metadata field for sample {sample['Sample Name']}")
            log_write("Troubleshooting:")
            log_write(" • Fill in the `Flex Probe Barcode IDs` column for every Flex sample -- the Trekker demultiplexer needs it to split the spatial reads")
            log_write(" • Use the probe barcode IDs from your Flex kit (e.g. BC001), separated by ',' or '|' when a sample carries more than one")
            log_write(" • If this sample is not actually Flex, correct its `Chemistry` column instead")
            log_write(f" • Metadata source: {METADATA_SRC}")
            sys.exit(1)
        if len(sample_barcodes) == 1:
            log_write(f"[WARNING]: only one probe barcode parsed from the `Flex Probe Barcode IDs` value for sample "
                      f"{sample['Sample Name']} ('{barcode_group}'); treating it as a single barcode")
            log_write(" • If this sample carries several probe barcodes, separate them with ',' or '|' (e.g. `BC001,BC002` or `BC001|BC002`)")
        sample_name = sample['Sample Name']
        for barcode in sample_barcodes:
            spatial_barcode = barcode.replace('BC', 'AB')
            demux_samples.append(f"TrekkerFX_{sample_name}_{spatial_barcode}")
            demux_barcodes.append(spatial_barcode)

    demux_samplesheet = pd.DataFrame({
        'samples': demux_samples,
        'barcodes': demux_barcodes
    })
    demux_samplesheet.to_csv(flex_samplesheets_path / 'Trekker_demux_samplesheet.csv', header=False, index=False)
    log_write(f"  Generated Trekker demultiplexing samplesheet: {flex_samplesheets_path / 'Trekker_demux_samplesheet.csv'}")

    # run trekker demultiplexer
    log_write(f"  Running Trekker demultiplexer... ", terminal=False, terminator="")
    with open(takara_pipeline_log, "w") as logfile:
        with console.status("  Running Trekker demultiplexer... "):
            proc = subprocess.Popen(
                [
                    'mamba', 'run', '-n', 'trekker', 'python',
                    takara_path / 'demultiplexing' / 'trekker_demux.py',
                    FLEX_R1_PATH,
                    FLEX_R2_PATH,
                    flex_samplesheets_path / 'Trekker_demux_samplesheet.csv',
                    flex_outputs_path / 'demux'
                ],
                stdout=logfile,
                stderr=subprocess.STDOUT,
                env=env
            )

            proc.wait()

    if proc.returncode == 0:
        log_write("Done.")
    else:
        job_crash("trekker_demux", proc.returncode, takara_pipeline_log)

    # retrieve bead barcode file
    tile_ids = [sanitize_path_component(t, "Puck ID") for t in metadata_df['Puck ID'].tolist()]
    puck_paths = []
    for tile_id in tile_ids:
        if not Path(flex_puck_path / f"{tile_id}_BeadBarcodes.csv").is_file():
            puck_path = retrieve_takara_bead_barcode_file(tile_id, flex_puck_path)
            with zipfile.ZipFile(puck_path, "r") as zf:
                zf.extractall(puck_path.parent)
            puck_path = str(puck_path).replace(".zip", ".txt")
            puck_file = pd.read_table(puck_path)
            puck_path = str(puck_path).replace(".txt", ".csv")
            puck_paths.append(puck_path)
            puck_file.to_csv(puck_path, header=False, index=False)
            log_write(f"  Downloaded puck file {puck_path} from Takeda website")
        else:
            puck_path = Path(flex_puck_path / f"{tile_id}_BeadBarcodes.csv")
            log_write(f"  Using existing puck file: {puck_path}")

    # extract experiment date
    date = BCL_ID.split('_')[0]
    try:
        datetime.strptime(date, "%Y%m%d")
    except ValueError:
        # find another workaround in the long run
        date = datetime.now().date()
        date = date.strftime("%Y%m%d")

    # generate pipeline samplesheet
    pipeline_samplesheet = {
        'sample': [],
        'sc_sample': [],
        'experiment_date': [],
        'barcode_file': [],
        'fastq_1': [],
        'fastq_2': [],
        'sc_outdir': [],
        'sc_platform': [],
        'profile': [],
        'subsample': [],
        'cores': []
    }
    for sample in demux_samples:
        spatial_barcode = sample.split('_')[-1]
        barcode = sample.replace('TrekkerFX_', '')
        # reset per iteration: a sample matching no metadata row must error here rather than
        # silently inherit the previous iteration's (or the download loop's leftover) puck_path
        puck_path = None
        for _, row in metadata_df.iterrows():
            if f"TrekkerFX_{row['Sample Name']}_{spatial_barcode}" == sample:
                puck_path = flex_puck_path / f"{sanitize_path_component(row['Puck ID'], 'Puck ID')}_BeadBarcodes.csv"
                break
        if puck_path is None:
            log_write(f"[ERROR]: no metadata row matches spatial barcode partition '{sample}'; cannot assign a puck file")
            log_write("Troubleshooting:")
            log_write(" • Partition names are built as TrekkerFX_<Sample Name>_<probe barcode with BC replaced by AB>, so this means the sample name or barcode changed mid-run")
            log_write(f" • Check the `Sample Name` and `Flex Probe Barcode IDs` columns are unchanged since this run started: {SUMMARY_PATH}")
            log_write(" • Re-run the spatial analysis stage so the partitions are rebuilt from the current metadata")
            log_write(f" • Metadata source: {METADATA_SRC}")
            sys.exit(1)
        demux_output = flex_outputs_path / 'demux'
        count_output = COUNT_OUTS / 'flex' / barcode.replace('AB', 'BC') / 'count' / 'sample_filtered_feature_bc_matrix'
        pipeline_samplesheet['sample'].append(sample)
        pipeline_samplesheet['sc_sample'].append(sample.replace('_AB', '_scRNAseq_AB'))
        pipeline_samplesheet['experiment_date'].append(date)
        pipeline_samplesheet['barcode_file'].append(puck_path)
        pipeline_samplesheet['sc_outdir'].append(count_output)
        pipeline_samplesheet['fastq_1'].append(demux_output / f"{sample}_R1.fastq.gz")
        pipeline_samplesheet['fastq_2'].append(demux_output / f"{sample}_R2.fastq.gz")
        pipeline_samplesheet['sc_platform'].append('TrekkerFX_FLEX')
        pipeline_samplesheet['profile'].append('conda')
        pipeline_samplesheet['subsample'].append('no')
        pipeline_samplesheet['cores'].append(str(NUM_THREADS))

    pipeline_samplesheet = pd.DataFrame(pipeline_samplesheet)
    pipeline_samplesheet.to_csv(flex_samplesheets_path / 'Trekker_flex_samplesheet.csv', index=False)
    log_write(f"Generated Trekker pipeline samplesheet: {flex_samplesheets_path / 'Trekker_flex_samplesheet.csv'}\n")

    # run the Takara pipeline on each partition
    (flex_outputs_path / 'trekker').mkdir(exist_ok=True)
    for _, sample in pipeline_samplesheet.iterrows():
        log_write(f"Processing sample [{sample['sample']}]... ")
        sub_samplesheet = pd.DataFrame(
            columns=pipeline_samplesheet.columns,
            data = [sample.tolist()]
        )

        (flex_samplesheets_path / 'samples').mkdir(exist_ok=True)
        sub_samplesheet_path = flex_samplesheets_path / 'samples' / f"{sample['sample']}.csv"
        sub_samplesheet.to_csv(sub_samplesheet_path, index=False)

        # run the pipeline script
        with open(LOG_PATH / 'takara_pipeline.log', "a") as logfile:
            proc = subprocess.Popen(
                [
                    'mamba', 'run', '-n', 'trekker', 'bash',
                    takara_path / 'profiling' / 'nuclei_locator_wrapper.sh',
                    sub_samplesheet_path,
                    takara_path / 'profiling',
                    flex_outputs_path / 'trekker',
                    trekker_env.parent
                ],
                stdout=logfile,
                stderr=subprocess.STDOUT,
                env=env
            )

            proc.wait()

        if proc.returncode == 0:
            log_write("Done.\n")
        else:
            job_crash("trekker_flex", proc.returncode, takara_pipeline_log)

    # generate the merge samplesheet
    merge_samplesheet = {
        'sample': [],
        'out_dir': []
    }
    for sample in demux_samples:
        merge_samplesheet['sample'].append(sample)
        trekker_output_path = flex_outputs_path / "trekker" / f"{date}_{sample}" / f"trekker_{sample}" / "output"
        merge_samplesheet['out_dir'].append(trekker_output_path)
    
    merge_samplesheet = pd.DataFrame(merge_samplesheet)
    merge_samplesheet.to_csv(flex_samplesheets_path / "Trekker_merge_samplesheet.csv", index=False)
    log_write(f"Generated Trekker merge samplesheet: {flex_samplesheets_path / 'Trekker_merge_samplesheet.csv'}\n")

    log_write('Running Trekker merge...')
    # merge partitions for each sample
    for _, sample in metadata_df.iterrows():
        log_write(f"  Processing sample [{sample['Sample Name']}]... ")
        # run the Takara merger
        with open(LOG_PATH / 'takara_pipeline.log', "a") as logfile:
            proc = subprocess.Popen(
                [
                    'mamba', 'run', '-n', 'trekker', 'bash',
                    takara_path / 'merging' / 'trekker_merger.sh',
                    flex_samplesheets_path / 'Trekker_merge_samplesheet.csv',
                    flex_outputs_path,
                    sample['Sample Name'],
                    'conda',
                    takara_path / 'merging',
                    str(trekker_env.parent)
                ],
                stdout=logfile,
                stderr=subprocess.STDOUT,
                env=env
            )

            proc.wait()

        if proc.returncode == 0:
            log_write("Done.\n")
        else:
            job_crash("trekker_merge", proc.returncode, takara_pipeline_log)