import yaml
import os
import sys
import time
import argparse
import tomllib
import subprocess
import tempfile

from datetime import datetime
from pathlib import Path
from rich.console import Console
from importlib.metadata import version, PackageNotFoundError


# initialize rich consoles for text formatting
err_console = Console(stderr=True)
console = Console()


def _timeout_from_env(name: str, default: int) -> int:
    """
    read a positive-integer timeout (in seconds) from the environment, falling back to `default`

    inputs:
     - name:    the name of the corresponding timeout environment variable
     - default: the default number of seconds for the timeout
    output:
     - the resolved/validated timeout limit in seconds
    """

    # retrieve the named environment variable and validate that it's a positive integer
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return default
    try:
        value = int(raw.strip())
    except ValueError:
        err_console.print(f"[yellow]\\[WARNING][/yellow]: {name} is not an integer ('{raw.strip()}'); using the default of {default}s")
        return default
    if value <= 0:
        err_console.print(f"[yellow]\\[WARNING][/yellow]: {name} must be a positive number of seconds (got {value}); using the default of {default}s")
        return default

    return value


# define timeout limits for GCS metadata and data transfer (in seconds)
GCS_METADATA_TIMEOUT = _timeout_from_env('SLIDR_GCS_METADATA_TIMEOUT', 5 * 60)
GCS_TRANSFER_TIMEOUT = _timeout_from_env('SLIDR_GCS_TRANSFER_TIMEOUT', 12 * 60 * 60)

# define timeout limits for Julia installation (in seconds)
INSTALLER_FETCH_TIMEOUT = _timeout_from_env('SLIDR_INSTALLER_FETCH_TIMEOUT', 5 * 60)
INSTALLER_RUN_TIMEOUT = _timeout_from_env('SLIDR_INSTALLER_RUN_TIMEOUT', 30 * 60)

# sentinel stored by `--fastqs` when it is passed without a directory, meaning "the FASTQs are in
# `paths.input_path`" - distinct from None, which means `--fastqs` was not passed at all
FASTQS_USE_INPUT_PATH = '<input_path>'

# host-specific locations an environment variable may override, e.g. slidr_slurm.sh pointing
# output_path at a workdir chosen at submit time (everything else comes only from config.yaml)
PATH_ENV_OVERRIDES = {
    'input_path': 'SLIDR_INPUT_PATH',
    'output_path': 'SLIDR_OUTPUT_PATH',
    'auth_key_path': 'SLIDR_AUTH_KEY_PATH',
    'software_path': 'SLIDR_SOFTWARE_PATH',
}

SETTINGS_ENV_OVERRIDES = {
    'gcs_download_dest': 'SLIDR_GCS_DOWNLOAD_DEST',
}

# define filesystem types that do not support the pipeline's logging pattern
# (opening, appending to, and closing the logfile repeatedly)
NETWORK_FSTYPES = frozenset({'cifs', 'smb3', 'smbfs'})

def get_version():
    """
    Retrieve project version
    """
    try:
        # attempt to read from installed package metadata
        return version("slidr")
    except PackageNotFoundError:
        # attempt to read directly from pyproject.toml if running uninstalled
        pyproject_path = Path(__file__).parent.parent.parent / "pyproject.toml"
        if pyproject_path.exists():
            try:
                # read package version from pyproject.toml
                with open(pyproject_path, 'rb') as f:
                    data = tomllib.load(f)
                    # safely look up "version" key in loaded dict
                    return data.get("project", {}).get("version", "unknown")
            except Exception as e:
                return "unknown"
        else:
            err_console.print("[ERROR]: pyproject.toml file not found in project root")
            err_console.print("Troubleshooting:")
            err_console.print(" • Run `git pull` to restore the missing pyproject.toml file")
            sys.exit(1)


def filesystem_type(path) -> str | None:
    """
    return the type of the filesystem backing `path`, or None if it cannot be determined

    inputs:
     - path: the path to identify the backing filesystem of
    output:
     - the filesystem type as the kernel names it ('ext4', 'cifs', 'nfs4', ...), or None where that
       cannot be read, which includes macOS and any platform without /proc
    """

    # resolve the path to absolute and follow symlinks
    try:
        target = Path(path).resolve()
    except OSError:
        return None

    # try to read the table of mounted filesystems
    mounts = []
    try:
        with open('/proc/self/mountinfo', 'r') as handle:
            for line in handle:
                # mountinfo splits at ' - ': the mount point is field 5 on the left of that
                # separator, the filesystem type the first field on the right
                left, separator, right = line.partition(' - ')
                if not separator:
                    continue
                left_fields, right_fields = left.split(), right.split()
                if len(left_fields) < 5 or not right_fields:
                    continue
                mounts.append((left_fields[4], right_fields[0]))
    except OSError:
        return None

    # find the longest mount point containing the path (e.g., `/mnt/share` must beat `/`);
    # `>=` keeps the later of two entries at one path, e.g. cifs mounted over its autofs trigger
    best_length, best_type = -1, None
    for mountpoint, fstype in mounts:
        try:
            target.relative_to(mountpoint)
        except ValueError:
            continue
        if len(mountpoint) >= best_length:
            best_length, best_type = len(mountpoint), fstype

    # return the filesystem type of the longest match
    return best_type


def resolve_log_staging(log_path: Path, bcl_id: str) -> tuple[Path, Path | None]:
    """
    decide where this run writes its logs, staging them on local disk when the run directory is on a
    network filesystem that cannot take the append pattern (see NETWORK_FSTYPES)

    inputs:
     - log_path: the run's real log directory, <OUT_PATH>/log
     - bcl_id:   this run's BCL ID, used to name the staging directory
    output:
     - (write_path, final_path): where the run should write logs, and where they must be copied at
       exit - final_path is None when no staging is happening, which is the local-disk case and
       leaves every existing path untouched
    """

    # user can turn off local staging function
    if os.environ.get('SLIDR_NO_LOG_STAGING') == '1':
        return log_path, None

    # identify the filesystem type of the mount at which output directory is located
    fstype = filesystem_type(log_path.parent)
    if fstype not in NETWORK_FSTYPES:
        return log_path, None

    # one staging directory per process, not per BCL
    # otherwise two runs of the same BCL can overlap and write to the same staged log
    override = os.environ.get('SLIDR_LOG_STAGE_DIR')
    base = Path(override) if override else Path(tempfile.gettempdir())
    staged = base / f'slidr-logs-{bcl_id}-{os.getpid()}'

    # try to create a staged logfile output directory and fall back to default option if can't
    try:
        staged.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        console.print(f"[yellow]\\[WARNING][/yellow]: {log_path.parent} is a {fstype} mount, where slidr's "
                      f"logging is unreliable, but the staging directory {staged} could not be created: {error}")
        console.print(f"  Logs will be written directly to {log_path}; a run may fail with "
                      f"'[Errno 9] Bad file descriptor'")
        console.print(f"  Set SLIDR_LOG_STAGE_DIR to a writable local directory to stage them instead")
        return log_path, None

    # check if the staging directory is itself on an unsupported network mount
    staged_fstype = filesystem_type(staged)
    if staged_fstype in NETWORK_FSTYPES:
        console.print(f"[yellow]\\[WARNING][/yellow]: the log staging directory {staged} is itself on a "
                      f"{staged_fstype} mount, so staging would not avoid the problem it exists for")
        console.print(f"  Logs will be written directly to {log_path}")
        console.print(f"  Point SLIDR_LOG_STAGE_DIR at local disk, or set SLIDR_NO_LOG_STAGING=1 to silence this")
        try:
            staged.rmdir()
        except OSError:
            pass
        return log_path, None

    return staged, log_path


def is_gcs_path(value) -> bool:
    """
    test whether a configured value explicitly names a Google Cloud Storage object

    inputs:
     - value: the configured value to check
    output:
     - True if the value is a gs:// URI
    """

    return isinstance(value, str) and value.strip().startswith('gs://')


def is_google_sheet(source) -> bool:
    """
    test whether a metadata source is a Google Sheet URL rather than a local file

    inputs:
     - source: the configured or overridden metadata source
    output:
     - True if the source is a Google Sheets URL
    """

    return isinstance(source, str) and source.strip().startswith("https://docs.google.com/spreadsheets")


def is_int(value) -> bool:
    """
    test whether a configuration value is a usable integer
    (to ensure parameter values are parsed correctly by PyYAML)

    inputs:
     - value: the configuration value to check
    output:
     - True if the value is an int and not a bool
    """

    return isinstance(value, int) and not isinstance(value, bool)


def is_number(value) -> bool:
    """
    whether a configuration value is a usable number (int or float), rejecting booleans

    inputs:
     - value: the configuration value to check
    output:
     - True if the value is an int or float and not a bool
    """

    return isinstance(value, (int, float)) and not isinstance(value, bool)


def bool_value_hint(value, field: str) -> str | None:
    """
    extra troubleshooting line for a numeric field that was given a boolean

    inputs:
     - value: the offending configuration value
     - field: dotted config field name, for the message
    output:
     - a troubleshooting bullet, or None when `value` is not a boolean
    """

    if not isinstance(value, bool):
        return None

    return (f" • `{field}` takes a number, not a yes/no switch - YAML reads the bare words"
            " true/false/yes/no/on/off as booleans, so quote the value if you meant a number")


def gcs_uri(path: Path | str) -> str:
    """
    normalize a configured Google Cloud Storage location to a `gs://` URI.

    inputs:
     - path: configured GCS location, with or without the scheme
    output:
     - the location as a `gs://`-prefixed URI with no trailing slash
    """

    text = str(path).strip().rstrip('/')
    if not text.startswith('gs://'):
        text = 'gs://' + text.lstrip('/')

    return text


def _first_token_line(text: str) -> str:
    """
    extract a secret from the contents of a token file
    (the first line that is neither blank nor a `#` comment)

    inputs:
     - text: full contents of the token file
    output:
     - the token, or '' if the contents hold none
    """

    lines = (line.strip() for line in text.splitlines())
    return next((line for line in lines if line and not line.startswith('#')), '')


def _read_gcs_token(uri: str) -> str:
    """
    read a secret out of a GCS object with `gcloud storage cat`, without writing it to disk

    inputs:
     - uri: gs:// URI of the object holding the token
    output:
     - the token
    """

    # assemble and run the gcloud command
    cmd = ['gcloud', 'storage', 'cat', uri]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=GCS_METADATA_TIMEOUT)
    except FileNotFoundError:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: `slack_token` is a gs:// URI but `gcloud` was not found on PATH, so it cannot be read")
        err_console.print("Troubleshooting:")
        err_console.print(" • Install the Google Cloud CLI: https://cloud.google.com/sdk/docs/install")
        err_console.print(" • On a cluster, check whether it needs loading first (e.g. `module load google-cloud-sdk`)")
        err_console.print(" • Or point `slack_token` at a local file, or disable alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: reading the Slack token from {uri} timed out after {GCS_METADATA_TIMEOUT}s")
        err_console.print("Troubleshooting:")
        err_console.print(" • Reading one small object should take seconds, so this usually means gcloud cannot reach the network at all")
        err_console.print(f" • Check the object is reachable by hand: `gcloud storage cat {uri}`")
        err_console.print(" • Raise the limit with SLIDR_GCS_METADATA_TIMEOUT=<seconds> if the link is genuinely this slow")
        err_console.print(" • Or point `slack_token` at a local file, or disable alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    # raise an error if gcloud process errored out; quote gcloud's stderr but never the object's
    # contents, so the secret cannot leak into a shared log
    if proc.returncode != 0:
        detail = (proc.stderr or '').strip().splitlines()
        err_console.print(f"[bold red]\\[ERROR][/bold red]: could not read the Slack token from {uri} (`gcloud storage cat` exited {proc.returncode})")
        for line in detail[-3:]:
            err_console.print(f"  {line}")
        err_console.print("Troubleshooting:")
        err_console.print(f" • Check the object exists: `gcloud storage ls {uri}`")
        err_console.print(" • Check the active account can read it: `gcloud auth print-access-token`")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    # raise an error if the returned token is empty or misconfigured (not in the first line of the file)
    token = _first_token_line(proc.stdout)
    if not token:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: the Slack token object {uri} contains no token")
        err_console.print("Troubleshooting:")
        err_console.print(" • The object should hold the bot token on its own line (lines starting with `#` are ignored)")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    return token


def _resolve_slack_token(value) -> str:
    """
    resolve the configured `slack_token` into an actual Slack bot token.

    inputs:
     - value: the raw `settings.slack_token` value
    output:
     - the token, stripped of surrounding whitespace
    """
    # exit if Slack alerts are enabled but no token is provided
    if not isinstance(value, str) or not value.strip():
        err_console.print(f"[bold red]\\[ERROR][/bold red]: `alerts` is enabled but `slack_token` is not set: {value!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Point `slack_token` at a file containing the bot token (preferred, so the token stays out of the config file)")
        err_console.print(" • Or set it to the token itself")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    value = value.strip()

    # stream the token value from a GCS bucket without saving it locally
    if is_gcs_path(value):
        return _read_gcs_token(value)

    # expand the ~ shorthand to an absolute directory path
    candidate = Path(value).expanduser()

    # if the provided value is a filename, try to read the first line
    if not candidate.is_file():
        looks_like_path = value.startswith(('/', '~', './', '../')) or '/' in value
        if looks_like_path:
            err_console.print(f"[bold red]\\[ERROR][/bold red]: `slack_token` looks like a file path, but no such file exists: {candidate}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Check the path for typos and that the file is readable by the user running the pipeline")
            err_console.print(" • On a --gcp or staged --slurm run the token file must exist on the compute node, not just your laptop")
            err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
            sys.exit(1)
        return value
    try:
        token = _first_token_line(candidate.read_text())
    except OSError as exc:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: could not read the Slack token from {candidate}: {exc}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check the file is readable by the user running the pipeline")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    # raise an error if the provided token file is empty or misconfigured
    if not token:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: the Slack token file {candidate} contains no token")
        err_console.print("Troubleshooting:")
        err_console.print(" • The file should hold the bot token on its own line (lines starting with `#` are ignored)")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    # raise a warning if slack token permissions are too broad
    mode = candidate.stat().st_mode
    if mode & 0o077:
        err_console.print(f"[bold yellow]\\[WARNING][/bold yellow]: the Slack token file {candidate} is readable by other users (mode {mode & 0o777:03o})")
        err_console.print(f" • Restrict it with `chmod 600 {candidate}`")

    return token


def _validate_resource_path(value, field: str, contents: str, required: bool = True):
    """
    check if the provided resource path is valid

    inputs:
     - value:     user-provided value for the resource path
     - field:     resource path setting name in the configfile
     - contents:  a string describing the required resource contents
     - required:  Boolean flag stating whether the resource is required for operation in a given mode
    outputs:
     - a `gs://` URI string for a bucket, a local Path otherwise, or None for an unset optional field

    """

    # check if the resource path is set in the configfile
    if value is None or (isinstance(value, str) and not value.strip()):
        if not required:
            return None
        err_console.print(f"[bold red]\\[ERROR][/bold red]: `{field}` is not set in the configuration file")
        err_console.print("Troubleshooting:")
        err_console.print(f" • Set `{field}` to a local directory containing {contents}")
        err_console.print(f" • Or to the GCS location holding it (e.g. gs://my-bucket/prefix), which is staged automatically")
        sys.exit(1)

    # validate the GCS path if the resource is a GCS bucket
    if is_gcs_path(value):
        return gcs_uri(value)

    # if the resource is not a GCS bucket, check if it's a local directory
    if not Path(value).is_dir():
        err_console.print(f"[bold red]\\[ERROR][/bold red]: `{field}` in the configuration file is not a directory: {value}")
        err_console.print("Troubleshooting:")
        err_console.print(f" • Check that `{field}` is set to a valid directory containing {contents}")
        err_console.print(f" • If this data lives in a bucket, write the full URL with its scheme: gs://{str(value).lstrip('/')}")
        err_console.print(" • A bucket is recognised by the `gs://` prefix alone -- there is no flag to pass")
        sys.exit(1)

    return Path(value)


def _parse_arguments() -> argparse.Namespace:
    """
    parse command-line arguments

    outputs:
     - an argparse.Namespace object containing command-line arguments

    """

    parser = argparse.ArgumentParser()

    # not required=True, which would break `--version`; presence is checked in _load() instead
    parser.add_argument(
        "--bcl",
        help="BCL ID of your run",
        type=str
    )

    parser.add_argument(
        "--config",
        "-cf",
        type=str,
        help="Path to a custom YAML configuration file"
    )

    # The path is optional: a bare --fastqs means "`paths.input_path` holds the FASTQs"
    parser.add_argument(
        "--fastqs",
        "-fq",
        type=str,
        nargs='?',
        const=FASTQS_USE_INPUT_PATH,
        help="Treat the input as already-demultiplexed FASTQ files and skip cellranger mkfastq. "
             "Takes an optional directory; with no value, `paths.input_path` is used"
    )

    parser.add_argument(
        "--metadata",
        "-md",
        type=str,
        help="Path to a .tsv/.csv metadata file, or a Google Sheet URL, overriding "
             "settings.metadata_source in the configuration file"
    )

    parser.add_argument(
        "--mkfastq",
        "-mf",
        action="store_true",
        help="Run cellranger mkfastq (if running individual pipeline components)"
    )

    parser.add_argument(
        "--count",
        "-ct",
        action="store_true",
        help="Run cellranger count (if running individual pipeline components)"
    )

    parser.add_argument(
        "--cellbender",
        "-cb",
        action="store_true",
        help="Run cellbender (if running individual pipeline components)"
    )

    parser.add_argument(
        "--no-cellbender",
        "-nb",
        action="store_true",
        help="Do not run cellbender when running the entire pipeline"
    )

    parser.add_argument(
        "--spatial-count",
        "-sc",
        action="store_true",
        help="Run spatial barcode counts program (if running individual pipeline components)"
    )

    parser.add_argument(
        "--spatial-analysis",
        "-sa",
        action="store_true",
        help="Run spatial analysis (if running individual pipeline components)"
    )

    parser.add_argument(
        "--run-all",
        "-ra",
        action="store_true",
        help="Run all pipeline components from start to finish"
    )

    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Overwrite existing outputs when running modules"
    )

    parser.add_argument(
        "--gcp",
        "-gc",
        action="store_true",
        help="Run the pipeline on Google Cloud Platform"
    )

    parser.add_argument(
        "--version",
        "-v",
        action="store_true",
        help="Print out software version"
    )

    return parser.parse_args()


def _load() -> tuple[argparse.ArgumentParser, dict]:
    """
    load configfile parameters, build project structure, and begin runtime logging

    outputs:
     - an argparse.Namespace object containing command-line arguments
     - a dictionary containing configfile parameters and global constants

    """

    # start timer
    START_TIME = time.time()

    # parse command-line arguments
    args = _parse_arguments()
    if args.version:
        print(f"slidr version {get_version()}")
        sys.exit(0)

    # --bcl is required for an actual run but deliberately not marked required in the parser, so that
    # --version works on its own; enforce it here now that the informational flags have been handled
    if args.bcl is None:
        err_console.print("[bold red]\\[ERROR][/bold red]: --bcl is required")
        err_console.print("Troubleshooting:")
        err_console.print(" • Pass the BCL run ID, e.g. `--bcl 20240101_RUNID`")
        err_console.print(" • It must match a value in the `BCL` column of your metadata sheet")
        err_console.print(" • Give the ID itself, not a path - the run folder is resolved as <input_path>/<BCL_ID>")
        sys.exit(1)


    # render the flags this run was given back into the spelling a user could retype; argparse uses
    # underscores in flag variable names - these are replaced with hyphens to give the actual flag spelling
    provided_args = []
    for arg, val in vars(args).items():
        if arg == "bcl" or val is None or val is False:
            continue
        flag = f"--{arg.replace('_', '-')}"
        provided_args.append(flag if val is True else f"{flag} {val}")

    # set config/config.yaml as the default configuration file
    ROOT_DIR = Path(__file__).parent.parent.parent
    config_path = ROOT_DIR / 'config' / 'config.yaml' if not args.config else Path(args.config)
    if not config_path.is_file():
        err_console.print(f"[bold red]\\[ERROR][/bold red]: {config_path} is not a file!")
        err_console.print("Troubleshooting:")
        err_console.print(" • Make sure a config/ directory exists inside the project root folder")
        err_console.print(" • Make sure a config.yaml file exists inside the config/ directory")
        sys.exit(1)

    # set script location as a global variable
    SCRIPT_PATH = ROOT_DIR / "workflow" / "scripts"
    if not SCRIPT_PATH.is_dir():
        err_console.print("[bold red]\\[ERROR][/bold red]: `scripts` subdirectory not found in project root!")
        err_console.print("Troubleshooting:")
        err_console.print(" • Restore the missing `scripts` subdirectory by running `git pull`")
        sys.exit(1)

    # set the software cache file path
    SOFTWARE_CACHE_FILE = ROOT_DIR / "software_cache.txt"

    # load configuration parameters
    with open(config_path, 'r') as file:
        try:
            config = yaml.safe_load(file)
        except Exception as exc:
            err_console.print(f"[ERROR]: could not load the configuration file config/config.yaml: {exc}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Validate the config file syntax with a YAML validator: https://www.yamllint.com")
            err_console.print(" • Make sure you have the read permissions on the file")       
            sys.exit(1)

        # retrieve individual configfile sections
        paths = config.get('paths', {})
        settings = config.get('settings', {})
        workflow = config.get('workflow', {})

        # allow environment variables to override certain paths to let Slurm jobs
        # write into their own work directories
        applied_overrides = []
        for section, overrides in ((paths, PATH_ENV_OVERRIDES), (settings, SETTINGS_ENV_OVERRIDES)):
            for field, env_var in overrides.items():
                override = os.environ.get(env_var)
                if override is not None and override.strip():
                    section[field] = override.strip()
                    applied_overrides.append(f"{field} <- {env_var}")

        # set global parameter values
        OUT_PATH = paths.get('output_path')
        INPUT_PATH = paths.get('input_path')
        SOFTWARE_PATH = paths.get('software_path')
        RAW_BARCODES_PATH = paths.get('raw_barcodes_path')
        AUTH_KEY_PATH = paths.get('auth_key_path')
        PUCK_PATH = paths.get('puck_path')
        REF_PATH = paths.get('reference_path')

        NUM_THREADS = settings.get('threads')
        MEM_SIZE = settings.get('memory')
        METADATA_SRC = settings.get('metadata_source')
        GCS_DOWNLOAD_DEST = settings.get('gcs_download_dest')
        OUTPUT_BUCKET = settings.get('output_bucket')
        REF_GENOME = settings.get('reference_genome')
        ALERTS = settings.get('alerts')
        SLACK_TOKEN = settings.get('slack_token')

        GENERATE_BAM = workflow.get('generate_bam')
        CELLBENDER_EPOCHS = workflow.get('cellbender_epochs')
        CELLBENDER_RATE = workflow.get('cellbender_learn_rate')
        CELLBENDER_DROPLETS = workflow.get('cellbender_total_droplets')
        CELLBENDER_CELLS = workflow.get('cellbender_estimated_cells')
        SPATIAL_DOWNSAMPLING = workflow.get('spatial_downsampling')
        PERCENT_UMI_FILTERING = workflow.get('top_n_percent_umi_filter')
        MAX_BEADS_PER_CELL = workflow.get('max_beads_per_cell')
        EMPTYDROPS_MIN_UMI = workflow.get('flex_emptydrops_minimum_umis')
        FLEX_PROBE_SET = workflow.get('flex_probe_set')
        FLEX_R1_PATH = workflow.get('flex_spatial_R1_path')
        FLEX_R2_PATH = workflow.get('flex_spatial_R2_path')
        FLEX_GEX_FASTQS = workflow.get('flex_gex_fastqs')
        TREKKER_SAMPLESHEET = workflow.get('trekker_samplesheet')
        SEURAT_FORMAT = workflow.get('seurat_format')

    # set the output GCS bucket destination from an environment variable and validate it
    OUTPUT_DEST = os.environ.get('SLIDR_OUTPUT_DEST')
    OUTPUT_DEST = OUTPUT_DEST.strip().rstrip('/') if OUTPUT_DEST and OUTPUT_DEST.strip() else None
    if OUTPUT_DEST is not None:
        if not is_gcs_path(OUTPUT_DEST):
            err_console.print(f"[bold red]\\[ERROR][/bold red]: SLIDR_OUTPUT_DEST must be a gs:// URI naming the folder to upload results to, but is: {OUTPUT_DEST!r}")
            err_console.print("Troubleshooting:")
            err_console.print(" • This is normally set for you by ./slidr --gcp; unset it to let the run choose its own folder under `settings.output_bucket`")
            err_console.print(" • Write the full URL including the scheme, e.g. gs://my-bucket/20240101_RUNID")
            sys.exit(1)
        applied_overrides.append("output destination <- SLIDR_OUTPUT_DEST")

    # override metadata source explicitly with the `--metadata` argument
    metadata_override = None
    if args.metadata is not None and args.metadata.strip():
        metadata_override = f"metadata_source <- --metadata"
        METADATA_SRC = args.metadata.strip()

    # validate configfile values:
    if not isinstance(METADATA_SRC, str) or not METADATA_SRC.strip():
        err_console.print(f"[bold red]\\[ERROR][/bold red]: no metadata source set: `settings.metadata_source` is {METADATA_SRC!r} and --metadata was not given")
        err_console.print("Troubleshooting:")
        err_console.print(" • Set `settings.metadata_source` to a .tsv/.csv path or a Google Sheet URL")
        err_console.print(" • Or pass one for this run only with `--metadata /path/to/metadata.tsv`")
        sys.exit(1)

    # validate the provided authentication key file/path
    if is_gcs_path(AUTH_KEY_PATH):
        AUTH_KEY_PATH = AUTH_KEY_PATH.strip().rstrip('/')
        if not AUTH_KEY_PATH[len('gs://'):]:
            err_console.print(f"[bold red]\\[ERROR][/bold red]: `auth_key_path` is a gs:// URI with no bucket or object: {AUTH_KEY_PATH!r}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Give the full object, e.g. gs://my-bucket/secrets/auth_key.json")
            sys.exit(1)
        if not is_google_sheet(METADATA_SRC):
            # nothing will read it, so do not spend a download on it
            err_console.print(f"[bold yellow]\\[WARNING][/bold yellow]: `auth_key_path` names a GCS object, but this run reads its metadata from a local file")
            err_console.print(" • The key is not needed and will not be downloaded")
            AUTH_KEY_PATH = None

    elif is_google_sheet(METADATA_SRC):
        if AUTH_KEY_PATH is None or not Path(AUTH_KEY_PATH).is_file():
            err_console.print(f"[bold red]\\[ERROR][/bold red]: `auth_key_path` must point at a service-account JSON key to read metadata from a Google Sheet, but is: {AUTH_KEY_PATH}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Set `paths.auth_key_path` to a Google service-account JSON key file, or to a gs:// object holding one")
            err_console.print(" • Create one under IAM & Admin → Service Accounts → Keys in the Google Cloud console, then share the sheet with its client_email")
            err_console.print(" • `gcloud auth login` is not a substitute for this file -- it is read directly, not through Application Default Credentials")
            err_console.print(" • Or point `settings.metadata_source` (or --metadata) at an exported .tsv/.csv instead, which needs no key at all")
            sys.exit(1)
        AUTH_KEY_PATH = Path(AUTH_KEY_PATH)

    elif AUTH_KEY_PATH is not None and str(AUTH_KEY_PATH).strip():
        # not needed for this run, but validate a configured value rather than silently ignoring a wrong path
        if not Path(AUTH_KEY_PATH).is_file():
            err_console.print(f"[bold yellow]\\[WARNING][/bold yellow]: `auth_key_path` is set but is not a file: {AUTH_KEY_PATH}")
            err_console.print(" • This run reads its metadata from a local file, so the key is not needed and the run will continue")
            err_console.print(" • Fix or clear the field before switching `metadata_source` to a Google Sheet")
            AUTH_KEY_PATH = None
        else:
            AUTH_KEY_PATH = Path(AUTH_KEY_PATH)

    # validate the paths for pre-compiled pucks and raw puck files if running spatial barcode module
    if args.spatial_count or args.run_all:
        PUCK_PATH = _validate_resource_path(
            PUCK_PATH, 'puck_path', 'pre-compiled barcode/spatial coordinate puck maps',
            required=False)
        RAW_BARCODES_PATH = _validate_resource_path(
            RAW_BARCODES_PATH, 'raw_barcodes_path',
            'raw barcode/spatial coordinate files for puck map generation', required=False)

    # validate the reference path if running cellranger count
    if args.count or args.spatial_analysis or args.run_all:
        REF_PATH = _validate_resource_path(
            REF_PATH, 'reference_path', 'reference genome directories')

    # validate the software path and check if it is a GCS location
    if is_gcs_path(SOFTWARE_PATH):
        SOFTWARE_PATH = gcs_uri(SOFTWARE_PATH)

        if not any(is_gcs_path(other) for other in (INPUT_PATH, REF_PATH, PUCK_PATH, RAW_BARCODES_PATH)):
            console.print(f"[bold yellow]\\[NOTE][/bold yellow]: `software_path` is a GCS location, so software will be downloaded from {SOFTWARE_PATH} on first use")
            console.print("  • Nothing else is staged for this run; set `software_path` to a local directory to avoid the download")
            console.print(f"  • Already-known executables in {SOFTWARE_CACHE_FILE} are used first, and skip it entirely")
    elif SOFTWARE_PATH is None or not Path(SOFTWARE_PATH).is_dir():
        err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for `software_path` in configuration file: {SOFTWARE_PATH}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check that `software_path` is set to a valid directory path containing cellranger/bcl2fastq executables")
        err_console.print(" • If the software lives in a bucket, write it as a full gs:// URL -- it is then downloaded on first use")
        err_console.print(f" • Or pin the executables directly by adding their paths to {SOFTWARE_CACHE_FILE}, which is consulted before this directory is scanned")
        sys.exit(1)
    else:
        SOFTWARE_PATH = Path(SOFTWARE_PATH)

    # validate allocated memory size value
    if not is_int(MEM_SIZE):
        err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for `memory` in configuration file: {MEM_SIZE!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check that `memory` is set to a integer value corresponding to the maximum amount of memory to be used by the pipeline (in GB)")
        hint = bool_value_hint(MEM_SIZE, 'settings.memory')
        if hint:
            err_console.print(hint)
        sys.exit(1)

    # validate allocated number of threads
    if not is_int(NUM_THREADS):
        err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for `threads` in configuration file: {NUM_THREADS!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check that `threads` is set to a integer value corresponding to the maximum number of CPU cores to be used by the pipeline")
        hint = bool_value_hint(NUM_THREADS, 'settings.threads')
        if hint:
            err_console.print(hint)
        sys.exit(1)

    # handle an ambiguous case where both the preprocessed FASTQs and raw BCLs are provided
    if args.fastqs is not None and args.mkfastq:
        err_console.print("[bold red]\\[ERROR][/bold red]: --fastqs and --mkfastq are mutually exclusive")
        err_console.print("Troubleshooting:")
        err_console.print(" • Drop --fastqs to demultiplex BCLs from `input_path` with cellranger mkfastq")
        err_console.print(" • Drop --mkfastq to run the downstream stages on the FASTQs you provided")
        sys.exit(1)

    # validate pipeline alerts settings
    if isinstance(ALERTS, str):
        if ALERTS.lower() in ['y', 'yes', 't', 'true']:
            ALERTS = True
        elif ALERTS.lower() in ['n', 'no', 'f', 'false']:
            ALERTS = False

    if ALERTS is None:
        ALERTS = False
    
    if not isinstance(ALERTS, bool):
        err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for the `alerts` field in the configfile: {ALERTS}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check that `alerts` is set to `yes`/`no` or `true`/`false`")
        sys.exit(1)

    # resolve the provided Slack token to a valid string
    if ALERTS:
        SLACK_TOKEN = _resolve_slack_token(SLACK_TOKEN)

    # validate BAM generation settings
    if isinstance(GENERATE_BAM, str):
        if GENERATE_BAM.lower() in ['f', 'false', 'n', 'no']:
            GENERATE_BAM = False
        elif GENERATE_BAM.lower() in ['t', 'true', 'y', 'yes']:
            GENERATE_BAM = True

    if GENERATE_BAM is None:
        GENERATE_BAM = False

    if not isinstance(GENERATE_BAM, bool):
        err_console.print(f'[ERROR]: unrecognized value for `generate_bam` field in configuration file (should be `true` or `false`): {GENERATE_BAM}')
        err_console.print(f"Troubleshooting:")
        err_console.print(f" • Edit the config file and set the 'generate_bam' field in the `workflow` section to")
        err_console.print(f"    - 'true' if you want cellranger count/multi to generate BAM genome alignment files OR")
        err_console.print(f"    - 'false' if you don't BAM files to be generated")
        sys.exit(1)

    # define metadata sheet specifications
    METADATA_SPECS = [
        'Run',
        'Email',
        'Sample Name',
        'BCL',
        'Species',
        'Chemistry',
        'RNA Index',
        'Lane',
        'SB Index',
        'SB Lane',
        'Puck ID'
    ]

    # if no output path is provided, the pipeline creates an 'outs' directory in the project root
    if OUT_PATH is None:
        OUT_PATH = ROOT_DIR / "outs"
        OUT_PATH.mkdir(exist_ok=True)

    # append the BCL ID to the output path, unless it already names the run directory
    BCL_ID = args.bcl
    OUT_PATH = Path(OUT_PATH)
    if OUT_PATH.name != BCL_ID:
        OUT_PATH = OUT_PATH / BCL_ID
    OUT_PATH.mkdir(exist_ok=True, parents=True)

    # define and create output subdirectories
    LOG_PATH = OUT_PATH / "log"
    METADATA_PATH = OUT_PATH / "metadata"
    OUTPUT_PATH = OUT_PATH / "output"
    TMP_PATH = OUT_PATH / "tmp"

    # on a CIFS/SMB run directory the logs are written to local disk for the duration of the run and
    # copied into <OUT_PATH>/log at exit, because slidr's per-line open/append/close cannot be done
    # reliably over SMB
    LOG_PATH, LOG_DEST = resolve_log_staging(LOG_PATH, BCL_ID)

    # define the remaining output files/subdirectories
    RUNTIME_LOG = LOG_PATH / "runtime.log"
    SUMMARY_PATH = METADATA_PATH / "metadata_summary.csv"
    SAMPLESHEET_PATH = METADATA_PATH / "samplesheet.csv"
    MKFASTQ_OUTS = OUTPUT_PATH / "mkfastq"
    COUNT_OUTS = OUTPUT_PATH / "count"
    CELLBENDER_OUTS = OUTPUT_PATH / "cellbender"
    SPATIAL_COUNT_OUTS = OUTPUT_PATH / "spatial_barcodes"
    SPATIAL_ANALYSIS_OUTS = OUTPUT_PATH / "spatial_analysis"
    FLEX_OUTS = OUTPUT_PATH / "flex"

    # create output directories
    for path in [LOG_PATH, LOG_DEST, METADATA_PATH, TMP_PATH, OUTPUT_PATH]:
        if path is not None:
            path.mkdir(exist_ok=True)

    # resolve where this run reads its input from
    INPUT_BUCKET = None
    STAGE_FASTQ_INPUT = False
    FASTQ_INPUT = args.fastqs
    input_is_bucket = is_gcs_path(INPUT_PATH)

    # validate the destination for GCS downloads if input is a Google Bucket
    if input_is_bucket:
        INPUT_BUCKET = gcs_uri(INPUT_PATH)
        if GCS_DOWNLOAD_DEST is None or not str(GCS_DOWNLOAD_DEST).strip():
            GCS_DOWNLOAD_DEST = OUT_PATH / "data"
        else:
            GCS_DOWNLOAD_DEST = Path(str(GCS_DOWNLOAD_DEST).strip())
        try:
            GCS_DOWNLOAD_DEST.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            err_console.print(f"[bold red]\\[ERROR][/bold red]: could not create the staging directory `gcs_download_dest`: {GCS_DOWNLOAD_DEST} ({exc})")
            err_console.print("Troubleshooting:")
            err_console.print(" • Check the parent directory exists and is writable from this machine")
            err_console.print(" • Set `settings.gcs_download_dest` to somewhere with room for the sequencing data (a BCL run is routinely hundreds of GB)")
            err_console.print(f" • Left unset it defaults to {OUT_PATH / 'data'}, inside the run's own directory")
            sys.exit(1)
        INPUT_PATH = GCS_DOWNLOAD_DEST

        # define the FASTQ input path for cases where data are staged from a Google Bucket
        if FASTQ_INPUT == FASTQS_USE_INPUT_PATH:
            FASTQ_INPUT = INPUT_PATH / BCL_ID
            STAGE_FASTQ_INPUT = True

    else:
        # a local input_path stages nothing
        GCS_DOWNLOAD_DEST = None

    # validate the FASTQ directory supplied with --fastqs
    if FASTQ_INPUT is not None and not STAGE_FASTQ_INPUT:
        if FASTQ_INPUT == FASTQS_USE_INPUT_PATH:
            if INPUT_PATH is None or not str(INPUT_PATH).strip():
                err_console.print("[bold red]\\[ERROR][/bold red]: --fastqs was given without a directory, so `paths.input_path` must point at the FASTQs -- but it is not set")
                err_console.print("Troubleshooting:")
                err_console.print(" • Set `paths.input_path` to the directory holding the .fastq.gz files")
                err_console.print(" • Or name the directory on the command line instead: `--fastqs /path/to/fastqs`")
                sys.exit(1)
            FASTQ_INPUT = INPUT_PATH
            source = "`paths.input_path`"
        else:
            source = "--fastqs"

        FASTQ_INPUT = Path(FASTQ_INPUT)
        if not FASTQ_INPUT.is_dir():
            err_console.print(f"[bold red]\\[ERROR][/bold red]: the FASTQ directory from {source} is not a directory or does not exist: {FASTQ_INPUT}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Check the path for typos and that it is readable")
            err_console.print(" • --fastqs expects a directory of .fastq.gz files, not an individual FASTQ file")
            sys.exit(1)
        if not any(FASTQ_INPUT.rglob('*.fastq.gz')):
            err_console.print(f"[bold red]\\[ERROR][/bold red]: no .fastq.gz files found in the FASTQ directory from {source}: {FASTQ_INPUT}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Check that the FASTQs are gzipped (cellranger requires .fastq.gz, not plain .fastq)")
            err_console.print(" • Omit --fastqs to demultiplex BCLs from `input_path` instead")
            sys.exit(1)

    # validate the data input path
    if INPUT_PATH is None or not Path(INPUT_PATH).is_dir():
        if FASTQ_INPUT is not None:
            INPUT_PATH = Path(FASTQ_INPUT)
        else:
            err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for the `input_path` field in the configfile: {INPUT_PATH}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Check that `input_path` is set to a valid directory path containing input BCL run folders")
            err_console.print(" • Alternatively, pass a directory of already-demultiplexed FASTQs with --fastqs")
            err_console.print(" • If the run folders live in a bucket, set `input_path` to that gs:// prefix")
            sys.exit(1)
    else:
        INPUT_PATH = Path(INPUT_PATH)

    # find current git branch name
    try:
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            text=True,
        ).strip()
        branch = f"slidr/{branch}"
    except:
        branch = "[unknown]"

    # set up summary file
    SUMMARY_LOG = LOG_PATH / f"{datetime.now().strftime('%y-%m-%d_%H:%M')}_summary.log"
    with open(SUMMARY_LOG, "w") as summary:
        summary.write("===================== [ RUN SUMMARY ] ======================\n\n")
        summary.write(f"  Date:                 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        summary.write(f"  BCL ID:               {BCL_ID}\n")
        summary.write(f"  Input path:           {INPUT_PATH}\n")
        if FASTQ_INPUT is not None:
            summary.write(f"  FASTQ input (--fastqs): {FASTQ_INPUT}\n")
        summary.write(f"  Output path:          {OUT_PATH}\n")
        summary.write(f"  Output bucket:        {OUTPUT_BUCKET}\n")
        if OUTPUT_DEST is not None:
            summary.write(f"  Upload to:            {OUTPUT_DEST}\n")
        if INPUT_BUCKET is not None:
            summary.write(f"  Input (GCS):          {INPUT_BUCKET}\n")
            summary.write(f"  Staged into:          {GCS_DOWNLOAD_DEST}\n")
        if is_gcs_path(REF_PATH):
            summary.write(f"  Reference (GCS):      {REF_PATH}\n")
        if is_gcs_path(PUCK_PATH):
            summary.write(f"  Pucks (GCS):          {PUCK_PATH}\n")
        if is_gcs_path(RAW_BARCODES_PATH):
            summary.write(f"  Raw barcodes (GCS):   {RAW_BARCODES_PATH}\n")
        if is_gcs_path(SOFTWARE_PATH):
            summary.write(f"  Software (GCS):       {SOFTWARE_PATH}\n")
        if applied_overrides:
            summary.write(f"  Config overrides:     {', '.join(applied_overrides)}\n")
        if metadata_override:
            summary.write(f"  Metadata override:    {metadata_override}\n")
        summary.write(f"  Metadata source:      {METADATA_SRC}\n")
        summary.write(f"  Pipeline version:     {get_version()}\n")
        summary.write(f"  GitHub branch:        {branch}\n")
        summary.write(f"  Memory:               {MEM_SIZE} GB\n")
        summary.write(f"  Cores:                {NUM_THREADS}\n")
        summary.write(f"  Arguments:\n")
        for arg in provided_args:
            summary.write(f"   {arg}\n")

    # display the stages requested by this invocation
    STAGE_FLAGS = [
        ('mkfastq', args.mkfastq),
        ('count', args.count),
        ('cellbender', args.cellbender and not args.no_cellbender),
        ('spatial-count', args.spatial_count),
        ('spatial-analysis', args.spatial_analysis),
    ]
    if args.run_all:
        REQUESTED_STAGES = [name for name, _ in STAGE_FLAGS if not (name == 'cellbender' and args.no_cellbender)]
    else:
        REQUESTED_STAGES = [name for name, requested in STAGE_FLAGS if requested]
    if args.fastqs is not None and 'mkfastq' in REQUESTED_STAGES:
        # providing FASTQs as input obviates the need to run mkfastq
        REQUESTED_STAGES.remove('mkfastq')
    SKIPPED_STAGES = [name for name, _ in STAGE_FLAGS if name not in REQUESTED_STAGES]

    # set up runtime logfile
    RUN_WIDTH = 74

    # auxiliary formatting functions
    def _banner(text: str, fill: str = '═') -> str:
        pad = max(0, RUN_WIDTH - len(text) - 2)
        left = pad // 2
        return f"{fill * left} {text} {fill * (pad - left)}"

    def _field(label: str, value) -> str:
        return f"  {label:<14} {value}"

    # write a run setup report to the logfile
    with open(RUNTIME_LOG, "w") as logfile:
        logfile.write(_banner(f"slidr {get_version()}") + "\n")
        logfile.write(_field("Run", f"{BCL_ID:<31} {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}") + "\n")
        logfile.write(_field("Output", OUT_PATH) + "\n")
        logfile.write(_field("Input", FASTQ_INPUT if FASTQ_INPUT is not None else INPUT_PATH) + "\n")
        if FASTQ_INPUT is not None:
            logfile.write(_field("", "(already demultiplexed FASTQs, via --fastqs)") + "\n")
        source = "Google Sheet" if is_google_sheet(METADATA_SRC) else METADATA_SRC
        logfile.write(_field("Metadata", source) + "\n")
        logfile.write(_field("Resources", f"{NUM_THREADS} cores · {MEM_SIZE} GB") + "\n")
        if INPUT_BUCKET is not None:
            logfile.write(_field("Staging", f"{INPUT_BUCKET} -> {GCS_DOWNLOAD_DEST}") + "\n")
        logfile.write(_field("Stages", ', '.join(REQUESTED_STAGES) or "none requested") + "\n")
        if SKIPPED_STAGES:
            logfile.write(_field("", f"(not requested: {', '.join(SKIPPED_STAGES)})") + "\n")
        if provided_args:
            logfile.write(_field("Arguments", ' '.join(provided_args)) + "\n")
        logfile.write("─" * RUN_WIDTH + "\n")

    console.print('[bold green]\\[SUCCESS][/bold green]: Setup complete')

    # optionally provide a Trekker samplesheet for running Flex analyses rather than re-deriving it
    if TREKKER_SAMPLESHEET is not None and str(TREKKER_SAMPLESHEET).strip():
        TREKKER_SAMPLESHEET = Path(str(TREKKER_SAMPLESHEET).strip()).expanduser()
        if not TREKKER_SAMPLESHEET.is_file():
            err_console.print(f"[bold red]\\[ERROR][/bold red]: `workflow.trekker_samplesheet` is set but is not a file: {TREKKER_SAMPLESHEET}")
            err_console.print("Troubleshooting:")
            err_console.print(" • Point it at a Trekker pipeline samplesheet (.csv or .tsv), or remove the field to have slidr generate one")
            err_console.print(" • The path is read on the machine the pipeline runs on, so a --gcp run needs it present on the VM")
            sys.exit(1)
    else:
        TREKKER_SAMPLESHEET = None

    # set the serialization format of the spatial-analysis Seurat object;
    # (qs2 is the default and replaces the deprecated qs format)
    SEURAT_FORMATS = ('qs2', 'qs')
    if SEURAT_FORMAT is None or (isinstance(SEURAT_FORMAT, str) and not SEURAT_FORMAT.strip()):
        SEURAT_FORMAT = 'qs2'
    elif not isinstance(SEURAT_FORMAT, str) or SEURAT_FORMAT.strip().lower() not in SEURAT_FORMATS:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: `workflow.seurat_format` must be one of {', '.join(SEURAT_FORMATS)}, but is {SEURAT_FORMAT!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • `qs2` writes seurat.qs2, read back with qs2::qs_read() -- the default, and what to use for new work")
        err_console.print(" • `qs` writes seurat.qs in the deprecated format, read back with qs::qread() -- for downstream code that still expects it")
        err_console.print(" • Remove the field entirely to accept the default of qs2")
        sys.exit(1)
    else:
        SEURAT_FORMAT = SEURAT_FORMAT.strip().lower()

    # package global constants into a dictionary
    cfg = {
        'root_path':                    ROOT_DIR,
        'run_path':                     OUT_PATH,
        'requested_stages':             REQUESTED_STAGES,
        'output_path':                  OUTPUT_PATH,
        'input_path':                   INPUT_PATH,
        'fastq_input':                  FASTQ_INPUT,
        'stage_fastq_input':            STAGE_FASTQ_INPUT,
        'script_path':                  SCRIPT_PATH,
        'metadata_path':                METADATA_PATH,
        'log_path':                     LOG_PATH,
        'log_dest':                     LOG_DEST,
        'tmp_path':                     TMP_PATH,
        'software_path':                SOFTWARE_PATH,
        'metadata_src':                 METADATA_SRC,
        'raw_barcodes_path':            RAW_BARCODES_PATH,
        'puck_path':                    PUCK_PATH,
        'ref_path':                     REF_PATH,
        'ref_genome':                   REF_GENOME,
        'software_cache_file':          SOFTWARE_CACHE_FILE,
        'auth_key_path':                AUTH_KEY_PATH,
        'input_bucket':                 INPUT_BUCKET,
        'gcs_download_dest':            GCS_DOWNLOAD_DEST,
        'output_bucket':                OUTPUT_BUCKET,
        'output_dest':                  OUTPUT_DEST,
        'num_threads':                  NUM_THREADS,
        'mem_size':                     MEM_SIZE,
        'alerts':                       ALERTS,
        'slack_token':                  SLACK_TOKEN,
        'metadata_specs':               METADATA_SPECS,
        'runtime_log':                  RUNTIME_LOG,
        'summary_log':                  SUMMARY_LOG,
        'summary_path':                 SUMMARY_PATH,
        'samplesheet_path':             SAMPLESHEET_PATH,
        'mkfastq_outs':                 MKFASTQ_OUTS,
        'count_outs':                   COUNT_OUTS,
        'cellbender_outs':              CELLBENDER_OUTS,
        'spatial_count_outs':           SPATIAL_COUNT_OUTS,
        'spatial_analysis_outs':        SPATIAL_ANALYSIS_OUTS,
        'flex_outs':                    FLEX_OUTS,
        'generate_bam':                 GENERATE_BAM,
        'cellbender_epochs':            CELLBENDER_EPOCHS,
        'cellbender_rate':              CELLBENDER_RATE,
        'cellbender_cells':             CELLBENDER_CELLS,
        'cellbender_droplets':          CELLBENDER_DROPLETS,
        'spatial_downsampling':         SPATIAL_DOWNSAMPLING,
        'percent_umi_filtering':        PERCENT_UMI_FILTERING,
        'max_beads_per_cell':           MAX_BEADS_PER_CELL,
        'emptydrops_min_umis':          EMPTYDROPS_MIN_UMI,
        'flex_probe_set':               FLEX_PROBE_SET,
        'flex_r1_path':                 FLEX_R1_PATH,
        'flex_r2_path':                 FLEX_R2_PATH,
        'flex_gex_fastqs':              FLEX_GEX_FASTQS,
        'trekker_samplesheet':          TREKKER_SAMPLESHEET,
        'seurat_format':                SEURAT_FORMAT,
        'start_time':                   START_TIME,
        'bcl_id':                       BCL_ID
    }

    return args, cfg


args, cfg = _load()