import yaml
import os
import sys
import time
import argparse
import tomllib
import subprocess

from datetime import datetime
from pathlib import Path
from rich.console import Console
from importlib.metadata import version, PackageNotFoundError


# initialize rich consoles for text formatting
err_console = Console(stderr=True)
console = Console()


# Timeouts for the `gcloud` subprocesses the pipeline shells out to. `subprocess.run` waits forever
# by default, which is how one wedged transfer could hang a whole run until the scheduler's wall
# clock killed it -- a failure mode that cost a 24h Slurm allocation with the analysis already
# finished and only the upload left to do.
#
# Two tiers, because the two kinds of call have honest durations orders of magnitude apart:
#  - metadata: `gcloud storage cat`/`ls` against a single small object. These take seconds, so a few
#    minutes is already pathological; failing fast turns a silent stall into a real error message.
#  - transfer: `gcloud storage cp` of a reference genome, a software tree or a whole BCL run folder,
#    routinely hundreds of GB. Hours here are legitimate, so this is deliberately generous. It is a
#    backstop against a process that will never finish, NOT a performance budget -- set it too tight
#    and it kills transfers that were going to succeed.
#
# Both are overridable, since link speed and run size vary enormously between sites.
def _timeout_from_env(name: str, default: int) -> int:
    """
    Read a positive-integer timeout (in seconds) from the environment, falling back to `default`.

    A malformed or non-positive value falls back rather than failing the run: these are backstops,
    and refusing to start over a typo in an optional tuning knob would be worse than ignoring it.
    """

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


GCS_METADATA_TIMEOUT = _timeout_from_env('SLIDR_GCS_METADATA_TIMEOUT', 5 * 60)
GCS_TRANSFER_TIMEOUT = _timeout_from_env('SLIDR_GCS_TRANSFER_TIMEOUT', 12 * 60 * 60)

# The same reasoning for the Julia toolchain slidr installs for itself when the software directory
# has no julia. Both steps reach the network and neither had a limit, so a dead proxy or a black-holed
# route left the run wedged before spatial counting rather than failing with something actionable.
# Split for the same reason as the GCS pair: fetching the installer is a few KB, while running it
# downloads and unpacks a whole toolchain.
INSTALLER_FETCH_TIMEOUT = _timeout_from_env('SLIDR_INSTALLER_FETCH_TIMEOUT', 5 * 60)
INSTALLER_RUN_TIMEOUT = _timeout_from_env('SLIDR_INSTALLER_RUN_TIMEOUT', 30 * 60)


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


def is_gcs_path(value) -> bool:
    """
    Whether a configured value explicitly names a Google Cloud Storage object.

    Only a literal `gs://` prefix counts. The `paths.*` directory fields also accept a bare
    `bucket/prefix` (normalized by gcs_uri()), but that form is indistinguishable from a relative
    local path, which is why that form is not accepted for those fields at all. A `gs://` URI needs no such
    gate: it cannot be mistaken for a local path, so it is safe to act on directly.

    Inputs:
     - value: the configured value to check
    Output:
     - True if the value is a gs:// URI
    """

    return isinstance(value, str) and value.strip().startswith('gs://')


def is_google_sheet(source) -> bool:
    """
    Whether a metadata source is a Google Sheet URL rather than a local file.

    Kept here so config.py's validation and helpers.load_metadata() agree on what counts as a Sheet:
    the service-account key is required for exactly the sources this returns True for, and demanding
    it for a local .tsv/.csv would force users who never touch Sheets to produce a key anyway.

    Inputs:
     - source: the configured or overridden metadata source
    Output:
     - True if the source is a Google Sheets URL
    """

    return isinstance(source, str) and source.strip().startswith("https://docs.google.com/spreadsheets")


def is_int(value) -> bool:
    """
    Whether a configuration value is a usable integer.

    `bool` is a subclass of `int`, so a bare isinstance(value, int) accepts True/False -- and PyYAML
    reads the unquoted words true/false/yes/no/on/off as booleans, so `threads: no` is an easy typo
    to make. Left unguarded those arrive as the numbers 1 and 0 and are used as real settings
    (a thread count of 1, a droplet cutoff of 0) instead of being rejected. Every integer config
    field goes through here so that rule lives in one place.

    Note that the equivalent float-only checks need no such guard: bool does NOT subclass float,
    so isinstance(True, float) is already False.

    Inputs:
     - value: the configuration value to check
    Output:
     - True if the value is an int and not a bool
    """

    return isinstance(value, int) and not isinstance(value, bool)


def is_number(value) -> bool:
    """
    Whether a configuration value is a usable number (int or float), rejecting booleans as is_int
    does -- for fields that legitimately accept either, such as a percentage.

    Inputs:
     - value: the configuration value to check
    Output:
     - True if the value is an int or float and not a bool
    """

    return isinstance(value, (int, float)) and not isinstance(value, bool)


def bool_value_hint(value, field: str) -> str | None:
    """
    Extra troubleshooting line for a numeric field that was given a boolean, or None if the value is
    not a boolean. YAML's bare-word booleans are the single most likely way to land on one, so the
    message names that rather than just restating the expected type.

    Inputs:
     - value: the offending configuration value
     - field: dotted config field name, for the message
    Output:
     - a troubleshooting bullet, or None when `value` is not a boolean
    """

    if not isinstance(value, bool):
        return None
    return (f" • `{field}` takes a number, not a yes/no switch -- YAML reads the bare words"
            " true/false/yes/no/on/off as booleans, so quote the value if you meant a number")


def gcs_uri(path: Path | str) -> str:
    """
    Normalize a configured Google Cloud Storage location to a `gs://` URI.

    Config fields may be written either as a full `gs://bucket/prefix` URI or as a bare
    `bucket/prefix`; both are accepted and normalized here so callers can join onto the result
    without worrying which form the user wrote. Trailing slashes are stripped so
    f'{gcs_uri(x)}/{name}' never produces a double slash.

    Inputs:
     - path: configured GCS location, with or without the scheme
    Output:
     - the location as a `gs://`-prefixed URI with no trailing slash
    """

    text = str(path).strip().rstrip('/')
    if not text.startswith('gs://'):
        text = 'gs://' + text.lstrip('/')
    return text


# Machine-local paths that an environment variable may override, taking precedence over the config
# file. These exist so one config.yaml can be shared across machines: the remote-run setup scripts
# (workflow/bash/slidr_gcp.sh, workflow/bash/slidr_slurm.sh) stage inputs into a working directory
# whose location they choose at run time, and cannot know it when the config was written. Only
# genuinely host-specific locations are overridable -- everything else stays declarative in the
# config file.
# Sentinel stored by `--fastqs` when it is passed without a directory, meaning "the FASTQs are in
# `paths.input_path`". A sentinel rather than None so "flag absent" stays distinguishable from "flag
# present with no value" -- the two mean opposite things (demultiplex BCLs vs. skip mkfastq).
FASTQS_USE_INPUT_PATH = '<input_path>'


PATH_ENV_OVERRIDES = {
    'input_path': 'SLIDR_INPUT_PATH',
    'output_path': 'SLIDR_OUTPUT_PATH',
    'auth_key_path': 'SLIDR_AUTH_KEY_PATH',
    'software_path': 'SLIDR_SOFTWARE_PATH',
}

# The same mechanism for the one `settings.*` field that names a machine-local location. It is in
# `settings` rather than `paths` because it is not an input the user points at -- it is where staged
# inputs are put -- but a cluster still needs to be able to redirect it per host without editing the
# shared config.
SETTINGS_ENV_OVERRIDES = {
    'gcs_download_dest': 'SLIDR_GCS_DOWNLOAD_DEST',
}


def _first_token_line(text: str) -> str:
    """
    Extract a secret from the contents of a token file: the first line that is neither blank nor a
    `#` comment. Shared by the local-file and GCS readers so both tolerate the same file shape.

    Inputs:
     - text: full contents of the token file
    Output:
     - the token, or '' if the contents hold none
    """

    lines = (line.strip() for line in text.splitlines())
    return next((line for line in lines if line and not line.startswith('#')), '')


def _read_gcs_token(uri: str) -> str:
    """
    Read a secret out of a GCS object with `gcloud storage cat`, without writing it to disk.

    Errors deliberately quote gcloud's stderr but never the object's contents: this function only ever
    handles secrets, and a failure message is the one place a token could accidentally end up in a log
    that gets shared.

    Inputs:
     - uri: gs:// URI of the object holding the token
    Output:
     - the token
    """

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
    Resolve the configured `slack_token` into an actual Slack bot token.

    The field accepts either a path to a file containing the token, or the token itself. The file form
    is preferred and documented as such: a config.yaml gets copied between machines, staged to a
    bucket for --gcp runs and committed by accident, and a token embedded in it travels along with it.
    Keeping the secret in a separate, permission-restricted file means sharing the config shares no
    credential.

    A literal token is still accepted so existing configs keep working. The two are told apart by
    whether the value names an existing file, not by pattern-matching the token format -- Slack has
    changed its prefixes before, and a path is unambiguous.

    Inputs:
     - value: the raw `settings.slack_token` value
    Output:
     - the token, stripped of surrounding whitespace
    """

    if not isinstance(value, str) or not value.strip():
        err_console.print(f"[bold red]\\[ERROR][/bold red]: `alerts` is enabled but `slack_token` is not set: {value!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Point `slack_token` at a file containing the bot token (preferred, so the token stays out of the config file)")
        err_console.print(" • Or set it to the token itself")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    value = value.strip()

    # A gs:// value is read straight out of the bucket. Unlike the service-account key -- which the
    # Google auth library insists on loading from a path, so it has to be written to disk first -- the
    # Slack token is only ever handed to WebClient as a string, so it never needs to touch the
    # filesystem at all. Streaming it keeps the secret out of the output tree entirely.
    if is_gcs_path(value):
        return _read_gcs_token(value)

    candidate = Path(value).expanduser()

    # A value that is not an existing file is taken as the token itself. Guard on the *shape* of the
    # value first: something that looks like a path but is not there is far more likely to be a typo
    # in a filename than a literal token, and treating it as a token would surface much later as an
    # opaque Slack auth failure.
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

    if not token:
        err_console.print(f"[bold red]\\[ERROR][/bold red]: the Slack token file {candidate} contains no token")
        err_console.print("Troubleshooting:")
        err_console.print(" • The file should hold the bot token on its own line (lines starting with `#` are ignored)")
        err_console.print(" • Or disable Slack alerts by setting `alerts` to `false`/`no`")
        sys.exit(1)

    # The point of the file is to keep the secret narrowly readable, so say something when it is not.
    # A warning rather than an error: the permissions may be owned by a shared-filesystem policy the
    # user cannot change, and refusing to run over it would be unhelpful.
    mode = candidate.stat().st_mode
    if mode & 0o077:
        err_console.print(f"[bold yellow]\\[WARNING][/bold yellow]: the Slack token file {candidate} is readable by other users (mode {mode & 0o777:03o})")
        err_console.print(f" • Restrict it with `chmod 600 {candidate}`")

    return token


def _parse_arguments() -> argparse.ArgumentParser:
    """
    Parse command-line arguments

    Outputs:
     - an argparse.Namespace object containing command-line arguments

    """

    parser = argparse.ArgumentParser()

    # Not required=True: argparse enforces that before any argument is inspected, which would make
    # `--version` unusable on its own (it would abort with "the following arguments are required:
    # --bcl", forcing callers to pass a meaningless --bcl just to read the version). Presence is
    # checked in _load() instead, immediately after the --version short-circuit.
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

    # The path is optional: a bare --fastqs means "`paths.input_path` holds the FASTQs", which saves
    # naming the same directory twice for the common case. nargs='?' with a const distinguishes the
    # three states that matter -- flag absent (None), flag with no value (FASTQS_USE_INPUT_PATH), and
    # flag with a path.
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
    Load configfile parameters, build project structure, and begin runtime logging

    Outputs:
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
        err_console.print(" • Give the ID itself, not a path -- the run folder is resolved as <input_path>/<BCL_ID>")
        sys.exit(1)


    # Render the flags this run was given back into the spelling a user could retype. vars(args) is
    # keyed by argparse *dest* names, which use underscores -- logging those verbatim recorded
    # "--no_cellbender" and "--spatial_count", neither of which is a real flag, so a logged command
    # line could not be copied. Values are included for the flags that take one.
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

        paths = config.get('paths', {})
        settings = config.get('settings', {})
        workflow = config.get('workflow', {})

        # apply host-specific path overrides from the environment before anything reads `paths`, so
        # a config.yaml shared from a bucket can be pointed at this machine's staging directories
        applied_overrides = []
        for section, overrides in ((paths, PATH_ENV_OVERRIDES), (settings, SETTINGS_ENV_OVERRIDES)):
            for field, env_var in overrides.items():
                override = os.environ.get(env_var)
                if override is not None and override.strip():
                    section[field] = override.strip()
                    applied_overrides.append(f"{field} <- {env_var}")

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
        EMPTYDROPS_MIN_UMI = workflow.get('flex_emptydrops_minimum_umis')
        FLEX_PROBE_SET = workflow.get('flex_probe_set')
        FLEX_R1_PATH = workflow.get('flex_spatial_R1_path')
        FLEX_R2_PATH = workflow.get('flex_spatial_R2_path')
        FLEX_GEX_FASTQS = workflow.get('flex_gex_fastqs')

    # The exact bucket folder this run's results go to, when the caller has already chosen one.
    #
    # `settings.output_bucket` says which bucket; the folder inside it is normally decided at the end of
    # the run by pipeline.unique_gcs_dest, which numbers around a name that is already taken. A --gcp run
    # cannot work that way: ./slidr has to resolve the folder *before* creating the VM, because that is
    # where it uploads the config.yaml the VM boots from. So it resolves it there and passes it down
    # here, and the upload uses it verbatim rather than re-deciding and landing somewhere else.
    #
    # Environment-only and deliberately not a config field: it names one specific run's destination, so
    # a value written into a config file would send every later run using that file to the same folder --
    # the exact collision the numbering exists to prevent.
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

    # --metadata overrides settings.metadata_source, so a one-off sheet or an exported TSV can be
    # used without editing the config file (and, on a --gcp/--slurm run, without re-staging it).
    # Applied before validation so the override is checked exactly like a configured value.
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

    # The service-account key is only ever read to authenticate against the Google Sheets and Drive
    # APIs (helpers.load_metadata), so it is required for a Sheet-backed run and irrelevant to one
    # driven by a local .tsv/.csv. Requiring it unconditionally forced users who never touch Sheets
    # to produce a key -- and, on a staged --gcp/--slurm run, to upload one -- purely to pass this
    # check. Note `gcloud auth login` is not a substitute: the key file is read directly rather than
    # through Application Default Credentials.
    # A gs:// key is left as a string here and downloaded on demand by helpers.load_metadata: only its
    # form can be checked at config time, since verifying the object exists would mean a network round
    # trip during startup. It follows the same per-value rule as the `paths.*`
    # directory fields -- a gs:// URI is unambiguous where their bare `bucket/prefix` form is not, and
    # fields: the scheme decides, and nothing about one field constrains another,
    # breaking a run whose only remote input is the key.
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
        # not needed for this run, but validate a configured value rather than silently ignoring a
        # path that is wrong -- the next Sheet-backed run would be the one to discover it
        if not Path(AUTH_KEY_PATH).is_file():
            err_console.print(f"[bold yellow]\\[WARNING][/bold yellow]: `auth_key_path` is set but is not a file: {AUTH_KEY_PATH}")
            err_console.print(" • This run reads its metadata from a local file, so the key is not needed and the run will continue")
            err_console.print(" • Fix or clear the field before switching `metadata_source` to a Google Sheet")
            AUTH_KEY_PATH = None
        else:
            AUTH_KEY_PATH = Path(AUTH_KEY_PATH)
    
    # Whether the reference genome, puck maps and raw barcodes are read off a local filesystem or
    # staged out of GCS is stated explicitly, never inferred from the shape of the configured path.
    # --gcp always stages (a fresh GCE VM has no access to the lab filesystem by definition); elsewhere
    # requests the same behaviour anywhere else, which is what lets a Slurm job on an unrelated
    # cluster pull its own inputs.

    # `input_path` is stageable too, and is resolved further down -- once the output tree is known,
    # since that is where a staged run downloads to by default. See "resolve where this run reads its
    # reads from" below.

    def _validate_resource_path(value, field: str, contents: str, required: bool = True):
        """
        Validate one of the stageable resource paths, returning it normalized: a `gs://` URI string
        when it names a bucket, a local `Path` otherwise. Returns None for an unset optional field.

        One rule, read off the value itself: a `gs://` URI names a bucket, anything else is a local
        directory. No flag selects between them, and the fields are independent, so sequencing data in
        a bucket with the reference genome on local disk is an ordinary configuration.

        The bare `bucket/prefix` form is deliberately not accepted. It cannot be told apart from a
        relative directory, and the flag that used to disambiguate it bought that one
        capability at the price of a mode: every path had to agree about being remote, mixed setups
        were unexpressible, and the same value meant different things depending on the command line. A
        scheme is three characters; requiring it makes every value mean exactly one thing.
        """
        if value is None or (isinstance(value, str) and not value.strip()):
            if not required:
                return None
            err_console.print(f"[bold red]\\[ERROR][/bold red]: `{field}` is not set in the configuration file")
            err_console.print("Troubleshooting:")
            err_console.print(f" • Set `{field}` to a local directory containing {contents}")
            err_console.print(f" • Or to the GCS location holding it (e.g. gs://my-bucket/prefix), which is staged automatically")
            sys.exit(1)

        # The bucket's contents cannot be checked without a network round trip per field, so only the
        # form is validated here; a missing object is reported when it is staged.
        if is_gcs_path(value):
            return gcs_uri(value)

        if not Path(value).is_dir():
            err_console.print(f"[bold red]\\[ERROR][/bold red]: `{field}` in the configuration file is not a directory: {value}")
            err_console.print("Troubleshooting:")
            err_console.print(f" • Check that `{field}` is set to a valid directory containing {contents}")
            err_console.print(f" • If this data lives in a bucket, write the full URL with its scheme: gs://{str(value).lstrip('/')}")
            err_console.print(" • A bucket is recognised by the `gs://` prefix alone -- there is no flag to pass")
            sys.exit(1)
        return Path(value)

    if args.spatial_count or args.run_all:
        # puck_path is validated when set but not *required* here, because whether this run needs it
        # depends on the chemistry -- and the chemistry comes from the metadata sheet, which is only
        # loaded after this function returns. A Flex run never reads puck_path at all (its spatial
        # barcodes come from Takara bead-barcode files via the Trekker path, and main.py skips
        # run_spatial_positioning entirely), so requiring it here would force Flex users to invent a
        # value purely to get past this check. run_spatial_positioning -- the only consumer, and one
        # Flex never reaches -- reports an unset puck_path itself.
        PUCK_PATH = _validate_resource_path(
            PUCK_PATH, 'puck_path', 'pre-compiled barcode/spatial coordinate puck maps',
            required=False)
        # raw_barcodes_path is optional: it is only consulted when a puck map is absent and has to
        # be generated, so an unset value is not an error until that actually happens
        RAW_BARCODES_PATH = _validate_resource_path(
            RAW_BARCODES_PATH, 'raw_barcodes_path',
            'raw barcode/spatial coordinate files for puck map generation', required=False)

    if args.count or args.spatial_analysis or args.run_all:
        REF_PATH = _validate_resource_path(
            REF_PATH, 'reference_path', 'reference genome directories')

    # `software_path` is stageable out of GCS like the three resource paths above, but the trigger is
    # deliberately narrower: only an explicit `gs://` URI is staged, and a local directory stays valid
    # even when other fields name buckets.
    #
    # Cellranger and bcl2fastq are routinely preinstalled on whatever machine runs the pipeline -- the
    # --gcp VM image is built exactly that way -- so treating every
    # staged run as "the software must live in a bucket" would break those runs while buying nothing.
    #
    # A `gs://` URI here is acted on by itself, exactly as `auth_key_path` above is. This field led the
    # exists to disambiguate the bare `bucket/prefix` form, which cannot be told apart from a relative
    # local directory -- and that form is not honoured here at all, so there is nothing left for the
    # flag to decide. Requiring it anyway made a legitimate configuration unexpressible: software in a
    # bucket with the data local errored without the flag, and errored *with* it too, because
    # _validate_resource_path then demands that reference_path/puck_path/raw_barcodes_path be GCS
    # locations as well. Staging is announced below rather than silent, which is what the requirement
    # was really protecting against; it is also lazy (a warm software_cache.txt skips it entirely) and
    # reused by later runs sharing the output tree.
    if is_gcs_path(SOFTWARE_PATH):
        SOFTWARE_PATH = gcs_uri(SOFTWARE_PATH)
        # Worth saying out loud when nothing *else* comes from a bucket: the download is then the one
        # remote step in an operation the caller thinks is entirely local. Keyed on the other fields
        # rather than on a global switch, which since paths decide for themselves cannot answer the
        # question -- a run with a gs:// input_path and no flag stages plenty, and the note would be
        # telling it otherwise.
        if not any(is_gcs_path(other) for other in (INPUT_PATH, REF_PATH, PUCK_PATH, RAW_BARCODES_PATH)):
            console.print(f"[bold yellow]\\[NOTE][/bold yellow]: `software_path` is a GCS location, so Cellranger/bcl2fastq will be downloaded from {SOFTWARE_PATH} on first use")
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


    if not is_int(MEM_SIZE):
        err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for `memory` in configuration file: {MEM_SIZE!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check that `memory` is set to a integer value corresponding to the maximum amount of memory to be used by the pipeline (in GB)")
        hint = bool_value_hint(MEM_SIZE, 'settings.memory')
        if hint:
            err_console.print(hint)
        sys.exit(1)

    if not is_int(NUM_THREADS):
        err_console.print(f"[bold red]\\[ERROR][/bold red]: invalid entry for `threads` in configuration file: {NUM_THREADS!r}")
        err_console.print("Troubleshooting:")
        err_console.print(" • Check that `threads` is set to a integer value corresponding to the maximum number of CPU cores to be used by the pipeline")
        hint = bool_value_hint(NUM_THREADS, 'settings.threads')
        if hint:
            err_console.print(hint)
        sys.exit(1)

    # --fastqs exists so the pipeline never has to *guess* what kind of input it was handed.
    # Inferring "BCL run folder vs. already-demultiplexed FASTQ directory" by validating
    # `input_path` against the Illumina run-folder schema is ambiguous by construction: a
    # validation failure means either "these are FASTQs" or "this is a malformed BCL directory",
    # and the two demand opposite responses (silently skip mkfastq vs. abort with an error).
    # An explicit --fastqs makes the user's intent unambiguous, which in turn lets every
    # check in the input resolution below and in pipeline.py be strict instead of best-effort.
    #
    # --mkfastq asks for BCLs to be demultiplexed; --fastqs asserts that already happened. Refuse the
    # combination rather than silently honouring one and dropping the other. Checked here, ahead of
    # everything else, because it is a pure flag conflict: no config value can make it valid.
    if args.fastqs is not None and args.mkfastq:
        err_console.print("[bold red]\\[ERROR][/bold red]: --fastqs and --mkfastq are mutually exclusive")
        err_console.print("Troubleshooting:")
        err_console.print(" • Drop --fastqs to demultiplex BCLs from `input_path` with cellranger mkfastq")
        err_console.print(" • Drop --mkfastq to run the downstream stages on the FASTQs you provided")
        sys.exit(1)

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

    if ALERTS:
        SLACK_TOKEN = _resolve_slack_token(SLACK_TOKEN)
    
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

    # Resolve the run's output directory as output_path/<BCL_ID>, accepting a configured path that
    # already ends in the BCL ID as pointing at the run directory itself. Without this, re-pointing
    # `output_path` at an existing run's directory (the natural thing to do when resuming or
    # inspecting one) silently nests a second copy at <BCL_ID>/<BCL_ID>/ and the run finds none of
    # its earlier outputs. This is the same convention resolve_bcl_dir() applies to `input_path`.
    BCL_ID = args.bcl
    OUT_PATH = Path(OUT_PATH)
    if OUT_PATH.name != BCL_ID:
        OUT_PATH = OUT_PATH / BCL_ID
    OUT_PATH.mkdir(exist_ok=True, parents=True)

    # define and create output sub-directories
    LOG_PATH = OUT_PATH / "log"
    METADATA_PATH = OUT_PATH / "metadata"
    OUTPUT_PATH = OUT_PATH / "output"
    TMP_PATH = OUT_PATH / "tmp"
    RUNTIME_LOG = LOG_PATH / "runtime.log"
    SUMMARY_PATH = METADATA_PATH / "metadata_summary.csv"
    SAMPLESHEET_PATH = METADATA_PATH / "samplesheet.csv"
    MKFASTQ_OUTS = OUTPUT_PATH / "mkfastq"
    COUNT_OUTS = OUTPUT_PATH / "count"
    CELLBENDER_OUTS = OUTPUT_PATH / "cellbender"
    SPATIAL_COUNT_OUTS = OUTPUT_PATH / "spatial_barcodes"
    SPATIAL_ANALYSIS_OUTS = OUTPUT_PATH / "spatial_analysis"
    FLEX_OUTS = OUTPUT_PATH / "flex"

    for path in [LOG_PATH, METADATA_PATH, TMP_PATH, OUTPUT_PATH]:
        path.mkdir(exist_ok=True)

    # ------------------------------------------------------------------------------------- #
    #                      resolve where this run reads its reads from                       #
    # ------------------------------------------------------------------------------------- #
    #
    # `paths.input_path` is dual-purpose, exactly like reference_path/puck_path/raw_barcodes_path: a
    # local directory of run folders normally, and the `gs://` prefix they are staged out of under
    # scheme. Which one it is comes from the value, and from nothing else.
    #
    # When staging, the configured value is kept as INPUT_BUCKET -- what to download -- and INPUT_PATH
    # is re-pointed at the local directory the download lands in, which is where every stage then
    # looks. Nothing downstream has to know the difference: INPUT_PATH is always a local Path holding
    # run folders, so resolve_bcl_dir, the Flex samplesheet and the log lines are identical either way.
    #
    # Resolved here, after the output tree exists, because an unset `settings.gcs_download_dest`
    # defaults into it.
    INPUT_BUCKET = None
    STAGE_FASTQ_INPUT = False
    FASTQ_INPUT = args.fastqs

    # Same rule as the resource paths above: a `gs://` input_path names a bucket on its own say-so,
    # and a bare `bucket/prefix` is not accepted at all, being indistinguishable from a
    # relative directory. Anything else is a local directory of run folders.
    #
    # Note what a bucket-backed input means for a `--slurm` run beyond staging itself: ./slidr
    # uses it to decide that the job must self-stage, and so submits workflow/bash/slidr_slurm.sh
    # (which checks gcloud, exports SLIDR_OUTPUT_PATH and syncs on the compute node) rather than
    # wrapping main.py directly. The launcher infers it from a gs:// input_path for that reason.
    input_is_bucket = is_gcs_path(INPUT_PATH)

    if input_is_bucket:
        INPUT_BUCKET = gcs_uri(INPUT_PATH)

        # Where the download lands. `settings.gcs_download_dest` names it explicitly; left unset it
        # defaults inside this run's own directory, which keeps a staged run self-contained -- on a VM
        # or a cluster that tree is the one place the pipeline is guaranteed to be able to write -- and
        # means no machine-specific configuration is needed to stage anywhere.
        #
        # Deliberately a sibling of output/ rather than output/data, which is where the other staged
        # resources live (output/reference, output/pucks, output/barcodes, output/software). main.py
        # uploads the whole of output/ to `settings.output_bucket` on success, and these are raw inputs
        # that came *out* of a bucket: putting them in there would re-upload the entire sequencing run,
        # routinely hundreds of GB, every time a run finished.
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

        # A bare --fastqs means "the reads are this run's input folder". Staged, that folder is
        # <gcs_download_dest>/<BCL_ID> -- the local copy of <input_path>/<BCL_ID> -- whether it holds
        # BCLs or FASTQs, so it resolves exactly as a run folder does. It cannot be checked here: it
        # is downloaded when the run needs it (main.py), which is also where it is then validated.
        if FASTQ_INPUT == FASTQS_USE_INPUT_PATH:
            FASTQ_INPUT = INPUT_PATH / BCL_ID
            STAGE_FASTQ_INPUT = True

    else:
        # a local input_path stages nothing, so the download destination names nowhere; keep it out of
        # the run's recorded config rather than reporting a directory this run will never write to
        GCS_DOWNLOAD_DEST = None

    # Validate the FASTQ directory supplied with --fastqs. Skipped for a staged bare --fastqs, whose
    # directory does not exist yet; everything else is checked here, before any stage runs.
    if FASTQ_INPUT is not None and not STAGE_FASTQ_INPUT:
        # A bare --fastqs makes `input_path` the only source of the reads, and therefore mandatory.
        # Report it here, where the reason is obvious, rather than letting the generic input_path
        # error below suggest passing the very flag that was just used.
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

    if INPUT_PATH is None or not Path(INPUT_PATH).is_dir():
        # with --fastqs the reads are already demultiplexed, so no BCL run folder is required;
        # point INPUT_PATH at the FASTQ directory so the remaining consumers of it (the Flex
        # cellranger-multi samplesheet, log messages) still resolve to the real input
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
        # Which inputs come from a bucket is now a per-field fact, so the summary records it per field
        # rather than as one yes/no, since the fields are independent:
        # whether a bare `bucket/prefix` was read as a bucket.
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

    # Which stages this invocation asked for. Computed here so the log can state the plan up front --
    # far more useful than four near-identical "Skipping X: --X not provided" lines discovered one at a
    # time while waiting. A stage may still be skipped later because its outputs are already present;
    # that is reported when it happens, since it is not knowable yet.
    STAGE_FLAGS = [
        ('mkfastq', args.mkfastq),
        ('count', args.count),
        ('cellbender', args.cellbender and not args.no_cellbender),
        ('spatial-count', args.spatial_count),
        ('spatial-analysis', args.spatial_analysis),
    ]
    if args.run_all:
        REQUESTED_STAGES = [name for name, _ in STAGE_FLAGS
                            if not (name == 'cellbender' and args.no_cellbender)]
    else:
        REQUESTED_STAGES = [name for name, requested in STAGE_FLAGS if requested]
    if args.fastqs is not None and 'mkfastq' in REQUESTED_STAGES:
        REQUESTED_STAGES.remove('mkfastq')
    SKIPPED_STAGES = [name for name, _ in STAGE_FLAGS if name not in REQUESTED_STAGES]

    # set up runtime logfile
    RUN_WIDTH = 74

    def _banner(text: str, fill: str = '═') -> str:
        """Centre `text` in a full-width rule."""
        pad = max(0, RUN_WIDTH - len(text) - 2)
        left = pad // 2
        return f"{fill * left} {text} {fill * (pad - left)}"

    def _field(label: str, value) -> str:
        return f"  {label:<14} {value}"

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
        # keyed on the input actually naming a bucket, not on the flag: a gs:// input_path stages
        # itself now, and reporting the flag would say nothing about where this run's reads come from
        if INPUT_BUCKET is not None:
            logfile.write(_field("Staging", f"{INPUT_BUCKET} → {GCS_DOWNLOAD_DEST}") + "\n")
        logfile.write(_field("Stages", ', '.join(REQUESTED_STAGES) or "none requested") + "\n")
        if SKIPPED_STAGES:
            logfile.write(_field("", f"(not requested: {', '.join(SKIPPED_STAGES)})") + "\n")
        if provided_args:
            logfile.write(_field("Arguments", ' '.join(provided_args)) + "\n")
        logfile.write("─" * RUN_WIDTH + "\n")

    console.print('[bold green]\\[SUCCESS][/bold green]: Setup complete')

    # package global constants into a dictionary
    cfg = {
        'root_path': ROOT_DIR,
        'run_path': OUT_PATH,
        'requested_stages': REQUESTED_STAGES,
        'output_path': OUTPUT_PATH,
        'input_path': INPUT_PATH,
        'fastq_input': FASTQ_INPUT,
        'stage_fastq_input': STAGE_FASTQ_INPUT,
        'script_path': SCRIPT_PATH,
        'metadata_path': METADATA_PATH,
        'log_path': LOG_PATH,
        'tmp_path': TMP_PATH,
        'software_path': SOFTWARE_PATH,
        'metadata_src': METADATA_SRC,
        'raw_barcodes_path': RAW_BARCODES_PATH,
        'puck_path': PUCK_PATH,
        'ref_path': REF_PATH,
        'ref_genome': REF_GENOME,
        'software_cache_file': SOFTWARE_CACHE_FILE,
        'auth_key_path': AUTH_KEY_PATH,
        # the two halves of a bucket-backed `paths.input_path`: the gs:// prefix to download from,
        # and the local directory it lands in (which is what 'input_path' above is set to). Both None
        # for a local run, where 'input_path' is the configured directory itself.
        'input_bucket': INPUT_BUCKET,
        'gcs_download_dest': GCS_DOWNLOAD_DEST,
        'output_bucket': OUTPUT_BUCKET,
        # set only when the caller already chose the exact folder to upload to (./slidr --gcp); None
        # otherwise, and the run picks a free folder under `output_bucket` itself
        'output_dest': OUTPUT_DEST,
        'num_threads': NUM_THREADS,
        'mem_size': MEM_SIZE,
        'alerts': ALERTS,
        'slack_token': SLACK_TOKEN,
        'metadata_specs': METADATA_SPECS,
        'runtime_log': RUNTIME_LOG,
        'summary_log': SUMMARY_LOG,
        'summary_path': SUMMARY_PATH,
        'samplesheet_path': SAMPLESHEET_PATH,
        'mkfastq_outs': MKFASTQ_OUTS,
        'count_outs': COUNT_OUTS,
        'cellbender_outs': CELLBENDER_OUTS,
        'spatial_count_outs': SPATIAL_COUNT_OUTS,
        'spatial_analysis_outs': SPATIAL_ANALYSIS_OUTS,
        'flex_outs': FLEX_OUTS,
        'generate_bam': GENERATE_BAM,
        'cellbender_epochs': CELLBENDER_EPOCHS,
        'cellbender_rate': CELLBENDER_RATE,
        'cellbender_cells': CELLBENDER_CELLS,
        'cellbender_droplets': CELLBENDER_DROPLETS,
        'spatial_downsampling': SPATIAL_DOWNSAMPLING,
        'percent_umi_filtering': PERCENT_UMI_FILTERING,
        'emptydrops_min_umis': EMPTYDROPS_MIN_UMI,
        'flex_probe_set': FLEX_PROBE_SET,
        'flex_r1_path': FLEX_R1_PATH,
        'flex_r2_path': FLEX_R2_PATH,
        'flex_gex_fastqs': FLEX_GEX_FASTQS,
        'start_time': START_TIME,
        'bcl_id': BCL_ID
    }

    return args, cfg


args, cfg = _load()