import gspread
import json
import time
import os
import re
import shutil
import sys
import csv
import h5py
import subprocess
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

from google.oauth2.service_account import Credentials
from datetime import datetime
from pathlib import Path
from rich.console import Console
from urllib.parse import urlparse, parse_qs

# import from sibling directories
sys.path.append(str(Path(__file__).parent))
from config import cfg, is_gcs_path, is_google_sheet


# unpack globals from config dictionary
TMP_PATH = cfg['tmp_path']
SOFTWARE_PATH = cfg['software_path']
INPUT_PATH = cfg['input_path']
FASTQ_INPUT = cfg['fastq_input']
RUN_PATH = cfg['run_path']
OUTPUT_PATH = cfg['output_path']
SCRIPT_PATH = cfg['script_path']
ROOT_PATH = cfg['root_path']
SOFTWARE_CACHE_FILE = cfg['software_cache_file']
AUTH_KEY_PATH = cfg['auth_key_path']
METADATA_SRC = cfg['metadata_src']
SUMMARY_PATH = cfg['summary_path']
METADATA_SPECS = cfg['metadata_specs']
RUNTIME_LOG = cfg['runtime_log']
ALERTS = cfg['alerts']
SLACK_TOKEN = cfg['slack_token']
BCL_ID = cfg['bcl_id']
START_TIME = cfg['start_time']
SUMMARY_LOG = cfg['summary_log']
STAGE_GCS = cfg['stage_gcs']

# Local copy of a `gs://` software_path, once staged. Memoized because the tree is large (a single
# cellranger release is a couple of GB) and test_and_install_software() is called once per sample per
# tool, so the download must happen at most once per run.
_SOFTWARE_DEST = None

WARN_THRESHOLD_GB = 100
CRITICAL_THRESHOLD_GB = 500

# cellranger/bcl2fastq FASTQ filenames: "<library>_S<n>_L<lane>_<R1|R2|I1|I2>_<chunk>.fastq.gz"
FASTQ_NAME_RE = re.compile(r'^(?P<lib>.+)_S\d+_L(?P<lane>\d+)_(?P<read>[RI][12])_\d+\.fastq\.gz$')

# Library naming convention required of an externally supplied FASTQ directory (--fastqs).
#
# mkfastq writes each library into its own mkfastq/<library>/ folder, so which library a FASTQ
# belongs to is unambiguous from its location. A directory handed in with --fastqs has no such
# structure, so the library has to be recoverable from the filename instead:
#
#   <sample name>_GEX_S1_L001_R1_001.fastq.gz   gene-expression library
#   <sample name>_SB_S1_L001_R1_001.fastq.gz    spatial-barcode library
#   <sample name>_S1_L001_R1_001.fastq.gz       gene-expression library (bare form; this is how
#                                               cellranger mkfastq itself names it)
#
# '-' works as the separator too, the mode token is case-insensitive (so mkfastq's own '_sb'
# spatial libraries are recognised), and trailing '_'/'-'-delimited extras are allowed.
FASTQ_MODES = ('GEX', 'SB')

if ALERTS:
    from slack_sdk import WebClient
    from slack_sdk.errors import SlackApiError


# Cellranger version used when a sample's `Cellranger` metadata column is blank or absent.
DEFAULT_CELLRANGER_VERSION = '8.0.1'

# Full Cellranger versions for the shorthand major versions users write in the metadata sheet ("V8",
# "V7"). test_and_install_software() matches an installation directory named
# `cellranger[-_v]<version>` exactly, so a bare major like "8" would not find `cellranger-8.0.1` --
# the major has to be expanded to the specific release that is installed.
#
# EDIT THIS TABLE to match the Cellranger releases available under `software_path`. A major version
# that is not listed here is passed through as written, so `V7.1.0` works without a table entry.
CELLRANGER_VERSIONS = {
    '6': '6.1.2',
    '7': '7.2.0',
    '8': '8.0.1',
    '9': '9.0.1',
}


def resolve_cellranger_version(
    value=None,
    default: str = DEFAULT_CELLRANGER_VERSION
) -> str:
    """
    Translate a `Cellranger` metadata value into a version string for test_and_install_software().

    Accepts the shorthand users actually write in the sheet -- "V8", "v8", "8" -- and expands it to
    the matching full release from CELLRANGER_VERSIONS. A value that already names a full version
    ("8.0.1", "V8.0.1") is used as written, as is an unrecognized major, so a release missing from the
    table can still be requested explicitly.

    Inputs:
     - value:   raw metadata value; None, NaN or blank selects `default`
     - default: version to use when the metadata declares none. The Flex path passes its own, since
                cellranger multi needs a newer release than the standard count path
    Output:
     - a version string such as '8.0.1'
    """

    # a blank cell arrives as NaN from pandas rather than as an empty string
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return default

    # strip the conventional leading "V"/"v" and any surrounding whitespace
    text = str(value).strip().lstrip('Vv').strip()
    if not text:
        return default

    # a bare major version needs expanding; anything already carrying a '.' is a full version
    return CELLRANGER_VERSIONS.get(text, text)


def _restore_exec_bits(root: Path) -> None:
    """
    Add the execute bit to everything under a staged software tree.

    `gcloud storage cp` does not carry POSIX permissions unless they were recorded at upload time, so
    a downloaded cellranger install arrives mode 644 and every executable in it is unrunnable -- which
    also makes `find -type f -executable` (below) match nothing, so the failure would surface as a
    misleading "software not found". Which of the thousands of files in a cellranger release are meant
    to be executable cannot be recovered from the objects, so the bit is set on all of them. This is a
    private, run-local copy under the output tree, so marking a data file executable has no
    consequence beyond the cosmetic.

    Inputs:
     - root: the staged software directory to fix up
    """

    for path in root.rglob('*'):
        try:
            path.chmod(path.stat().st_mode | 0o111)
        except OSError:
            # a single unreadable entry should not abort the run; if it was one of the executables
            # the caller needs, the search below reports it as missing
            continue


def software_dir() -> Path:
    """
    Local directory to scan for external executables.

    `paths.software_path` is normally a local directory and is returned as-is. When it names a
    `gs://` location (only permitted alongside --stage-gcs, see config.py) the tree is downloaded into
    the run's output directory on first use and the local copy is returned instead, because the
    executable search is a `find` over a real filesystem and cannot run against a bucket.

    Staging is lazy and memoized: it happens on the first call that actually needs to scan -- which a
    warm software_cache.txt avoids entirely -- and at most once per run.

    Output:
     - a directory that can be scanned for executables
    """

    global _SOFTWARE_DEST

    if not is_gcs_path(SOFTWARE_PATH):
        return Path(SOFTWARE_PATH)

    if _SOFTWARE_DEST is not None:
        return _SOFTWARE_DEST

    dest = Path(OUTPUT_PATH) / "software"
    # reuse a copy a previous run in this output tree already brought down; the software is far larger
    # than the sequencing data staging bothers to skip, and it does not change between runs
    if dest.is_dir() and any(dest.iterdir()):
        log_write(f"  Software already staged at {dest}")
    else:
        dest.mkdir(parents=True, exist_ok=True)
        log_write(f"  Staging software from {SOFTWARE_PATH} (this can take a while)... ", terminator="")
        # copy the prefix's *contents* into dest: `cp -r gs://p dest` would nest them under
        # dest/<last path segment>, which the executable search would still find but which makes the
        # staged tree's shape depend on how the prefix happens to be named
        stage_from_gcs(
            f'{SOFTWARE_PATH}/*',
            dest,
            recursive=True,
            description='software installations'
        )
        _restore_exec_bits(dest)
        log_write("Done.")

    _SOFTWARE_DEST = dest
    return dest


def test_and_install_software(
    software: str,
    version: str = ""
) -> Path:
    """
    Search for an external executable and try to install it if not found locally

    Inputs:
     - software: name of the executable to search for
    Output:
     - a pathlib.Path object containing the path to specified executable
    """

    # create the software cache file if does not exist
    Path(SOFTWARE_CACHE_FILE).touch(exist_ok=True)

    # search the local software cache for matching paths
    with open(SOFTWARE_CACHE_FILE, 'r') as file:
        for line in file:
            line = line.strip()
            if line.endswith(software) and version in line:
                # if a matching software path is an executable, return the path
                if Path(line).is_file() and os.access(Path(line), os.X_OK):
                    print(f"  Using {software} at {line}")
                    return Path(line)
    
    # if slidr cannot find path to required software, it is installed automatically in $HOME/.local
    install_path = Path.home() / ".local" / "slidr" / "bin" / software

    # check if software has been installed automatically before
    if install_path.is_file() and os.access(install_path, os.X_OK):
        with open(SOFTWARE_CACHE_FILE, 'a') as file:
            file.write(str(install_path) + '\n')
        print(f"  Using {software} at {install_path}")
        return install_path

    # Resolve the directory to scan, downloading it first if `software_path` names a bucket. Deferred
    # until here on purpose: the cache lookups above are what a warm software_cache.txt hits, and they
    # must not pay for a multi-GB staging download to answer from a file they already have.
    search_root = software_dir() if SOFTWARE_PATH is not None else None

    # format the versioned name if version is provided
    if version:
        version_escaped = re.escape(version)
        # `[-_v][-_v]*` rather than the more obvious `[-_v]+`: find's -regex takes a *basic* regex on
        # BSD/macOS, where `+` is a literal plus rather than a quantifier, so `[-_v]+` silently matched
        # nothing there while working under GNU find on Linux. `*` is a quantifier in both dialects, as
        # is the `\.` that re.escape() puts in the version, so this form resolves on either.
        regex = rf".*/{software}[-_v][-_v]*{version_escaped}$"
        cmd = [
            "find",
            str(search_root),
            "-regex",
            regex,
            "-type",
            "d"
        ]
    else:
        # `-name` (not a bare `software`, which find would read as a second directory to search) and
        # `-perm -u+x` rather than GNU's `-executable`, which BSD find rejects outright. Both matter
        # more for a staged tree than a local one: everything staged out of GCS has the execute bit
        # restored wholesale (see _restore_exec_bits), so an unfiltered search would happily return the
        # first data file it walked past.
        cmd = [
            "find",
            str(search_root),
            "-name",
            software,
            "-type",
            "f",
            "-perm",
            "-u+x"
        ]

    # try searching the entire provided software directory for matching software (slow)
    if search_root is not None and Path(search_root).is_dir():
        print(f"Scanning {search_root} for {software} executable. This may take a while...")
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            # if found, add path to the software cache
            execs = result.stdout.strip().splitlines()
            if len(execs) > 0:
                executable_path = Path(execs[0])
                if version:
                    executable_path = executable_path / "bin" / software
                    if not executable_path.is_file() or not os.access(executable_path, os.X_OK):
                        log_write(f"[ERROR]: found {software} at {executable_path}, but it is not an executable file")
                        log_write("Troubleshooting:")
                        log_write(f" • Check the file exists and is executable: `ls -l {executable_path}`")
                        log_write(f" • Make it executable if the permission bit is missing: `chmod +x {executable_path}`")
                        log_write(f" • An extracted-but-incomplete {software} install often looks like this; re-extract the release tarball")
                        log_write(f" • Or pin a known-good path by adding it to {SOFTWARE_CACHE_FILE}")
                        sys.exit(1)
                with open(SOFTWARE_CACHE_FILE, 'a') as file:
                    file.write(str(executable_path) + '\n')
                print(f"  Using {software} at {executable_path}")
                return Path(executable_path)
        except Exception as exc:
            log_write(f"[WARNING]: Could not search {search_root} for {software}: {exc}")

    # software not found locally, attempting to install automatically
    match software:
        case "cellranger":
            # proprietary software, cannot install automatically
            log_write("[ERROR]: Cellranger software not found on the system.\n")
            log_write("Troubleshooting:")
            log_write(" • Follow the steps on the 10x Genomics website to install Cellranger: https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links")
            sys.exit(1)

        case "bcl2fastq2":
            # proprietary software, cannot install automatically
            log_write("[ERROR]: bcl2fastq2 software not found on the system.\n")
            log_write("Troubleshooting:")
            log_write(" • Follow the steps on the Illumina website to install bcl2fastq2: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html")
            sys.exit(1)

        case "julia":
            # checking for system Julia
            julia_path = shutil.which("julia")
            if julia_path is not None:
                with open(SOFTWARE_CACHE_FILE, 'a') as file:
                    file.write(str(julia_path) + '\n')
                
                print(f'Using Julia at {julia_path}')
                return Path(julia_path)

            # installing Julia automatically
            install_path.mkdir(parents=True, exist_ok=True)
            print("Installing Julia...")
            env = os.environ.copy()
            env["JULIAUP_DEPOT_PATH"] = str(install_path)

            # download and execute the Julia installer
            try:
                curl_proc = subprocess.run(
                    ["curl", "-fsSL", "https://install.julialang.org"],
                    capture_output=True,
                    check=True,
                    env=env
                )
                sh_result = subprocess.run(
                    ["sh", "-s", "--", "-y"],
                    input=curl_proc.stdout,
                    env=env
                )
                # verify the result
                if sh_result.returncode == 0:
                    julia_path = install_path / "bin" / "julia"
                    with open(SOFTWARE_CACHE_FILE, 'a') as file:
                        file.write(str(julia_path) + '\n')
                    print(f"Installation successful - using {software} at {julia_path}")
                    return Path(julia_path)
                else:
                    log_write(f"[ERROR]: the Julia installer exited with code {sh_result.returncode}")
                    log_write("Troubleshooting:")
                    log_write(" • Install Julia yourself with juliaup: https://julialang.org/downloads/")
                    log_write(f" • Then add the path to the `julia` binary to {SOFTWARE_CACHE_FILE} so slidr uses it directly")
                    log_write(f" • Check there is free space and write permission on {install_path}")
                    sys.exit(1)
            except Exception as e:
                log_write(f"[ERROR]: Julia installation failed: {e}")
                log_write("Troubleshooting:")
                log_write(" • The installer is downloaded with curl, so this is usually no network access or a proxy blocking https://install.julialang.org")
                log_write(" • Install Julia yourself with juliaup: https://julialang.org/downloads/")
                log_write(f" • Then add the path to the `julia` binary to {SOFTWARE_CACHE_FILE} so slidr uses it directly")
                sys.exit(1)
        case _:
            log_write(f'[INTERNAL ERROR]: Unrecognized software name: {software}')
            log_write("Troubleshooting:")
            log_write(" • This is a bug in slidr, not a problem with your setup or data")
            log_write(" • test_and_install_software() only knows `cellranger`, `bcl2fastq2` and `julia`")
            log_write(" • Please report it at https://github.com/kabanovskyd/slidr/issues, quoting the software name above")
            sys.exit(1)


def read_service_account_key(auth_key_path) -> dict:
    """
    Return the parsed service-account key `paths.auth_key_path` names, without ever writing it to disk.

    A `gs://` value is read with `gcloud storage cat` and parsed in memory, so the credential exists
    only for the life of the process. This replaces downloading it into the run's tmp directory: that
    copy had to be created with the right mode, removed afterwards and kept out of the uploaded output
    tree, and every one of those is a chance to leave a private key sitting on a shared filesystem or in
    a bucket. Not writing it at all removes the class of mistake rather than guarding each instance, and
    it is why neither launcher script stages a key file any more.

    A local path is read as before -- a key already on the machine is not made safer by refusing to
    open it.

    Errors quote gcloud's stderr but never the object's contents: this function only ever handles a
    credential, and an error message is the one place one could end up in a log that gets shared.

    Inputs:
     - auth_key_path: the configured value -- a local Path, or a gs:// URI string
    Output:
     - the key as a dict, ready for Credentials.from_service_account_info
    """

    if not is_gcs_path(auth_key_path):
        source = str(auth_key_path)
        try:
            payload = Path(auth_key_path).read_text()
        except OSError as exc:
            log_write(f"[ERROR]: could not read the service account key at {auth_key_path}: {exc}")
            log_write("Troubleshooting:")
            log_write(" • Check the file exists and this user can read it")
            log_write(" • Or set `paths.auth_key_path` to a gs:// object holding the key, read without being downloaded")
            sys.exit(1)
    else:
        source = str(auth_key_path)
        cmd = ['gcloud', 'storage', 'cat', source]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            log_write(f"[ERROR]: `auth_key_path` is a gs:// URI but `gcloud` was not found on PATH, so the key cannot be read")
            log_write("Troubleshooting:")
            log_write(" • Install the Google Cloud CLI: https://cloud.google.com/sdk/docs/install")
            log_write(" • On a cluster, check whether it needs loading first (e.g. `module load google-cloud-sdk`)")
            log_write(" • Or point `paths.auth_key_path` at a local copy of the key instead")
            sys.exit(1)

        if proc.returncode != 0:
            detail = (proc.stderr or '').strip().splitlines()
            log_write(f"[ERROR]: could not read the service account key from {source} (`gcloud storage cat` exited {proc.returncode})")
            for line in detail[-3:]:
                log_write(f"  {line}")
            log_write("Troubleshooting:")
            log_write(f" • Check the object exists and is a single file, not a prefix: `gcloud storage ls {source}`")
            log_write(" • Check the active account can read it: `gcloud auth print-access-token`")
            log_write(" • Or point `paths.auth_key_path` at a local copy of the key instead")
            sys.exit(1)

        payload = proc.stdout

    try:
        info = json.loads(payload)
    except ValueError as exc:
        log_write(f"[ERROR]: the service account key at {source} is not valid JSON: {exc}")
        log_write("Troubleshooting:")
        log_write(" • The key must be the JSON file downloaded from IAM & Admin → Service Accounts → Keys")
        log_write(" • A .p12 key will not work here; create a JSON one instead")
        sys.exit(1)

    if not isinstance(info, dict) or 'client_email' not in info:
        log_write(f"[ERROR]: the file at {source} is not a service-account key (no `client_email` field)")
        log_write("Troubleshooting:")
        log_write(" • An OAuth client-secret file will not work here -- the key must be a *service account* key")
        log_write(" • Create one under IAM & Admin → Service Accounts → Keys in the Google Cloud console")
        sys.exit(1)

    return info


def load_metadata(
    input_source: str
) -> pd.DataFrame:
    """
    Read a provided input source (Google Sheet or TSV) and return raw data

    Inputs:
     - input_source: URL to a google sheet OR path to TSV containing sample metadata
    Output:
     - a pd.DataFrame object containing sample metadata loaded from the provided source

    """
    METADATA_OPT = [
        'Merge RNA From BCL',
        'Merge Spatial From BCL',
        'Add RNA Index',
        'Add SB Index',
        'Add Puck ID',
        'Cellranger',
        'Flex Probe Barcode IDs'
    ]

    # download data from a Google Sheet. Uses the same predicate as config.py's auth_key_path check,
    # so the key is demanded for exactly the sources that reach this branch.
    if is_google_sheet(input_source):
        # define the interaction scope
        SCOPES = [
            "https://www.googleapis.com/auth/spreadsheets",
            "https://www.googleapis.com/auth/drive"
        ]

        # Authenticate and create a Google Sheet client. The key is parsed in memory rather than read
        # from a file gspread is handed, so a gs:// key never lands on this machine's disk.
        key_info = read_service_account_key(AUTH_KEY_PATH)
        service_account = key_info.get('client_email', '(unknown)')
        try:
            creds = Credentials.from_service_account_info(key_info, scopes=SCOPES)
            client = gspread.Client(auth=creds)
        except Exception as exc:
            log_write(f"[ERROR]: could not authorize Google Cloud account with credentials from {AUTH_KEY_PATH}: {exc}")
            log_write("Troubleshooting:")
            log_write(f" • The key parsed as a service-account key for {service_account}, so it is the credential itself that was rejected")
            log_write(" • Check the key has not been deleted or disabled under IAM & Admin → Service Accounts → Keys")
            log_write(" • Check the Google Sheets and Google Drive APIs are enabled for the service account's project")
            sys.exit(1)

        # open the Google Sheet
        try:
            sheet = client.open_by_url(input_source)
        except Exception as exc:
            log_write(f"[ERROR]: could not open Google Sheet with provided URL ({input_source}): {exc}")
            log_write("Troubleshooting:")
            log_write(f" • Share the sheet with the service account's email address ({service_account}), with at least Viewer access")
            log_write(" • Check `settings.metadata_source` is the full sheet URL, copied from the browser address bar")
            log_write(" • Check the Google Sheets and Google Drive APIs are enabled for the service account's project")
            sys.exit(1)

        # extract the worksheet (tab) specified by the URL's gid, defaulting to the first tab
        parsed_url = urlparse(input_source)
        gid = parse_qs(parsed_url.fragment).get('gid') or parse_qs(parsed_url.query).get('gid')
        gid = int(gid[0]) if gid else None

        try:
            worksheet = sheet.get_worksheet_by_id(gid) if gid is not None else sheet.get_worksheet(0)
        except Exception as exc:
            log_write(f"[ERROR]: worksheet (gid={gid}) not found in the Google Sheet at {input_source}: {exc}")
            log_write("Troubleshooting:")
            log_write(" • The worksheet tab is taken from the `#gid=` fragment of the URL; open the tab you want in the browser and re-copy the URL")
            log_write(" • Check the tab has not been deleted or renamed since the URL was saved")
            log_write(" • Drop the `#gid=...` fragment from `settings.metadata_source` to use the first tab instead")
            sys.exit(1)

        # read entries as a DataFrame
        raw = worksheet.get_all_records(numericise_ignore=['all'])
        if not raw:
            log_write(f"[ERROR]: worksheet '{worksheet.title}' in {input_source} contains no rows")
            log_write("Troubleshooting:")
            log_write(" • Check the `#gid=` fragment in the URL points at the tab that actually holds the sample metadata")
            log_write(" • The first row must be the column headers, with at least one data row beneath it")
            log_write(" • A tab whose header row is blank reads as empty even when data is present further down")
            sys.exit(1)
        records = pd.DataFrame(raw)

        for field in METADATA_OPT:
            if records[field].notna().any():
                METADATA_SPECS.append(field)

        # subset the DataFrame to the required columns
        if not set(METADATA_SPECS).issubset(records.columns):
            missing_cols = set(METADATA_SPECS) - set(records.columns)
            log_write(f"[ERROR]: Google Sheet input missing required columns: {', '.join(sorted(missing_cols))}")
            log_write("Troubleshooting:")
            log_write(f" • Add the missing column(s) to the header row of worksheet '{worksheet.title}'")
            log_write(f" • Column names are matched exactly, including case and spaces; the full required set is: {', '.join(METADATA_SPECS)}")
            log_write(" • Check the header row is the *first* row of the tab -- a title or blank row above it shifts the headers out of view")
            sys.exit(1)

        return records[METADATA_SPECS]

    # read input metadata
    elif Path(input_source).is_file():
        # parse input TSV
        if input_source.endswith('.tsv'):
            metadata = pd.read_table(input_source, dtype=str)
            for field in METADATA_OPT:
                if field in metadata.columns and metadata[field].notna().any():
                    METADATA_SPECS.append(field)
            if not set(METADATA_SPECS).issubset(metadata.columns):
                missing_cols = set(METADATA_SPECS) - set(metadata.columns)
                log_write(f"[ERROR]: {input_source} missing required columns: {', '.join(sorted(missing_cols))}")
                log_write("Troubleshooting:")
                log_write(f" • Add the missing column(s) to the header line of {input_source}")
                log_write(f" • Column names are matched exactly, including case and spaces; the full required set is: {', '.join(METADATA_SPECS)}")
                log_write(" • A .tsv file must be tab-separated -- a comma-separated file with a .tsv extension reads as one single column")
                sys.exit(1)
            return metadata[METADATA_SPECS]

        # parse input CSV
        elif input_source.endswith('.csv'):
            metadata = pd.read_csv(input_source, dtype=str)
            for field in METADATA_OPT:
                if field in metadata.columns and metadata[field].notna().any():
                    METADATA_SPECS.append(field)
            log_write('[WARNING]: metadata passed in as CSV, please make sure that comma-containing values in `Lane`/`SB Lane` columns are safely wrapped in double quotes')
            log_write(' • An unquoted "1,2" in a lane column splits into two CSV fields and shifts every later column by one')
            log_write(' • Passing the metadata as .tsv avoids the problem entirely')
            if not set(METADATA_SPECS).issubset(metadata.columns):
                missing_cols = set(METADATA_SPECS) - set(metadata.columns)
                log_write(f"[ERROR]: {input_source} missing required columns: {', '.join(sorted(missing_cols))}")
                log_write("Troubleshooting:")
                log_write(f" • Add the missing column(s) to the header line of {input_source}")
                log_write(f" • Column names are matched exactly, including case and spaces; the full required set is: {', '.join(METADATA_SPECS)}")
                log_write(" • Check no unquoted comma in an earlier row shifted the columns (see the warning above)")
                sys.exit(1)
            if 'Flex' in metadata['Chemistry'].values:
                METADATA_SPECS.append('Flex Probe Barcode IDs')
            return metadata[METADATA_SPECS]
        else:
            log_write(f"[ERROR]: unrecognized metadata file format (should be .tsv or .csv): {input_source}")
            log_write("Troubleshooting:")
            log_write(" • Rename the file so its extension is .tsv (tab-separated) or .csv (comma-separated) -- the format is chosen by extension")
            log_write(" • An .xlsx spreadsheet must be exported first (File → Download → Tab-separated values)")
            log_write(" • Or set `settings.metadata_source` to a Google Sheet URL instead")
            sys.exit(1)
    else:
        log_write(f"[ERROR]: Cannot find experimental metadata at {input_source}.")
        log_write("Troubleshooting:")
        log_write(" • Set `settings.metadata_source` to either an absolute path to a .tsv/.csv file, or a full Google Sheet URL")
        log_write(" • A Google Sheet URL is only recognized when it starts with https://docs.google.com/spreadsheets")
        log_write(" • Check the path for typos and that the file is readable from this machine")
        log_write(" • On a --gcp or staged --slurm run the path must exist on the compute node, not on your laptop")
        sys.exit(1)


# Per-module output layout and dependency graph used by need_run_module().
#   dir             : subdirectory under OUTPUT_PATH that holds each sample's <sample>/ output
#                     folder (these are the real on-disk names from config.py's *_OUTS, which do
#                     NOT all match the module key -- e.g. spatial_positioning writes to
#                     spatial_barcodes/).
#   outputs         : filenames that must all be present for a sample's output to count as
#                     complete. None means "validate FASTQ presence instead" (mkfastq).
#   mtime_inputs    : upstream modules whose outputs, when present on disk, must be OLDER than
#                     this module's outputs for the outputs to be considered up to date.
#   required_inputs : upstream modules whose outputs MUST be present before this module can run;
#                     a hard error is raised if they're missing when the module needs to run.
#                     (count / spatial_positioning are intentionally empty: they can also run
#                     from FASTQs supplied directly via --fastqs, and cellbender is optional for
#                     spatial_analysis so it is an mtime input but not a required one.)
_MODULE_IO = {
    'mkfastq': {
        'dir': 'mkfastq',
        'outputs': None,
        'mtime_inputs': [],
        'required_inputs': [],
    },
    'count': {
        'dir': 'count',
        'outputs': ['filtered_feature_bc_matrix.h5', 'molecule_info.h5', 'metrics_summary.csv'],
        'mtime_inputs': ['mkfastq'],
        'required_inputs': [],
    },
    'cellbender': {
        'dir': 'cellbender',
        'outputs': ['cellbender_output_filtered.h5'],
        'mtime_inputs': ['count'],
        'required_inputs': ['count'],
    },
    'spatial_positioning': {
        'dir': 'spatial_barcodes',
        'outputs': ['SBcounts.h5'],
        'mtime_inputs': ['mkfastq'],
        'required_inputs': [],
    },
    'spatial_analysis': {
        'dir': 'spatial_analysis',
        'outputs': ['spatial_metadata.json', 'summary.pdf', 'cb_whitelist.txt',
                    'coords.csv', 'matrix.csv.gz', 'seurat.qs'],
        'mtime_inputs': ['count', 'cellbender', 'spatial_positioning'],
        'required_inputs': ['count', 'spatial_positioning'],
    },
}


def _module_sample_dir(module: str, sample_name: str) -> Path:
    """Path to a single sample's output folder for a module."""
    return OUTPUT_PATH / _MODULE_IO[module]['dir'] / sample_name


def _outputs_complete(module: str, sample_name: str) -> bool:
    """True if every expected output file for (module, sample) is present on disk."""
    sample_dir = _module_sample_dir(module, sample_name)
    if not sample_dir.is_dir():
        return False
    outputs = _MODULE_IO[module]['outputs']
    if outputs is None:
        # mkfastq: require at least one R1 and one R2 demultiplexed FASTQ in the sample folder
        fastqs = [f.name for f in sample_dir.glob('*.fastq.gz')]
        return any('_R1_' in n for n in fastqs) and any('_R2_' in n for n in fastqs)
    return all((sample_dir / f).exists() for f in outputs)


def _output_files(module: str, sample_name: str) -> list[Path]:
    """Existing expected output-file Paths for (module, sample); the FASTQs for mkfastq."""
    sample_dir = _module_sample_dir(module, sample_name)
    if not sample_dir.is_dir():
        return []
    outputs = _MODULE_IO[module]['outputs']
    if outputs is None:
        return list(sample_dir.glob('*.fastq.gz'))
    return [sample_dir / f for f in outputs if (sample_dir / f).exists()]


def _input_files(module: str, input_module: str, sample_name: str) -> list[Path]:
    """Existing files from an upstream module that feed `module` for this sample."""
    # with --fastqs, mkfastq never runs and so has no outputs to compare mtimes against; the
    # externally supplied FASTQs are the real upstream input, so use those instead. Without this,
    # count/spatial_positioning outputs would always look 'fresh' and a newer set of input FASTQs
    # would be silently ignored on re-run.
    if input_module == 'mkfastq' and FASTQ_INPUT is not None:
        return find_sample_fastqs(FASTQ_INPUT, sample_name, 'SB' if module == 'spatial_positioning' else 'GEX')

    # spatial_positioning consumes the spatial-barcode library, which mkfastq demultiplexes into
    # a separate <sample>_sb/ folder rather than the sample's own <sample>/ folder
    if module == 'spatial_positioning' and input_module == 'mkfastq':
        return _output_files('mkfastq', f'{sample_name}_sb')
    return _output_files(input_module, sample_name)


def _sample_output_state(module: str, sample_name: str) -> str:
    """
    Classify a sample's outputs for a module as one of:
     - 'missing' : one or more expected output files are absent
     - 'stale'   : outputs are all present but an existing upstream input is newer than them
     - 'fresh'   : outputs are all present and no existing input is newer (nothing to redo)
    A module with no comparable upstream inputs present is treated as 'fresh' when its outputs
    exist, so already-computed results are never needlessly overwritten.
    """
    if not _outputs_complete(module, sample_name):
        return 'missing'

    output_mtimes = [f.stat().st_mtime for f in _output_files(module, sample_name)]
    if not output_mtimes:
        return 'fresh'

    input_mtimes = []
    for input_module in _MODULE_IO[module]['mtime_inputs']:
        input_mtimes += [f.stat().st_mtime for f in _input_files(module, input_module, sample_name)]
    if not input_mtimes:
        return 'fresh'

    # outputs are up to date only if the newest input predates the oldest output
    return 'fresh' if max(input_mtimes) < min(output_mtimes) else 'stale'


def _relevant_sample_names(module: str, metadata_df: pd.DataFrame) -> list[str]:
    """Sample names a module produces per-sample outputs for."""
    names = metadata_df['Sample Name'].astype(str).tolist()
    # spatial-barcode (_sb) libraries have no independent per-sample outputs for the modules that
    # run after mkfastq; skip them if they ever appear as their own metadata rows. Match the exact
    # `_sb` suffix, not a substring: a real sample whose name merely contains "_sb" (e.g.
    # "lung_sb_rep1") must NOT be dropped, or its stage would be silently skipped.
    if module in ('count', 'cellbender', 'spatial_positioning', 'spatial_analysis'):
        names = [n for n in names if not n.endswith('_sb')]
    return names


def row_declares(sample_row, column: str) -> bool:
    """
    True if the metadata row carries a non-empty value in `column` (optional columns may be absent).

    Public because `create_samplesheet` decides which libraries to write from the same columns this
    file's completeness check reads them back with. Both have to agree on what a declared library is:
    a blank primary index means the samplesheet writer emits no primary library, and so the
    completeness check must not require its read pair (see `_library_complete`). One predicate rather
    than two keeps a blank, a NaN and a whitespace-only cell from counting as declared in one place and
    not the other.
    """
    return (column in sample_row.index
            and pd.notna(sample_row[column])
            and str(sample_row[column]).strip() != '')


def _is_flex(sample_row) -> bool:
    """True if the sample uses Flex chemistry (matched exactly, as everywhere else in the pipeline)."""
    return 'Chemistry' in sample_row.index and str(sample_row['Chemistry']).strip() == 'Flex'


def _dir_has_read_pair(
    sample_dir: Path,
    require_token: str | None = None,
    exclude_token: str | None = None
) -> bool:
    """
    True if `sample_dir` holds at least one R1 and one R2 FASTQ, optionally filtered by a filename
    token: only files containing `require_token` are considered, and files containing
    `exclude_token` are ignored.
    """
    if not sample_dir.is_dir():
        return False
    names = [f.name for f in sample_dir.glob('*.fastq.gz')]
    if require_token is not None:
        names = [n for n in names if require_token in n]
    if exclude_token is not None:
        names = [n for n in names if exclude_token not in n]
    return any('_R1_' in n for n in names) and any('_R2_' in n for n in names)


def _library_complete(
    library_dir: Path,
    sample_row,
    index_column: str,
    merge_column: str,
    split_token: str
) -> bool:
    """
    True if every BCL contributing one of a sample's libraries has been demultiplexed into
    `library_dir`.

    Which read pairs to expect is stated by the metadata rather than guessed from the directory: a
    primary `index_column` means mkfastq produces a pair whose filenames lack `split_token`, and a
    `merge_column` means it produces one carrying it. A sample may declare either or both -- a
    library sequenced entirely on a merged-in run has a blank primary index (see
    `Merge ... From BCL` / `Add ... Index` in the metadata reference) -- so demanding a pair the
    metadata never asked for would leave the sample permanently incomplete, re-running mkfastq on
    every invocation with no amount of demultiplexing able to satisfy it.

    A row declaring neither has no library of this kind to demultiplex at all, and is reported as
    nothing-to-do rather than never-satisfiable. Nothing downstream would name it: no module lists
    `mkfastq` in its `required_inputs`, and run_count's own guard only checks that the mkfastq tree is
    non-empty overall, not that a given sample has reads in it. So the case is rejected up front instead
    -- `create_samplesheet` refuses a blank `RNA Index` with no `Merge RNA From BCL` to explain it, and
    warns for the spatial equivalent, both against the metadata that caused it. By the time a row
    reaches here it therefore declares at least one library, and the two branches above cover it.

    Inputs:
     - library_dir:  mkfastq output folder for this library (mkfastq/<sample>[_sb]/)
     - sample_row:   the sample's metadata row
     - index_column: metadata column naming the primary-BCL index ('RNA Index' / 'SB Index')
     - merge_column: metadata column naming the merged-in BCL
     - split_token:  filename token create_samplesheet gives the merged-in library
    Output:
     - True if every declared contributing BCL's read pair is present
    """

    if row_declares(sample_row, index_column) and not _dir_has_read_pair(library_dir, exclude_token=split_token):
        return False
    if row_declares(sample_row, merge_column) and not _dir_has_read_pair(library_dir, require_token=split_token):
        return False
    return True


def _mkfastq_sample_complete(sample_row) -> bool:
    """
    True only if every BCL that contributes to this sample has been demultiplexed.

    A sample's gene-expression FASTQs live in mkfastq/<sample>/ and its spatial-barcode FASTQs in
    mkfastq/<sample>_sb/. When a sample merges reads from an additional BCL, create_samplesheet
    names that extra library <sample>_split_rna / <sample>_split_sb, so the merged FASTQs carry a
    '_split_rna' / '_split_sb' token in their filenames while the primary BCL's do not. Requiring
    each read pair the metadata declares, in both the GEX and spatial folders, ensures a
    split/multi-BCL sample is not treated as done as soon as just one BCL's FASTQs appear -- see
    _library_complete for why the primary pair is required only when a primary index is declared.
    """
    name = str(sample_row['Sample Name'])
    rna_dir = OUTPUT_PATH / 'mkfastq' / name
    sb_dir = OUTPUT_PATH / 'mkfastq' / f'{name}_sb'

    # gene-expression library
    if not _library_complete(rna_dir, sample_row, 'RNA Index', 'Merge RNA From BCL', '_split_rna'):
        return False

    # spatial-barcode library: skipped for Flex chemistry, whose spatial barcodes are supplied as
    # external FASTQs (flex_spatial_R1/R2_path -> Trekker demux) rather than demultiplexed by
    # mkfastq, so mkfastq/<sample>_sb/ is never produced for Flex samples
    if not _is_flex(sample_row):
        if not _library_complete(sb_dir, sample_row, 'SB Index', 'Merge Spatial From BCL', '_split_sb'):
            return False

    return True


def need_run_module(
    module: str,
    metadata_df: pd.DataFrame
) -> bool:
    """
    Decide whether a pipeline module needs to run, to avoid overwriting up-to-date outputs.

    A module is skipped (returns False) only when, for every sample it is responsible for, all of
    its expected output files are present AND no existing upstream input file is newer than those
    outputs. Otherwise it runs (returns True); if it must run but a required upstream module's
    outputs are missing, a clear error is logged and the process exits.

    Inputs:
     - module:      pipeline module name (a key of _MODULE_IO)
     - metadata_df: DataFrame of sample metadata (must contain a 'Sample Name' column)
    Output:
     - True if the module should run, False if it can be skipped
    """

    if module not in _MODULE_IO:
        log_write(f'[INTERNAL ERROR]: unrecognized module passed to need_run_module(): `{module}`')
        log_write("Troubleshooting:")
        log_write(" • This is a bug in slidr, not a problem with your setup or data")
        log_write(f" • The known modules are: {', '.join(_MODULE_IO)}")
        log_write(" • Please report it at https://github.com/kabanovskyd/slidr/issues, quoting the module name above")
        sys.exit(1)

    # _MODULE_IO models the standard (non-Flex) on-disk layout. Flex chemistry uses a different
    # layout (cellranger multi writes count outputs under count/flex/..., and the spatial side goes
    # through the Takara/Trekker path), none of which this skip-logic models. Rather than
    # mis-classify Flex outputs as missing (which would both re-run needlessly and hard-exit on the
    # required-inputs check below), always run the requested stage for Flex and let it handle its
    # own I/O.
    if 'Chemistry' in metadata_df.columns and (metadata_df['Chemistry'].astype(str).str.strip() == 'Flex').any():
        return True

    # classify each sample's outputs once
    if module == 'mkfastq':
        # mkfastq has no upstream file outputs to compare mtimes against, so a sample is simply
        # complete ('fresh') or not ('missing'). A split/multi-BCL sample counts as complete only
        # when every contributing BCL's FASTQs are present (see _mkfastq_sample_complete), so a
        # failed/absent merge BCL correctly forces a re-run rather than being silently skipped.
        states = {}
        for _, row in metadata_df.iterrows():
            name = str(row['Sample Name'])
            states[name] = 'fresh' if _mkfastq_sample_complete(row) else 'missing'
    else:
        samples = _relevant_sample_names(module, metadata_df)
        states = {sample: _sample_output_state(module, sample) for sample in samples}

    missing = [s for s, state in states.items() if state == 'missing']
    stale = [s for s, state in states.items() if state == 'stale']

    # every sample is present and up to date -> nothing to do
    if not missing and not stale:
        return False

    # the module must run; verify its hard-required upstream outputs exist first
    required_inputs = _MODULE_IO[module]['required_inputs']
    if required_inputs:
        missing_inputs = {}
        for sample in samples:
            absent = [im for im in required_inputs if not _outputs_complete(im, sample)]
            if absent:
                missing_inputs[sample] = absent
        if missing_inputs:
            log_write(f"[ERROR]: cannot run {module}: required upstream outputs are missing.")
            for sample, mods in missing_inputs.items():
                log_write(f"  • {sample}: missing {', '.join(mods)} output(s)")

            # name the exact flags that regenerate what is missing, rather than "the upstream stage(s)"
            stage_flags = {
                'mkfastq': '--mkfastq',
                'count': '--count',
                'cellbender': '--cellbender',
                'spatial_positioning': '--spatial-count',
                'spatial_analysis': '--spatial-analysis',
            }
            absent_modules = sorted({m for mods in missing_inputs.values() for m in mods})
            flags = ' '.join(stage_flags[m] for m in absent_modules if m in stage_flags)

            log_write("Troubleshooting:")
            if flags:
                log_write(f" • Generate them by re-running the upstream stage(s): `./slidr --bcl {cfg['bcl_id']} {flags}`")
            log_write(" • Or run the pipeline end to end with --run-all")
            log_write(f" • If the outputs do exist, check they are complete -- every expected file must be present under {OUTPUT_PATH}")
            for absent in absent_modules:
                expected = _MODULE_IO.get(absent, {}).get('outputs')
                if expected:
                    log_write(f"   - {absent} must contain: {', '.join(expected)}")
            log_write(" • If the reads were processed elsewhere, pass the FASTQ directory with --fastqs so mkfastq can be skipped")
            sys.exit(1)

    # log why the module is running
    if missing:
        log_ts(f"{module} needs to run — no output yet for {len(missing)} sample(s)",
               [f"• {sample}" for sample in missing])
    if stale:
        log_ts(f"{module} needs to run — inputs newer than outputs for {len(stale)} sample(s)",
               [f"• {sample}" for sample in stale])

    return True


def format_multi_samplesheet(
    gene_expression: dict,
    libraries: list[dict],
    samples: list[dict],
    outdir: Path | str,
):
    """
    Format a samplesheet for running cellranger multi based on provided inputs

    Inputs:
    - gene_expression: 
    - libraries:
    - samples:                list of samples to analyze
    - outdir:                 path to output samplesheet directory

    Outputs:
    - outpath:                path to formatted samplesheet

    """

    # aux function for writing individual samplesheet sections as CSV
    def write_section(writer, header, data):
        # write header
        writer.writerow([header])

        # write dictionary values
        if isinstance(data, dict):
            for k, v in data.items():
                writer.writerow([k, v])

        # write list values
        elif isinstance(data, list):
            writer.writerow(list(data[0].keys()))
            for row in data:
                writer.writerow(list(row.values()))

        # blank line between sections
        writer.writerow([])

    # save samplesheet at specified path
    outpath = Path(outdir) / "multi_samplesheet.csv"

    # write the samplesheet to destination
    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f)
        write_section(writer, "[gene-expression]", gene_expression)
        write_section(writer, "[libraries]", libraries)
        write_section(writer, "[samples]", samples)

    return outpath


def parse_cellranger_html(
    html_path: str,
    plateau_buffer: int = 5_000
) -> dict:
    """
    Parse a cellranger count web_summary.html and return cellbender parameters
 
    Inputs:
    - html_path:               path to the cellranger ``web_summary.html`` file.
    - plateau_buffer:          number of barcodes to add beyond the cliff base to place the cutoff
                               safely inside the empty-droplet plateau
    Outputs:
    - A dictionary of key cellbender metrics:
        - expected_cells
        - total_droplets_included
        - cliff_base_rank
        - cellranger_estimated_cells
    """

    def unparseable(detail: str) -> None:
        """
        Report that web_summary.html could not be mined for cellbender parameters, and exit.

        Every failure in this function has the same cause (the report is truncated, or was written
        by a cellranger version whose summary JSON is laid out differently) and the same two ways
        out (regenerate the report, or set the parameters by hand), so they share one message
        instead of repeating it per parse step.
        """

        log_write(f"[ERROR]: Could not parse {html_path} to load cellbender parameters automatically: {detail}")
        log_write("Troubleshooting:")
        log_write(" • Set `workflow.cellbender_total_droplets` and `workflow.cellbender_estimated_cells` in the config file to skip auto-detection entirely")
        log_write(" • Or re-run the count stage with --count --force to regenerate the report, then re-run cellbender")
        log_write(" • A truncated report is the usual cause; check the file opens in a browser and shows a barcode-rank plot")
        log_write(" • Cellbender's own guidance on choosing these values: https://cellbender.readthedocs.io/en/latest/usage/index.html")
        sys.exit(1)

    # read the provided HTML file and search for key string
    html = Path(html_path).read_text(errors="ignore")
    match = re.search(r"const data = ", html)
    if not match:
        unparseable("`const data = ` not found in file")

    # parse the JSON elements in the file
    decoder = json.JSONDecoder()
    try:
        data, _ = decoder.raw_decode(html[match.end():])
    except json.JSONDecodeError as exc:
        unparseable(str(exc))

    # extract the number of cells estimated by cellranger
    try:
        summary_tab = data["summary"]["summary_tab"]
        cells_metric = summary_tab["filtered_bcs_transcriptome_union"]["metric"]
        cellranger_estimated_cells = int(cells_metric.replace(",", ""))
    except (KeyError, TypeError, ValueError) as exc:
        unparseable(f"could not read the estimated cell count: {exc}")

    # extract the knee plot data
    try:
        knee_plot = data["summary"]["summary_tab"]["cells"]["barcode_knee_plot"]
        plot_traces = knee_plot["data"]
    except (KeyError, TypeError) as exc:
        unparseable(f"could not read the barcode rank (knee) plot: {exc}")
 
    cells_100_max_ranks = []
    cells_all_max_ranks = []
 
    # iterate over the traces in the knee plot
    for trace in plot_traces:
        if trace.get("name") != "Cells":
            continue
        x = trace.get("x", [])
        if not x:
            continue
 
        # add the maximum rank for the "Cells" trace
        cells_all_max_ranks.append(max(x))
 
        # user regex to find a trace with "100%" label
        text_list = trace.get("text", [])
        first_text = text_list[0] if text_list else ""
        pct_match = re.match(r"(\d+)%\s*Cells", first_text)
        if pct_match and int(pct_match.group(1)) == 100:
            # save the maximum rank from the 100% trace
            cells_100_max_ranks.append(max(x))
 
    # exit if cannot extract barcode ranks from traces
    if not cells_100_max_ranks:
        unparseable("no '100% Cells' trace in the barcode rank plot, so expected_cells cannot be determined")
    if not cells_all_max_ranks:
        unparseable("no 'Cells' trace in the barcode rank plot, so cliff_base_rank cannot be determined")
 
    # extract the maximum ranks from each and compute the total droplets value
    expected_cells = max(cells_100_max_ranks)
    cliff_base_rank = max(cells_all_max_ranks)
    total_droplets_included = round((cliff_base_rank + plateau_buffer) / 1_000) * 1_000
 
    return {
        "expected_cells": expected_cells,
        "total_droplets_included": total_droplets_included,
        "cliff_base_rank": cliff_base_rank,
        "cellranger_estimated_cells": cellranger_estimated_cells,
    }


def ensure_conda_env(
    env_name: str,
    environment_yml: str | None = 'envs/conda.yml'
) -> Path:
    """
    Ensure a named conda environment exists, creating it from a YAML file if absent.

    An environment is located by name across every directory mamba knows about, so one installed
    outside the default location (a shared, prebuilt env on a cluster, say) is picked up as long as
    its parent directory is listed under `envs_dirs` in ~/.condarc.

    Inputs:
     - env_name:              name of the conda environment to locate or create
     - environment_yml:       environment definition file used if the environment does not already
                              exist, resolved relative to the repository root (default:
                              envs/conda.yml). Pass None for an environment slidr cannot build
                              itself, in which case a missing environment is reported rather than
                              created.
    Output:
     - Path to the bin/ directory of the resolved conda environment
    """

    # Resolve the definition file against the repository root rather than the current directory.
    # ./slidr happens to cd there first, but a run launched any other way (a bare `uv run
    # workflow/main.py` from elsewhere, a Slurm job whose --chdir differs) would otherwise look for
    # envs/conda.yml relative to wherever it started and report it missing.
    env_file = None if environment_yml is None else Path(ROOT_PATH) / environment_yml

    # exit if a mamba installation is not found on the system
    if shutil.which('mamba') is None:
        log_write(f"[ERROR]: 'mamba' not found on PATH, so the `{env_name}` environment cannot be created or located")
        log_write("Troubleshooting:")
        log_write(" • Launch the pipeline through ./slidr, which installs Miniforge automatically if it is missing")
        log_write(" • Or activate an existing mamba/miniforge installation first (e.g. `source ~/miniforge3/bin/activate`)")
        log_write(" • Install Miniforge manually: https://github.com/conda-forge/miniforge#install")
        log_write(" • On a cluster, check whether conda/mamba needs loading first (e.g. `module load miniforge`)")
        sys.exit(1)

    # list existing envs and see if the target one exists. The result is parsed defensively: a failed
    # `mamba env list` would otherwise hand empty/non-JSON output to json.loads and surface as a bare
    # JSONDecodeError that says nothing about what actually went wrong.
    def list_envs() -> list[str]:
        result = subprocess.run(
            ['mamba', 'env', 'list', '--json'],
            cwd=SCRIPT_PATH,
            capture_output=True,
            text=True
        )
        try:
            return json.loads(result.stdout)['envs']
        except (json.JSONDecodeError, KeyError, TypeError) as exc:
            log_write(f"[ERROR]: could not read the list of conda environments from mamba: {exc}")
            log_write("Troubleshooting:")
            log_write(" • Run `mamba env list --json` by hand -- its own error names the cause")
            detail = (result.stderr or '').strip().splitlines()
            for line in detail[-3:]:
                log_write(f"   {line}")
            log_write(" • A broken ~/.condarc (bad YAML, or an `envs_dirs` entry that cannot be read) is a common cause")
            sys.exit(1)

    envs = list_envs()
    match = next((e for e in envs if Path(e).name == env_name), None)

    # if not found, create an env from environment file
    if match is None:
        # Some environments are not slidr's to build -- the Trekker/Flex env has no definition in
        # this repository and is expected to be provided by the site. Say that plainly, and point at
        # envs_dirs, rather than blaming a definition file the user was never given.
        if env_file is None:
            log_write(f"[ERROR]: the `{env_name}` conda environment was not found, and slidr has no definition file to build it from")
            log_write("Troubleshooting:")
            log_write(f" • This environment is not created by slidr; it has to exist before the run starts")
            log_write(" • If it lives outside your default conda location, add that directory to `envs_dirs` in ~/.condarc:")
            log_write("     envs_dirs:")
            log_write("       - /path/to/shared/envs")
            log_write(" • Check what mamba can currently see: `mamba env list`")
            log_write(f" • Environments found: {', '.join(sorted(Path(e).name for e in envs)) or '(none)'}")
            sys.exit(1)

        if not env_file.is_file():
            log_write(f'[ERROR]: the `{env_name}` conda environment does not exist and its definition file {env_file} is not a valid file')
            log_write("Troubleshooting:")
            log_write(f" • This file is checked into the repository; restore it with `git checkout {environment_yml}`")
            log_write(f" • Or create the environment yourself: `mamba create -n {env_name} ...`")
            log_write(f" • Environments mamba can currently see: {', '.join(sorted(Path(e).name for e in envs)) or '(none)'}")
            sys.exit(1)

        subprocess.run(['mamba', 'env', 'create', '-f', str(env_file)], check=True)

        # re-query to get the env path after creation
        envs = list_envs()
        match = next((e for e in envs if Path(e).name == env_name), None)
        if match is None:
            log_write(f"[ERROR]: the `{env_name}` conda environment is still not present after creating it from {env_file}")
            log_write("Troubleshooting:")
            log_write(f" • Check the `name:` field inside {env_file} is `{env_name}` -- mamba names the environment from the file, not from slidr")
            log_write(f" • Create it by hand to see the full error: `mamba env create -f {env_file}`")
            log_write(f" • If a broken half-built environment is in the way, remove it first: `mamba env remove -n {env_name}`")
            sys.exit(1)

    return Path(match) / 'bin'


def slack_message(message: str) -> None:
    """
    Send a Slack alert to the user, tagged with this run's BCL ID.

    The BCL ID is prepended here rather than at each call site so every alert carries it: a user
    running several BCLs at once (or the same BCL on two machines) otherwise receives a stream of
    identical-looking "started"/"error" DMs with nothing to tell them apart.

    Inputs:
     - message:            content of the message to send
    """

    # instantiate the WebClient object by the specific Slack token
    try:
        client = WebClient(token=SLACK_TOKEN)
    except Exception as exc:
        log_write(f"[WARNING]: cannot set up Slack WebClient: {exc}")
        log_write("Disabling Slack alerts for the duration of this run.")
        log_write(" • Check `settings.slack_token` is a valid bot token (it should start with `xoxb-`)")
        log_write(" • Set `settings.alerts` to false to silence this warning")
        return

    # extract user's email from the metadata sheet
    if not SUMMARY_PATH.is_file():
        log_write(f"[WARNING]: cannot read {SUMMARY_PATH} as it does not exist. User email could not be extracted, so Slack alerts have been disabled.")
        log_write(" • This file is written once the metadata loads, so an alert raised before that point cannot be delivered")
        log_write(" • The run itself is unaffected; check the runtime log for the underlying error")
        return
    metadata_df = pd.read_csv(SUMMARY_PATH)
    email=metadata_df['Email'].tolist()[0]

    # look up user ID by email
    try:
        result = client.users_lookupByEmail(email=email)
    except Exception as exc:
        log_write(f"[WARNING]: could not look up user email {email}: {exc}")
        log_write(" • Check the `Email` metadata column holds the address the user's Slack account is registered under")
        log_write(" • The Slack bot needs the `users:read.email` scope to resolve an address to a user")
        return
    user_id = result["user"]["id"]

    try:
        # send a DM to the user through the client object, with the BCL ID as a bold first line so
        # the run is identifiable from the DM preview without opening the message
        client.chat_postMessage(
            channel=user_id,
            text=f"*{BCL_ID}*\n{message}"
        )
    except SlackApiError as e:
        print(f"[WARNING]: failed to send a Slack message: {e.response['error']}")
        print(" • The bot needs the `chat:write` scope, and must not be restricted from DMing this user")
        print(" • Set `settings.alerts` to false to disable Slack alerts for future runs")


def create_tmp_dir(module: str) -> Path:
    tmp_dir = TMP_PATH / module / 'run_0'
    if tmp_dir.is_dir():
        suffix = 1
        while tmp_dir.is_dir():
            tmp_dir = TMP_PATH / module / f'run_{suffix}'
            suffix += 1
    tmp_dir.mkdir(parents=True)

    tmp_dir_size = sum(f.stat().st_size for f in TMP_PATH.rglob("*") if f.is_file()) / 1024**3
    if tmp_dir_size > WARN_THRESHOLD_GB:
        log_write(f'[WARNING]: the size of your {TMP_PATH} directory is {tmp_dir_size:.1f} GB. Please consider deleting it after this run is completed.')
        log_write(' • Each stage writes into a fresh run_N/ subdirectory, so scratch from earlier runs accumulates here')
        log_write(f' • Once this run has finished successfully the whole directory is safe to remove: `rm -rf {TMP_PATH}`')
    if tmp_dir_size > CRITICAL_THRESHOLD_GB:
        log_write(f'[ERROR]: the size of your {TMP_PATH} directory has reached the critical threshold of {CRITICAL_THRESHOLD_GB} GB ({tmp_dir_size:.1f} GB)')
        log_write("Troubleshooting:")
        log_write(f" • Delete the scratch directory and re-run: `rm -rf {TMP_PATH}`")
        log_write(" • Nothing in it is a pipeline output; outputs live under the `output/` directory alongside it")
        log_write(f" • Check the free space on this filesystem: `df -h {TMP_PATH}`")
        log_write(" • Point `paths.output_path` at a larger scratch filesystem if runs of this size are routine")
        sys.exit(1)

    return tmp_dir


def validate_bcl_dir(bcl_path: Path | str) -> bool:
    """
    Scan a directory and check whether it conforms to the standard Illumina BCL
    run-folder layout required by cellranger mkfastq, rather than e.g. an already
    demultiplexed FASTQ directory.

    Checks for:
     - RunInfo.xml, parseable and containing <Flowcell> and <Reads> elements
     - RunParameters.xml (NovaSeq/NextSeq) or runParameters.xml (HiSeq/MiSeq)
     - a Data/Intensities/BaseCalls directory containing basecall files
       (*.cbcl for NovaSeq/NextSeq, or *.bcl/*.bcl.gz for older instruments)

    Inputs:
     - bcl_path: path to the candidate BCL run directory
    Output:
     - True if the directory conforms to the expected BCL run-folder structure,
       False otherwise (details of what's missing are written to the runtime log)
    """

    bcl_path = Path(bcl_path)
    issues = []

    if not bcl_path.is_dir():
        return False

    # RunInfo.xml is required to extract the flowcell ID and read structure for mkfastq
    run_info = bcl_path / "RunInfo.xml"
    if not run_info.is_file():
        issues.append(f"missing RunInfo.xml at {run_info}")
    else:
        try:
            root = ET.parse(run_info).getroot()
            if root.find(".//Flowcell") is None:
                issues.append(f"{run_info} does not contain a <Flowcell> element")
            if root.find(".//Reads") is None:
                issues.append(f"{run_info} does not contain a <Reads> element")
        except ET.ParseError as exc:
            issues.append(f"could not parse {run_info}: {exc}")

    # RunParameters.xml naming differs by instrument (NovaSeq/NextSeq vs. HiSeq/MiSeq)
    if not any((bcl_path / name).is_file() for name in ("RunParameters.xml", "runParameters.xml")):
        issues.append(f"missing RunParameters.xml/runParameters.xml in {bcl_path}")

    # basecall files must actually be present under Data/Intensities/BaseCalls
    basecalls = bcl_path / "Data" / "Intensities" / "BaseCalls"
    if not basecalls.is_dir():
        issues.append(f"missing Data/Intensities/BaseCalls directory in {bcl_path}")
    elif not any(basecalls.rglob("*.cbcl")) and not any(basecalls.rglob("*.bcl*")):
        issues.append(f"no .cbcl or .bcl/.bcl.gz basecall files found under {basecalls}")

    if issues:
        log_write(f"[WARNING]: {bcl_path} does not conform to the expected BCL run directory structure:", terminal=False)
        for issue in issues:
            log_write(f"  • {issue}", terminal=False)
        return False

    return True


def stage_from_gcs(
    source: str,
    dest: Path | str,
    recursive: bool = False,
    required: bool = True,
    description: str | None = None
) -> bool:
    """
    Copy an object or prefix out of Google Cloud Storage onto local disk with `gcloud storage cp`.

    This is the single place remote inputs are staged, so every resource the pipeline pulls down
    (reference genome, puck maps, raw barcodes) reports failures the same way and, unlike the
    per-call-site subprocess blocks this replaces, surfaces gcloud's own stderr instead of sending
    it to /dev/null -- without which a staging failure is undiagnosable from the log.

    Inputs:
     - source:      GCS source URI (gs://...)
     - dest:        local destination file or directory
     - recursive:   pass -r, for copying a prefix rather than a single object
     - required:    exit on failure; when False a failure is a warning and the caller decides what
                    to do (a missing puck map, for instance, is regenerated from raw barcodes)
     - description: short human-readable name of what is being staged, used in log messages
    Output:
     - True if the copy succeeded, False if it failed and `required` was False
    """

    label = description or source
    cmd = ['gcloud', 'storage', 'cp'] + (['-r'] if recursive else []) + [source, str(dest)]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        log_write("[ERROR]: `gcloud` was not found on PATH, so inputs cannot be staged from GCS")
        log_write("Troubleshooting:")
        log_write(" • Install the Google Cloud CLI: https://cloud.google.com/sdk/docs/install")
        log_write(" • On a cluster, check whether it needs loading first (e.g. `module load google-cloud-sdk`)")
        log_write(" • Drop --stage-gcs to read these inputs off the local filesystem instead")
        sys.exit(1)

    if proc.returncode == 0:
        return True

    # gcloud writes the actionable part of its diagnostics to stderr; keep the tail of it
    detail = (proc.stderr or proc.stdout or '').strip().splitlines()
    detail = detail[-3:] if detail else ['(no output)']

    severity = 'ERROR' if required else 'WARNING'
    log_write(f"[{severity}]: could not stage {label} from GCS (`{' '.join(cmd)}` exited {proc.returncode})")
    for line in detail:
        log_write(f"  {line}")

    if not required:
        return False

    log_write("Troubleshooting:")
    log_write(f" • Check that {source} exists: `gcloud storage ls {source}`")
    log_write(" • Check that you are authenticated: `gcloud auth print-access-token`")
    log_write(" • Check that the active account can read the bucket")
    sys.exit(1)


def resolve_bcl_dir(bcl_id: str) -> Path:
    """
    Resolve the Illumina run folder for a BCL ID from the run's local input directory.

    That directory is `paths.input_path` for a local run, and the `settings.gcs_download_dest` the
    bucket was staged into for a `--stage-gcs`/`--gcp` one -- config.py re-points INPUT_PATH at the
    latter, so this function (and everything else downstream) sees one local root either way.

    It is documented as the root directory holding BCL run folders, so the run folder is normally
    input_path/<BCL_ID>; a path that already points at the run folder itself (its basename is this
    run's BCL ID) is accepted, in which case the root holding every run folder is its parent.
    That second case matters for split-BCL runs: the extra run folders a sample merges reads from sit
    beside the primary one, not inside it. This mirrors how run_mkfastq derives each --run argument,
    so validating the directory returned here checks what mkfastq will actually be handed -- unlike
    validating `input_path` itself, which for the documented layout is a directory *of* run folders
    and so never satisfies the BCL run-folder schema.

    Inputs:
     - bcl_id: BCL run ID -- this run's (args.bcl) or one named by a `Merge ... From BCL` column
    Output:
     - path to the run folder for that BCL ID
    """

    input_path = Path(INPUT_PATH)
    root = input_path.parent if input_path.name == BCL_ID else input_path
    return root / bcl_id


def merge_bcls(metadata_df: pd.DataFrame) -> list[str]:
    """
    The additional BCL run IDs a split-BCL run merges reads from.

    A library is sometimes sequenced across more than one run: the `Merge RNA From BCL` /
    `Merge Spatial From BCL` columns name the extra run folder a sample's gene-expression or spatial
    reads are topped up from. Those folders are inputs to this run exactly as the primary BCL is --
    mkfastq demultiplexes each of them in turn -- so both the demultiplexing loop and GCS staging
    need the same list, derived in one place.

    Inputs:
     - metadata_df: the run's metadata table
    Output:
     - de-duplicated merge-from BCL IDs in metadata order, excluding this run's own BCL and any
       blank/NaN entry
    """

    bcls = []
    for column in ('Merge RNA From BCL', 'Merge Spatial From BCL'):
        if column not in metadata_df:
            continue
        for value in metadata_df[column]:
            if value is None or pd.isna(value):
                continue
            bcl = str(value).strip()
            if bcl and bcl != BCL_ID and bcl not in bcls:
                bcls.append(bcl)
    return bcls


def declared_lanes(sample_row, column: str) -> list[str]:
    """
    Lane tokens ('L001', 'L002', ...) a sample declares in a metadata lane column, zero-padded to
    the form cellranger/bcl2fastq embed in FASTQ filenames.

    Inputs:
     - sample_row: a metadata DataFrame row
     - column:     lane column to read ('Lane' for gene expression, 'SB Lane' for spatial barcodes)
    Output:
     - list of lane tokens, or an empty list meaning "all lanes" (a blank value, a '*' wildcard, or
       an unparseable entry), in which case callers must not filter by lane
    """

    if column not in sample_row.index or pd.isna(sample_row[column]):
        return []

    lanes = []
    for lane in str(sample_row[column]).split(','):
        lane = lane.strip()
        if not lane or lane == '*':
            return []
        try:
            lanes.append(f"L{int(lane):03d}")
        except ValueError:
            log_write(f"[WARNING]: unrecognized lane value '{lane}' in the `{column}` metadata column for "
                      f"{sample_row.get('Sample Name', 'sample')}; matching FASTQs from all lanes instead")
            log_write(" • Lane values must be plain integers, comma-separated for several lanes (e.g. `1,2`), or `*` for all lanes")
            log_write(" • Leaving the column blank also means all lanes, so set it only if this sample really is lane-restricted")
            return []

    return lanes


# Separators accepted between probe barcode IDs in the `Flex Probe Barcode IDs` metadata column.
# A comma is the documented form, but a pipe is common in exports from instrument software and in
# spreadsheets where a comma would otherwise need quoting. Both are accepted so "BC001|BC002" is not
# silently taken as one barcode literally named "BC001|BC002" -- which cellranger multi then rejects
# with an opaque error about an unknown probe barcode.
PROBE_BARCODE_SEPARATORS = ',|'


def split_probe_barcodes(value) -> list[str]:
    """
    Split a `Flex Probe Barcode IDs` metadata value into individual probe barcode IDs.

    Accepts ',' or '|' as the separator (including a mixture of the two), drops empty fields left by
    a trailing or doubled separator, and strips surrounding whitespace.

    Inputs:
     - value: raw metadata value for the column
    Output:
     - list of probe barcode IDs in the order written; empty when the value holds none
    """

    # a blank cell arrives as NaN from pandas rather than as an empty string
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return []

    pattern = rf'[{re.escape(PROBE_BARCODE_SEPARATORS)}]'
    return [bc.strip() for bc in re.split(pattern, str(value)) if bc.strip()]


def library_matches_sample(library: str, sample_name: str, mode: str) -> bool:
    """
    Whether a FASTQ library prefix names the given sample's given library, per the convention
    documented on FASTQ_MODES.

    The sample name is anchored to the start of the library and must be followed by the mode token
    (not merely be a substring of the library), so a sample never claims the FASTQs of another
    sample whose name it is a prefix of -- 'tumor' must not pick up 'tumor_ln_GEX'. This is the
    same exact-match discipline run_mkfastq's route_map uses when routing demultiplexed FASTQs.
    Only the mode token is matched case-insensitively; the sample name is matched as written.

    Inputs:
     - library:     library prefix parsed out of a FASTQ filename by FASTQ_NAME_RE
     - sample_name: sample name from the metadata
     - mode:        library mode token, 'GEX' or 'SB'
    Output:
     - True if the library belongs to that sample's library of that mode
    """

    pattern = rf'{re.escape(sample_name)}[_-]?(?i:{re.escape(mode)})(?:[_-].*)?'
    if re.fullmatch(pattern, library):
        return True

    # a library named after the sample alone is its gene-expression library
    return mode == 'GEX' and library == sample_name


def find_sample_fastqs(
    fastq_path: Path | str,
    sample_name: str,
    mode: str,
    lanes: list[str] | None = None,
    reads: tuple[str, ...] | None = None
) -> list[Path]:
    """
    Find the FASTQs in an externally supplied directory (--fastqs) that belong to one library of
    one sample.

    Inputs:
     - fastq_path:  directory of already-demultiplexed FASTQs (searched recursively)
     - sample_name: sample name from the metadata
     - mode:        library mode token, 'GEX' or 'SB'
     - lanes:       lane tokens to restrict to (from declared_lanes); None/empty means all lanes
     - reads:       read tokens to restrict to, e.g. ('R1',); None means any read
    Output:
     - sorted list of matching FASTQ paths
    """

    matches = []
    for file in sorted(Path(fastq_path).rglob('*.fastq.gz')):
        parsed = FASTQ_NAME_RE.match(file.name)
        if parsed is None:
            continue
        if not library_matches_sample(parsed.group('lib'), sample_name, mode):
            continue
        if reads is not None and parsed.group('read') not in reads:
            continue
        if lanes and f"L{parsed.group('lane')}" not in lanes:
            continue
        matches.append(file)

    return matches


def missing_fastqs(
    fastq_path: Path | str,
    metadata_df: pd.DataFrame,
    modes: tuple[str, ...] = FASTQ_MODES
) -> dict:
    """
    Check that a directory supplied with --fastqs holds every FASTQ the run needs: an R1 and an R2
    for each sample's requested libraries, in each lane the metadata declares.

    This is the completeness guardrail --fastqs exists to make possible. Because the user stated
    explicitly that these are demultiplexed FASTQs, a failure here is unambiguously "your input is
    incomplete" and is reported as an error, rather than being one of two indistinguishable
    readings of a failed schema check.

    Inputs:
     - fastq_path:  directory of already-demultiplexed FASTQs
     - metadata_df: sample metadata (needs 'Sample Name', 'Lane', 'SB Lane', 'Chemistry')
     - modes:       library modes to require, a subset of FASTQ_MODES; callers narrow this so a
                    stage is never blocked by a library it does not consume (a count-only run must
                    not fail because the spatial-barcode FASTQs are absent)
    Output:
     - {sample name: {mode: [missing read tokens]}}; empty when nothing is missing
    """

    fastq_path = Path(fastq_path)
    missing = {}
    lane_columns = {'GEX': 'Lane', 'SB': 'SB Lane'}

    for _, sample in metadata_df.iterrows():
        sample_name = str(sample['Sample Name'])

        required = [mode for mode in FASTQ_MODES if mode in modes]
        # Flex spatial barcodes are supplied separately via the workflow.flex_spatial_R1_path /
        # flex_spatial_R2_path config fields and demultiplexed by the Takara/Trekker path, so no
        # spatial-barcode library is expected alongside the GEX FASTQs for a Flex sample
        if _is_flex(sample):
            required = [mode for mode in required if mode != 'SB']

        for mode, lane_column in ((mode, lane_columns[mode]) for mode in required):
            lanes = declared_lanes(sample, lane_column)
            found = find_sample_fastqs(fastq_path, sample_name, mode, lanes)
            present = {FASTQ_NAME_RE.match(file.name).group('read') for file in found}
            for read in ('R1', 'R2'):
                if read not in present:
                    missing.setdefault(sample_name, {}).setdefault(mode, []).append(read)

    return missing


def stage_spatial_fastqs(
    fastq_path: Path | str,
    sample_name: str,
    lanes: list[str] | None = None
) -> Path:
    """
    Collect one sample's spatial-barcode FASTQs out of an externally supplied directory into a
    dedicated staging directory of symlinks, and return that directory.

    spatial_count.jl takes a *directory* and selects files by substring-matching the sample ID, so
    handing it a mixed --fastqs directory would let the sample's gene-expression reads match the
    same ID and be counted as spatial reads. Staging just the spatial-barcode library keeps that
    selection exact, and symlinking avoids duplicating FASTQs that are routinely hundreds of GB.

    Inputs:
     - fastq_path:  directory of already-demultiplexed FASTQs
     - sample_name: sample name from the metadata
     - lanes:       lane tokens to restrict to (from declared_lanes); None/empty means all lanes
    Output:
     - path to the staging directory holding only this sample's spatial-barcode FASTQs
    """

    fastq_path = Path(fastq_path)
    matches = find_sample_fastqs(fastq_path, sample_name, 'SB', lanes)
    if not matches:
        log_write(f"[ERROR]: no spatial-barcode FASTQs found for sample {sample_name} in {fastq_path}")
        log_write("Troubleshooting:")
        log_write(f" • FASTQs passed with --fastqs must contain the sample name and an `SB` token, "
                  f"e.g. {sample_name}_SB_S1_L001_R1_001.fastq.gz")
        log_write(" • Check that the lanes in the `SB Lane` metadata column match the lanes in the filenames")
        sys.exit(1)

    # rebuild the staging directory from scratch so a re-run can never mix in links left over from
    # a previous run's (possibly different) --fastqs directory
    stage_dir = TMP_PATH / 'spatial_fastqs' / sample_name
    if stage_dir.is_dir():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True)

    for file in matches:
        link = stage_dir / file.name
        if link.exists():
            log_write(f"[ERROR]: two spatial-barcode FASTQs for sample {sample_name} share the filename "
                      f"{file.name} in different subdirectories of {fastq_path}")
            log_write("Troubleshooting:")
            log_write(" • spatial_count.jl processes a flat directory, so filenames must be unique")
            log_write(" • Flatten the directory passed to --fastqs, or rename the duplicate files")
            sys.exit(1)
        link.symlink_to(file.resolve())

    return stage_dir


def sanitize_path_component(value: str, field_name: str = "value") -> str:
    """
    Restrict a metadata-sourced string (e.g. Puck ID, Sample Name, BCL ID) to
    safe filename characters before it is used to build a filesystem path.

    Rejects path separators, '..', and any other character outside a
    conservative allowlist, so a value like '../../etc' or an absolute path
    can't escape the intended output directory.

    Inputs:
     - value:              the raw metadata value to check
     - field_name:         human-readable name used in the error message
    Outputs:
     - the value, unchanged, if it passes validation
    """

    value = str(value).strip()

    if not re.fullmatch(r'[A-Za-z0-9._-]+', value) or value in {'.', '..'}:
        raise ValueError(
            f"Invalid {field_name}: {value!r}. Only letters, digits, '.', '_', and '-' are allowed "
            "(no path separators or '..').\n"
            "Troubleshooting:\n"
            f" • This value is used to build a filename, so correct the `{field_name}` metadata column to a plain identifier\n"
            " • A stray space, '/' or newline pasted into the cell is the usual cause\n"
            " • Give the bare identifier here, not a path -- the directory it lives in comes from the config file"
        )

    return value


def retrieve_takara_bead_barcode_file(
    tile_id: str,
    outpath: str | Path
) -> Path:
    """
    Programmatically download a bead barcodes file from the Takara website by ID

    Inputs:
     - tile_id:           ID of the puck (tile) you're trying to download
     - outpath:           Path to directory in which the barcodes file should be saved
    Outputs:
     - save_path:         Path of the saved barcodes zip file
    """

    tile_id = sanitize_path_component(tile_id, "Puck ID")

    # The download drives a real headless browser, which needs two separate things installed: the
    # `playwright` Python package, and the Chromium *binary* it controls (a separate ~150MB download
    # into ~/.cache/ms-playwright that `pip install playwright` does NOT fetch). Both are reported
    # here with the command that fixes them, because the native errors are unhelpful: a missing
    # package raises a bare ImportError, and a missing browser raises a long Playwright message that
    # buries the one line that matters. Callers only reach this function when the puck CSV is absent,
    # so supplying the file by hand is always a valid way out.
    manual_hint = (f" • Or skip the download entirely by placing {tile_id}_BeadBarcodes.csv in "
                   f"{outpath} yourself -- an existing file is used as-is")

    try:
        from playwright.sync_api import Error as PlaywrightError, sync_playwright
    except ImportError as exc:
        log_write(f"[ERROR]: the `playwright` package is required to download puck {tile_id}'s bead barcodes, but is not installed: {exc}")
        log_write("Troubleshooting:")
        log_write(" • Install it and its browser: `uv add playwright && uv run playwright install chromium`")
        log_write(manual_hint)
        sys.exit(1)

    with sync_playwright() as p:
        # open a chromium browser and navigate to the Takara page
        try:
            browser = p.chromium.launch()
        except PlaywrightError as exc:
            log_write(f"[ERROR]: could not launch a headless Chromium browser to download puck {tile_id}'s bead barcodes")
            log_write("Troubleshooting:")
            log_write(" • The Chromium binary is a separate download from the Python package; install it with `uv run playwright install chromium`")
            log_write(" • It is cached in ~/.cache/ms-playwright, which must exist and be readable by the user running the pipeline")
            log_write(" • On a headless VM you may also need the system libraries: `uv run playwright install-deps chromium`")
            log_write(manual_hint)
            log_write("  Playwright reported:")
            for line in str(exc).strip().splitlines()[:5]:
                log_write(f"   {line}")
            sys.exit(1)
        page = browser.new_page()
        page.goto("https://c2ack773.caspio.com/dp/7391b0003c9dcd640f5348308dc9")

        # identify the puck ID submission form and enter the ID
        page.locator("input[name='InsertRecordSearchedTileID']").fill(tile_id)

        # submit form
        page.locator("input[type='submit']").click()

        # click the download button and download the zipped barcode file
        with page.expect_download() as download_info:
            page.locator("a[href$='BeadBarcodes.zip']").click()

        download = download_info.value

        outpath = Path(outpath)
        if not outpath.is_dir():
            outpath.mkdir()

        # save the barcodes file and close the browser
        save_path = outpath / f"{tile_id}_BeadBarcodes.zip"
        if save_path.resolve().parent != outpath.resolve():
            raise ValueError(
                f"Refusing to save outside of {outpath}: resolved to {save_path.resolve()}\n"
                "Troubleshooting:\n"
                f" • The puck ID {tile_id!r} resolved to a path outside the intended output directory\n"
                " • Correct the `Puck ID` metadata column to a plain identifier with no path separators\n"
                f" • Check {outpath} is not a symlink pointing somewhere unexpected"
            )
        download.save_as(save_path)
        browser.close()

        return save_path


def downsample_spatial(
    in_path: Path | str,
    out_path: Path | str,
    rate: float = 0.2
):
    seed = 42
    chunk_size = 50_000_000   # rows per chunk; adjust down if still memory-tight

    with h5py.File(in_path, "r") as fin:
        n_rows = fin["matrix/cb_index"].shape[0]
        rng = np.random.default_rng(seed)

        # Generate a boolean keep-mask via chunked random draws (avoids allocating
        # one full-size bool array if n_rows is enormous; 1.29B bools = ~1.3GB, so
        # actually fine to do in one shot on a highmem box, but chunked is safer)
        n_keep = int(n_rows * rate)
        keep_idx = rng.choice(n_rows, size=n_keep, replace=False)
        keep_idx.sort()   # critical: sorted indices make HDF5 fancy-indexing far faster

        with h5py.File(out_path, "w") as fout:
            # Copy small groups untouched
            fin.copy("lists", fout)
            fin.copy("puck", fout)
            fin.copy("metadata", fout)

            # Chunked copy of the big parallel arrays using the shared keep_idx
            reads_kept_total = 0
            for key in ["matrix/cb_index", "matrix/sb_index", "matrix/umi", "matrix/reads"]:
                dtype = fin[key].dtype
                out_ds = fout.create_dataset(key, shape=(n_keep,), dtype=dtype,
                                            chunks=True, compression="gzip")
                write_pos = 0
                for start in range(0, n_rows, chunk_size):
                    end = min(start + chunk_size, n_rows)
                    # find which keep_idx fall in this chunk
                    lo = np.searchsorted(keep_idx, start)
                    hi = np.searchsorted(keep_idx, end)
                    local_idx = keep_idx[lo:hi] - start
                    if len(local_idx) == 0:
                        continue
                    chunk_data = fin[key][start:end][local_idx]
                    out_ds[write_pos:write_pos+len(chunk_data)] = chunk_data
                    write_pos += len(chunk_data)
                    log_write(f"{key}: wrote rows {write_pos}/{n_keep}", end="\r")
                if key == "matrix/reads":
                    reads_kept_total = int(out_ds[:].sum())
            print()

            # load_matrix.R enforces an exact bookkeeping identity:
            #   reads_total == reads_noumi + reads_noup + reads_nosb + sum(matrix/reads)
            # reads_noumi/noup/nosb are derived from metadata/UP_matching and
            # metadata/SB_matching count arrays, which -- like reads_total -- describe
            # reads that never made it into the matrix at all, so they can't be
            # recovered from keep_idx. We rescale them by the same `frac` used for the
            # matrix and then force reads_total to exactly balance, rather than
            # independently rescaling reads_total (which would only balance the
            # identity approximately and still fail the stopifnot in R).
            def decode(arr):
                return [x.decode() if isinstance(x, bytes) else x for x in arr]

            up_type = decode(fin["metadata/UP_matching/type"][:])
            up_count = fin["metadata/UP_matching/count"][:]
            sb_type = decode(fin["metadata/SB_matching/type"][:])
            sb_count = fin["metadata/SB_matching/count"][:]

            scaled_up_count = np.round(up_count * rate).astype(up_count.dtype)
            scaled_sb_count = np.round(sb_count * rate).astype(sb_count.dtype)
            fout["metadata/UP_matching/count"][:] = scaled_up_count
            fout["metadata/SB_matching/count"][:] = scaled_sb_count

            noumi_idx = [i for i, t in enumerate(up_type) if t in ("umi_N", "umi_homopolymer")]
            noup_idx = [i for i, t in enumerate(up_type) if t in ("none", "GG")]
            # must match load_matrix.R's reads_nosb subset exactly: none + HD1ambig + SB2delambig
            # (the second-half-deletion ambiguous class). Omitting SB2delambig here rescales the
            # downsampled num_reads off by that count and breaks load_matrix.R's stopifnot.
            nosb_idx = [i for i, t in enumerate(sb_type) if t in ("none", "HD1ambig", "SB2delambig")]
            reads_noumi = int(scaled_up_count[noumi_idx].sum())
            reads_noup = int(scaled_up_count[noup_idx].sum())
            reads_nosb = int(scaled_sb_count[nosb_idx].sum())

            new_reads_total = reads_noumi + reads_noup + reads_nosb + reads_kept_total
            fout["metadata/num_reads"][()] = new_reads_total

    log_write(f"Done. Kept {n_keep:,} / {n_rows:,} rows ({rate:.0%})")
    log_write(f"Rescaled metadata/num_reads -> {new_reads_total:,} to preserve load_matrix.R's read-accounting invariant")


# Width of the timestamp column every progress line is written in ("HH:MM:SS" plus two spaces).
# Continuation lines are indented to match, so the log reads as one aligned column of times with
# wrapped detail hanging beneath each.
TS_WIDTH = 10


def format_duration(seconds: float) -> str:
    """
    Render an elapsed time compactly: '45s', '12m34s', '2h07m'.

    Per-stage durations are the main thing anyone reads a pipeline log for -- they say which stage to
    optimize and how to size the next allocation -- so they are written out rather than left implicit
    in a pair of timestamps the reader has to subtract.

    Inputs:
     - seconds: elapsed seconds
    Output:
     - a short human-readable duration
    """

    seconds = int(round(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}h{minutes:02d}m"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def run_relative(path: Path | str) -> str:
    """
    Render a path relative to this run's directory when it sits inside it, and unchanged otherwise.

    The run root is stated once in the log header, so repeating it on every subsequent line costs
    width and adds nothing -- a full absolute path is preserved only when it points somewhere else
    (an input directory, a bucket, a software install), which is exactly when it carries information.

    Inputs:
     - path: any path or path-like value
    Output:
     - the path as a string, shortened to a run-relative one where possible
    """

    try:
        return str(Path(path).relative_to(RUN_PATH))
    except (ValueError, TypeError):
        return str(path)


def log_ts(
    message: str,
    details: list[str] | None = None,
    **kwargs
) -> None:
    """
    Write a progress line stamped with the time of day, plus optional detail lines aligned beneath it.

    Only the time is shown, not the date: the date is in the run header, and repeating it on every
    line pushes the message itself off the edge of a terminal.

    Inputs:
     - message: the line to write
     - details: extra lines, indented to the timestamp column
    """

    stamp = datetime.now().strftime('%H:%M:%S')
    lines = [f"{stamp}  {message}"]
    for detail in details or []:
        lines.append(f"{'':<{TS_WIDTH}}{detail}")
    log_write(lines, **kwargs)


def log_detail(details: list[str] | str, **kwargs) -> None:
    """
    Write continuation lines aligned under the timestamp column, with no timestamp of their own.

    Inputs:
     - details: line or lines to write
    """

    if isinstance(details, str):
        details = [details]
    log_write([f"{'':<{TS_WIDTH}}{detail}" for detail in details], **kwargs)


def log_section(title: str, width: int = 74) -> None:
    """
    Write a lightweight rule introducing a stage.

    Deliberately quieter than the run header: a stage boundary should be findable when skimming
    without competing with the run's own opening and closing banners. Only stages that actually run
    get one -- a banner for something that was skipped is noise.

    Inputs:
     - title: stage name
     - width: total rule width
    """

    prefix = f"── {title} "
    log_write(["", f"{prefix}{'─' * max(0, width - len(prefix))}"])


def log_run_started(message: str) -> None:
    """
    Log the start of the run and send the corresponding Slack alert.

    The alert is fired from here rather than inferred from the message text inside log_write, so it
    happens exactly once per run and cannot be lost by rewording a log line.

    Inputs:
     - message: the line to log
    """

    log_ts(message)
    if ALERTS:
        slack_message(f"🚀 {message}")


def notify_complete(message: str) -> None:
    """
    Send the Slack alert for a completed run.

    Explicit for the same reason log_run_started() is: the alert used to depend on log_write spotting
    a '[DONE]' substring, which forced a machine-readable marker into the banner a human reads.

    Inputs:
     - message: the alert text
    """

    if ALERTS:
        slack_message(f"✅ {message}")


def log_skipped(label: str, reason: str) -> None:
    """
    Note a stage that was requested but turned out not to need running.

    Written as a single timestamped line with no section rule: a stage that did nothing should not
    look, at a glance, like one that did. Stages that were never requested are not reported here at
    all -- the run header already lists them.

    Inputs:
     - label:  stage name
     - reason: why it was skipped
    """

    log_ts(f"skipped {label} — {reason}")


def run_banner(text: str, fill: str = '═', width: int = 74) -> str:
    """
    Centre `text` in a full-width rule, for the run's opening and closing lines.

    Inputs:
     - text:  the text to centre
     - fill:  rule character
     - width: total width
    Output:
     - the banner line
    """

    pad = max(0, width - len(text) - 2)
    left = pad // 2
    return f"{fill * left} {text} {fill * (pad - left)}"


def log_stage_start(label: str, details: list[str] | None = None) -> float:
    """
    Announce a stage and return its start time for log_stage_done().

    Inputs:
     - label:   stage name, as it should read in the log
     - details: extra lines to write beneath the announcement
    Output:
     - the stage's start time, to be handed back to log_stage_done()
    """

    log_section(label)
    started = time.time()
    log_ts("started", details)
    return started


def log_stage_done(
    started: float,
    outlog: Path | str | None = None,
    outdir: Path | str | None = None
) -> None:
    """
    Close out a stage with a tick and its elapsed time, plus where its output went.

    Inputs:
     - started: the value returned by log_stage_start()
     - outlog:  file holding the stage's own console output, if any
     - outdir:  directory the stage's outputs were written to, if any
    """

    details = []
    if outlog is not None:
        details.append(f"log     → {run_relative(outlog)}")
    if outdir is not None:
        details.append(f"outputs → {run_relative(outdir)}")
    log_ts(f"✓ finished in {format_duration(time.time() - started)}", details)


def log_write(
    messages: list[str] | str = "",
    logpath: Path = RUNTIME_LOG,
    mode: str = "a",
    terminal: bool = True,
    terminator: str = "\n"
):
    """
    Write provided messages to a logfile; optionally prints to stdout.

    Inputs:
     - messages: string or list of strings to write
     - logpath:  destination logfile (default: RUNTIME_LOG)
     - mode:     file open mode ("w" or "a")
     - terminal:   also print to stdout/stderr (default: True)
    """

    # instantiate the Console objects for stdout and stderr
    console = Console()
    error_console = Console(stderr=True)

    if isinstance(messages, str):
        messages = [messages]

    # open the corresponding logfile
    with open(logpath, mode) as logfile:
        for message in messages:
            # write message to logfile
            logfile.write(message + terminator)
            if terminal:
                # if terminal printing is enabled, print message to terminal
                if 'ERROR' in message:
                    if '[ERROR]' in message:
                        if ALERTS:
                            # alert the user if Slack alerts are enabled
                            slack_message(f"🚨 {message}\nPlease see the logfile for details on the crash.")
                        # format the message with color printing codes
                        message = message.replace('[ERROR]', '[bold red]\\[ERROR][/bold red]')
                    elif '[INTERNAL ERROR]' in message:
                        message = message.replace('[INTERNAL ERROR]', '[bold red]\\[INTERNAL ERROR][/bold red]')
                        if ALERTS:
                            slack_message(f"🚨 {message}\nPlease see the logfile for details on the crash.")

                    # print to the stderr console and flush
                    error_console.print(message, end=terminator)
                    error_console.file.flush()

                    # log total runtime and result
                    end_time = time.time()
                    elapsed_time = end_time - START_TIME
                    hours, remainder = divmod(elapsed_time, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    log_write(f"  Total runtime:        {int(hours):02d}:{int(minutes):02d}:{seconds:02.0f}", SUMMARY_LOG, terminal=False)
                    log_write(f"  Result:               FAILURE", SUMMARY_LOG, terminal=False)

                    continue

                # format messages with color printing codes
                elif '[WARNING]' in message:
                    message = message.replace('[WARNING]', '[bold yellow]\\[WARNING][/bold yellow]')
                elif '[DONE]' in message:
                    if ALERTS:
                        # alert user if Slack alerts are enabled. The newline stripping is done
                        # before the f-string rather than inside it: a backslash in an f-string
                        # expression is a syntax error before Python 3.12, and pyproject.toml
                        # still supports 3.11.
                        alert_text = message.replace('\n', '')
                        slack_message(f"✅ {alert_text}")
                    message = message.replace('[DONE]', '[bold green]\\[DONE][/bold green]')
                # NB: there is deliberately no "run started" trigger here. It used to fire on the
                # substring 'started\n', which alerted once per stage that announced itself, and would
                # have stopped firing silently the moment a message was reworded or lost its trailing
                # newline. log_run_started() now sends that alert explicitly, once per run.

                # print to console and flush
                console.print(message, end=terminator)
                console.file.flush()


# Stage-specific troubleshooting hints for job_crash(). Keyed on the `job_name` each caller passes
# in, so a crash report names the things that actually go wrong in *that* tool rather than offering
# the same generic advice for every stage. A job with no entry here still gets the common hints
# below, so adding a new stage never has to touch this table.
_CRASH_HINTS = {
    'cellranger mkfastq': [
        "Cellranger refuses to write into an existing output directory -- re-run with --force to clear it",
        "An index collision or an index that matches no reads is the usual cause; check the `RNA Index`/`SB Index` and `Lane`/`SB Lane` metadata columns against the run",
        "`bcl2fastq` must be on PATH for mkfastq to work; check it was found in the software scan above",
        "Too little memory or disk also shows up here -- raise `settings.memory` or free space under `paths.output_path`",
    ],
    'cellranger count': [
        "Cellranger refuses to write into an existing output directory -- re-run with --force to clear it",
        "Check the `Chemistry` metadata column matches how the library was prepared; a wrong chemistry fails during read parsing",
        "Check the reference genome matches the sample's species (`settings.reference_genome` / the `Species` column)",
        "Out-of-memory kills are common here; raise `settings.memory` (currently the value passed as --localmem)",
    ],
    'cellbender': [
        "A CUDA out-of-memory error means the GPU is too small for this matrix -- lower `workflow.cellbender_total_droplets`",
        "Check `workflow.cellbender_total_droplets` exceeds `workflow.cellbender_estimated_cells`; CellBender rejects the reverse",
        "A failure to converge usually means the droplet/cell estimates are wrong for this sample -- see the barcode-rank plot in the count stage's web_summary.html",
    ],
    'generate_puck_csv.jl': [
        "Check the puck's raw barcode directory holds both BeadBarcodes.txt and BeadLocations.txt",
        "Check `paths.raw_barcodes_path` points at the directory *containing* the per-puck subdirectories",
        "Or supply a pre-built <puck id>.csv in `paths.puck_path` to skip generation entirely",
    ],
    'spatial_count.jl': [
        "Check the spatial-barcode FASTQs are complete and not truncated -- a mismatched R1/R2 read count fails here",
        "Check the `Puck ID` metadata column names a puck whose CSV has 3 columns (barcode, x, y) and no header",
        "Out-of-memory kills are common on deeply sequenced samples; set `workflow.spatial_downsampling` (e.g. 0.5) to thin the reads first",
    ],
    'run_spatial.R': [
        "Check the count (or cellbender) outputs for every sample are present and complete",
        "Check the spatial-barcode counting stage produced an SBcounts.h5 for every sample",
        "A sample with too few positioned cells fails here; check the per-sample bead counts in the log",
    ],
    'trekker_demux': [
        "Check `workflow.flex_spatial_R1_path` and `workflow.flex_spatial_R2_path` point at the Flex spatial FASTQs",
        "A duplicate probe barcode across two samples fails here; check the `Flex Probe Barcode IDs` column for reuse",
    ],
    'trekker_flex': [
        "Check the demultiplexed FASTQs for this partition exist under the run's flex/demux/ directory",
        "Check the puck barcode file for this sample downloaded successfully into flex/pucks/",
    ],
    'trekker_merge': [
        "Check every partition of this sample completed the profiling step above",
        "A partition that produced no positioned nuclei has no output to merge; check the per-partition logs",
    ],
}


def job_crash(
    job_name: str,
    return_code: int,
    logfile: str
) -> None:
    """
    Log a subprocess failure with its exit code and logfile location, then exit

    Inputs:
     - job_name:    human-readable name of the failed job (used in the error message)
     - return_code: non-zero exit code returned by the subprocess
     - logfile:     path to the file containing the subprocess stdout/stderr output
    """

    # log crash
    log_ts(f"[ERROR]: {job_name} returned exit code {return_code}",
           [f"see {run_relative(logfile)} for the tool's own output"])

    log_write("Troubleshooting:")
    log_write(f" • Read the tail of that log first -- it names the real cause: `tail -50 {logfile}`")

    # exit code 137 is SIGKILL, which on a shared machine or a cgroup-limited job is nearly always
    # the OOM killer rather than anything wrong with the inputs; say so instead of sending the user
    # to look for a configuration mistake that isn't there
    if return_code in (137, -9):
        log_write(" • Exit code 137 means the job was killed with SIGKILL, which almost always means it ran out of memory")
        log_write(" • Raise `settings.memory`, or give the run a bigger machine (--machine for --gcp, a larger allocation for --slurm)")
    elif return_code in (139, -11):
        log_write(" • Exit code 139 is a segmentation fault inside the tool itself; check you are on a supported version")

    for hint in _CRASH_HINTS.get(job_name, []):
        log_write(f" • {hint}")

    log_write(" • Once the cause is fixed, re-run this stage alone rather than the whole pipeline; add --force to overwrite the partial outputs")
    log_write("Shutting down... x_x")

    # exit
    sys.exit(1)


def job_success(
    job_name: str,
    outlog: Path | str,
    outdir: Path | str,
    started: float | None = None
) -> None:
    """
    Log the successful completion of a stage, with its elapsed time and where its output went.

    The Slack alert is sent explicitly rather than by log_write noticing a '[DONE]' substring, which
    lets the log line itself stay clean -- the reader gets a tick and a duration, not a marker that
    only ever existed to drive alerting.

    Inputs:
     - job_name: human-readable name of the completed job
     - outlog:   path to the file containing the subprocess stdout/stderr output
     - outdir:   path to the directory where output files were written
     - started:  the stage's start time from log_stage_start(), if known, for the duration
    """

    suffix = f" in {format_duration(time.time() - started)}" if started is not None else ""
    log_ts(f"✓ {job_name} finished{suffix}", [
        f"log     → {run_relative(outlog)}",
        f"outputs → {run_relative(outdir)}",
    ])

    if ALERTS:
        slack_message(f"✅ {job_name} finished successfully{suffix}")
