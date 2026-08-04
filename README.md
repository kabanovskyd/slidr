<p align="center">
  <img src="assets/logo_small.png" alt="slidr logo" width="420"/>
</p>

# `slidr` — The Slide-Tag Analysis Pipeline

A pipeline for processing Slide-Tag spatial transcriptomics data. Starting from BCL files or pre-demultiplexed FASTQs, `slidr` runs demultiplexing, gene expression quantification, ambient RNA removal, and spatial barcode assignment to produce analysis-ready outputs.



## Pipeline stages

| Stage | Tool | Flag |
|---|---|---|
| BCL → FASTQ demultiplexing | Cellranger mkfastq | `--mkfastq` |
| Gene expression quantification | Cellranger count | `--count` |
| Ambient RNA removal | CellBender | `--cellbender` |
| Spatial barcode counting | Julia (spatial_count.jl) | `--spatial-count` |
| Spatial analysis | R / Seurat | `--spatial-analysis` |



## Requirements

**Must be installed manually** (proprietary licences):
- [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/downloads) — 8.0.1 by default;
  a per-sample `Cellranger` metadata column can select another installed release, and Flex chemistry
  needs 9.x for `cellranger multi`
- [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) ≥ 2.20

**Installed automatically by `./slidr`:**
- [uv](https://docs.astral.sh/uv/) — Python package management
- [Miniforge / conda](https://github.com/conda-forge/miniforge) — R environment management
- [Julia](https://julialang.org/) — spatial barcode processing

**Python dependencies** are managed by `uv` via `pyproject.toml`.  
**R dependencies** are managed by `conda` via `envs/conda.yml`.

---

## Installation

```bash
git clone https://github.com/kabanovskyd/slidr
cd slidr
```

No further setup is required — `./slidr` installs missing dependencies on first run.

---

## Configuration

Edit `config/config.yaml` before running:

```yaml
paths:
  output_path: /data/slidetag           # where outputs are written
  input_path: /mnt/sequencer            # root directory containing BCL run folders
  software_path: /mnt/lab_software      # directory scanned for Cellranger / bcl2fastq.
                                        # May also be a gs:// prefix, staged with --stage-gcs
  raw_barcodes_path: /mnt/barcodes      # raw puck barcode files (BeadBarcodes.txt / BeadLocations.txt)
  puck_path: /mnt/pucks                 # processed puck CSV files (generated if absent)
  reference_path: /mnt/reference        # directory containing reference genome subdirectories
  auth_key_path: /path/to/auth_key.json # Google service account key; only needed for Sheet metadata.
                                        # May also be a gs:// object

settings:
  memory: 50                            # GB allocated to Cellranger
  threads: 16                           # cores allocated to Cellranger and Julia
  metadata_source: https://docs.google.com/spreadsheets/d/...#gid=0   # Sheet URL, or a .tsv/.csv path
  input_bucket:                         # GCS prefix remote runs stage their inputs from
  output_bucket:                        # GCS bucket outputs are uploaded to after a run
  reference_genome: refdata-gex-GRCm39-2024-A   # directory name under reference_path
  alerts: false                         # send Slack alerts on errors/completion
  slack_token:                          # file path, gs:// object, or literal token (if alerts: true)

workflow:
  generate_bam: false
  cellbender_total_droplets:            # null to auto-detect from the barcode-rank plot
  cellbender_estimated_cells:           # null to let CellBender estimate
  cellbender_epochs: 160
  cellbender_learn_rate: 0.5
  spatial_downsampling:                 # optional float; thins spatial reads before counting
  top_n_percent_umi_filter:             # optional float 0-100; bead UMI filtering percentile
  flex_emptydrops_minimum_umis: 100     # Flex chemistry only
  flex_probe_set:                       # Flex chemistry only
  flex_spatial_R1_path:                 # Flex chemistry only
  flex_spatial_R2_path:                 # Flex chemistry only
  flex_gex_fastqs:                      # Flex chemistry only; list of GEX FASTQ prefixes
```

Numeric fields must be given as plain numbers. Note that YAML reads the bare words `yes`/`no`/`on`/`off`
as booleans, so `threads: no` is rejected rather than silently treated as `0`.

When staging inputs from Google Cloud Storage (`--gcp`, or `--slurm --stage-gcs`), point `reference_path`,
`puck_path` and `raw_barcodes_path` at `gs://` locations instead of local directories.

`software_path` may also name a `gs://` location, in which case the Cellranger/bcl2fastq tree is
downloaded to `output/software/` before it is scanned. Unlike the three above it is *not* required to
be a bucket on a staged run — a local directory stays valid, because the machine running the pipeline
usually has the software installed already (the `--gcp` VM image does). Two consequences of a staged
software tree worth knowing:

- Only a full `gs://…` URL triggers it. The bare `bucket/prefix` form the other paths accept is not
  honoured here, since it cannot be told apart from a relative local directory and guessing wrong means
  a multi-GB download instead of an error. `--stage-gcs` is still required.
- The download happens once per run, lazily — a `software_cache.txt` that already pins the executables
  skips it entirely — and is reused by later runs sharing the output tree. Because GCS does not carry
  POSIX permissions, execute bits are restored on the staged copy afterwards.

### `settings.input_bucket`

The GCS prefix a run stages its inputs from. This is the only place the location is configured —
there is no command-line equivalent. Include the `gs://` scheme.

```yaml
settings:
  input_bucket: gs://slidr_data/inputs
```

What the prefix must hold depends on where the run executes, because only a VM lacks a checkout to read
its configuration from:

| Object | `--gcp` | `--slurm --stage-gcs` |
|---|---|---|
| `<BCL_ID>/` run folder | required | required |
| `config.yaml` | required | **not used** — the job runs with `config/config.yaml` from the checkout it was submitted from |
| `auth_key.json` | required | optional; staged if present, otherwise `paths.auth_key_path` from the local config applies |

Upload whichever apply before launching — the pipeline never puts them there for you. `./slidr` reads
this field locally and hands the location to the VM or compute node, which then downloads what it needs.

A Slurm job runs from a real checkout on a cluster you administer, so `config/config.yaml` is already
there and is the file you just edited; staging a copy written for another machine on top of it only
invited the two to disagree. In practice this means a config change takes effect on your next `sbatch`
with nothing to re-upload, and the CPU/memory `./slidr` requests from Slurm are read from the same file
the job itself will use.

`--gcp` always stages, since a fresh VM has nothing else to read. A `--slurm` run stages only when
asked with `--stage-gcs`; without it the job reads everything off the cluster's own filesystem. Local
runs ignore the field entirely.

### Google Sheets authentication

If you are reading input metadata from a Google Sheet, you will first need to download an authentication key that will allow your program to interface with it. You can do so by following the steps below:

1. Navigate to your project within Google Cloud Console (top dropdown)
2. Navigate to the left sidebar -> APIs & Services -> Enabled APIs & Services
3. Check that the list of enabled services includes `Google Sheets API` and `Google Drive API`; if not, add them with the "Enable APIs and Services" button
4. Navigate to the "Credentials" tab in the "APIs & Services" sidebar
5. Click on "Create credentials" button and follow the steps to create a service account
6. Under the "Service Accounts" section on the same page, click on the email of the account you just created, navigate to the "Keys" tab, and click on "Add keys" -> "Create new key" -> "JSON"
7. Connect your service account to the Google Sheet: navigate to the sheet URL, click the "Share" button, paste the service account email, and give it "Editor" permissions
8. Add the path to the downloaded JSON key to the `paths.auth_key_path` field in the configfile and **keep it secure** - __never share, commit, or otherwise expose your key to people outside your organization__. The field also accepts a `gs://` object, so the key can live in a bucket rather than on every machine that runs the pipeline.

### Software cache initialization

`Slidr` will automatically search for matching software executables in the provided software directory and populate a cache file (`software_cache.txt`) in the project root to save executable paths for future runs. However, this step make take a while - to skip this (or specify a specific executable version if you have several versions installed), create the software cache file and add absolute paths on individual lines:

```
/path/to/lab_software/cellranger-8.0.1/bin/cellranger
/path/to/lab_software/bcl2fastq2_v2.20.0/bin/bcl2fastq
/home/unix/your_username/.juliaup/bin/julia
...
```

### Metadata format

Sample metadata must be supplied as a Google Sheet or a `.tsv`/`.csv` file, and must contain all of the
following columns:

| Column | Description |
|---|---|
| `Run` | `YES` to include this sample, anything else to skip |
| `Email` | User email, used for Slack notifications |
| `Sample Name` | Unique sample identifier |
| `BCL` | BCL run ID (matched to `--bcl`); the run folder is resolved as `input_path/<BCL>` |
| `Species` | `Human` or `Mouse` (selects the reference genome automatically) |
| `Chemistry` | Sequencing chemistry, e.g. `3Pv3`, `5P`, `Flex` |
| `RNA Index` | Cellranger index for the RNA library |
| `Lane` | Sequencing lane(s), comma-separated |
| `SB Index` | Cellranger index for the spatial barcode library |
| `SB Lane` | Lane(s) for the spatial barcode library, comma-separated |
| `Puck ID` | Identifier matching a puck CSV in `puck_path` |

These columns are optional, and only read when non-empty:

| Column | Description |
|---|---|
| `Merge RNA From BCL` | Additional BCL to merge RNA FASTQs from |
| `Merge Spatial From BCL` | Additional BCL to merge spatial FASTQs from |
| `Add RNA Index` | Index for the merged-in RNA library |
| `Add SB Index` | Index for the merged-in spatial library |
| `Add Puck ID` | Puck override for a merged sample |
| `Cellranger` | Per-sample Cellranger version override, e.g. `V8` or `8.0.1` (default `8.0.1`; `9.0.1` for Flex) |
| `Flex Probe Barcode IDs` | Probe barcode IDs (required for `Flex` chemistry), separated by `,` or `\|` |

For a Google Sheet, the worksheet tab is taken from the `#gid=` fragment of the URL. Column names are
matched exactly, including case and spacing.

> A ready-to-edit `example_metadata.tsv` is included in the repository root, carrying every column
> above. Its rows show a two-sample 3' run, a skipped row (`Run: NO`), a cross-BCL RNA merge, a
> `Cellranger` override and a Flex pool. Optional columns you leave empty are ignored.

---

## Usage

`./slidr` is the single entry point for local, GCP and Slurm runs. It installs any missing
dependencies, then either runs the pipeline in your current shell or dispatches it to a VM or the
scheduler:

```bash
./slidr --bcl <BCL_ID> [options]
```

### Run the full pipeline

```bash
# locally
./slidr --bcl RUN_20240101 --run-all

# on a GCP VM (streams logs back; see settings.input_bucket above)
./slidr --bcl RUN_20240101 --gcp --project my-project --run-all

# as a Slurm job
./slidr --bcl RUN_20240101 --slurm --run-all
```

### Run individual stages

```bash
# Demultiplexing only
./slidr --bcl RUN_20240101 --mkfastq

# Gene expression quantification (supply FASTQs directly)
./slidr --bcl RUN_20240101 --count --fastqs /path/to/fastqs

# Ambient RNA removal
./slidr --bcl RUN_20240101 --cellbender

# Spatial barcode counting
./slidr --bcl RUN_20240101 --spatial-count

# Spatial analysis
./slidr --bcl RUN_20240101 --spatial-analysis

# Everything except ambient RNA removal, overwriting existing outputs
./slidr --bcl RUN_20240101 --run-all --no-cellbender --force
```

### All options

Run `./slidr --help` for the authoritative list.

**Stage selection and general options**

| Flag | Description |
|---|---|
| `--bcl BCL_ID` | BCL run ID (**required**, except with `--version`/`--help`) |
| `--run-all` | Run all pipeline stages end-to-end (the default if no stage flag is given) |
| `--mkfastq` | Run Cellranger mkfastq only |
| `--count` | Run Cellranger count only |
| `--cellbender` | Run Cellbender only |
| `--no-cellbender` | Skip Cellbender when running the full pipeline |
| `--spatial-count` | Run spatial barcode counting only |
| `--spatial-analysis` | Run spatial analysis only |
| `--fastqs [PATH]` | Run on already-demultiplexed FASTQs; skips mkfastq (mutually exclusive with `--mkfastq`). PATH is optional — with no value the FASTQs are read from `paths.input_path` |
| `--force` | Overwrite existing outputs and re-run |
| `--version`, `-v` | Print the pipeline version and exit |
| `--help`, `-h` | Print usage and exit |

**Where the run executes**

| Flag | Description |
|---|---|
| `--gcp` | Run on a GCP VM, streaming logs back to your terminal |
| `--slurm` | Submit as a Slurm batch job |
| *(no flag)* | Where inputs are staged from is configuration, not a flag — set `settings.input_bucket` in `config/config.yaml`. `--stage-gcs` is what opts a `--slurm` run into using it; `--gcp` always does |

**GCP options** (only with `--gcp`)

| Flag | Default | Description |
|---|---|---|
| `--project PROJECT` | none — pass explicitly | GCP project ID |
| `--zone ZONE` | `us-central1-a` | Compute zone |
| `--machine TYPE` | `n1-standard-16` | Machine type |
| `--disk SIZE` | `200GB` | Boot disk size |
| `--gpu [TYPE]` | off; `nvidia-tesla-t4` when bare | Attach a GPU; TYPE is an optional GCE accelerator name (e.g. `nvidia-l4`) |

**Slurm options** (only with `--slurm`; CPUs and memory come from `settings.threads`/`settings.memory`)

| Flag | Default | Description |
|---|---|---|
| `--partition NAME` | cluster default | Partition to submit to |
| `--time LIMIT` | `24:00:00` | Job time limit |
| `--workdir PATH` | derived — see below | Where to stage inputs and write outputs |

For a staged run the workdir is normally derived, so you rarely pass `--workdir`. It resolves in order:

1. `--workdir`, if given
2. the parent of `paths.output_path` in the config the job will run with
3. `$SCRATCH/slidr`, if `$SCRATCH` is set
4. `<repo>/slidr-work/<BCL_ID>` — a last resort, and usually the wrong filesystem for hundreds of GB

Whichever applies is printed at submit time along with where it came from, so a large download never
lands somewhere unexpected.

Two things a staged run needs are set up outside slidr rather than through flags:

- **GCS credentials.** Run `gcloud auth login` once on the submit node; it persists to
  `~/.config/gcloud`, which compute nodes see on a shared home. On a headless cluster use
  `gcloud auth activate-service-account --key-file=/path/to/key.json` instead, which persists the same
  way. If your home is *not* shared with compute nodes, point `CLOUDSDK_CONFIG` at somewhere that is.
- **`software_path` for this cluster**, if you would rather not edit `config/config.yaml`. Export the
  override before submitting — `sbatch` passes the environment through:
  ```bash
  export SLIDR_SOFTWARE_PATH=/opt/cellranger
  ./slidr --bcl 20240101_RUNID --slurm --stage-gcs
  ```
  The same works for `SLIDR_INPUT_PATH` and `SLIDR_OUTPUT_PATH`, though the staging script sets those two
  itself, and for `SLIDR_AUTH_KEY_PATH`, which it sets only when the bucket actually had a key.
  Alternatively, point `software_path` at a `gs://` location and let the job download the software.

`workflow/main.py` additionally accepts `--config PATH` to run against a config file outside `config/`,
and `--stage-gcs` to stage `reference_path`/`puck_path`/`raw_barcodes_path` (and a `gs://`
`software_path`) from GCS without using `--gcp`. `--stage-gcs` is what the Slurm payload passes; neither
is exposed on `./slidr` itself.

---

## Output structure

```
<output_path>/                       # `paths.output_path` from the config file
└── <BCL_ID>/                        # one directory per run, named for --bcl
    ├── log/                         # a record of how the run went
    │   ├── runtime.log              # pipeline execution log
    │   ├── <timestamp>_summary.log  # short overview of execution settings / details
    │   ├── mkfastq.log              # console output from Cellranger mkfastq
    │   ├── count.log                # console output from Cellranger count (and multi, for Flex)
    │   ├── cellbender.log           # console output from CellBender
    │   ├── spatial_barcodes.log     # console output from spatial_count.jl / generate_puck_csv.jl
    │   ├── spatial_analysis.log     # console output from the R / Seurat analysis scripts
    │   └── takara_pipeline.log      # Flex only: Trekker demux / profiling / merge output
    ├── metadata/                    # the resolved inputs describing what this run processes
    │   ├── metadata_summary.csv     # the metadata rows selected for this run
    │   ├── samplesheet.csv          # Cellranger sample sheet generated from them
    │   └── multi_samplesheet.csv    # Flex only: sample sheet for `cellranger multi`
    ├── output/                      # pipeline outputs organized by the module that produced them
    │   ├── mkfastq/                 # demultiplexed FASTQs
    │   ├── count/                   # Cellranger count outputs per sample
    │   │   └── flex/                # Flex only: one folder per probe barcode, from `cellranger multi`
    │   ├── cellbender/              # CellBender filtered matrices per sample
    │   ├── spatial_barcodes/        # SBcounts.h5 per sample
    │   ├── spatial_analysis/        # Seurat objects, coordinates, summary PDFs
    │   ├── flex/                    # Flex only: Takara / Trekker spatial intermediates
    │   │   ├── demux/               # spatial reads split per probe barcode by Trekker
    │   │   ├── pucks/               # bead barcode files downloaded from Takara
    │   │   ├── samplesheets/        # generated Trekker demux / profiling / merge sheets
    │   │   └── trekker/             # per-partition nuclei positioning output
    │   ├── reference/               # staged only: local copy of the reference genome
    │   ├── pucks/                   # staged only: puck coordinate CSVs, downloaded or generated
    │   └── barcodes/                # staged only: raw bead barcodes for puck generation
    └── tmp/                         # temporary Cellranger output directory
```

No single run contains every entry: the `flex/` tree, `multi_samplesheet.csv` and `takara_pipeline.log`
appear only for Flex chemistry, and `reference/`/`pucks/`/`barcodes/` only when staging from GCS. The
run directory is `<output_path>/<BCL_ID>/`, unless `output_path` already ends in the BCL ID — in which
case it is taken as the run directory itself rather than nested again.

### Key output files per sample

| File | Stage | Description |
|---|---|---|
| `count/<sample>/outs/filtered_feature_bc_matrix.h5` | count | Cellranger filtered expression matrix |
| `cellbender/<sample>/cellbender_output_filtered.h5` | cellbender | Ambient-corrected expression matrix |
| `spatial_barcodes/<sample>/SBcounts.h5` | spatial-count | Cell barcode → spatial coordinate mapping |
| `spatial_analysis/<sample>/seurat.qs` | spatial-analysis | Seurat object with spatial embeddings |
| `spatial_analysis/<sample>/coords.csv` | spatial-analysis | Cell coordinates |
| `spatial_analysis/<sample>/summary.pdf` | spatial-analysis | QC summary plots |

---

## Troubleshooting / known issues
- **Cellranger refusing to run**
  - If you're attempting to re-run the mkfastq and count stages of the analysis, you need to **manually delete the output directories for those modules**, as cellranger will refuse to overwrite existing files and will crash.

## Julia installation

If no Julia interpreter is found on `PATH` or in `software_cache.txt`, slidr downloads the official
installer from `install.julialang.org` and pipes it to `sh`. **This download is not checksum-verified.**

To avoid the automatic install entirely — recommended on shared or network-restricted systems — install
Julia yourself and pin the interpreter by adding its absolute path to `software_cache.txt` in the
project root:

```
/home/unix/your_username/.juliaup/bin/julia
```

slidr uses a cached path as-is and skips the download.
