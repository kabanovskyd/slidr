# slidr — Slide-Tag Analysis Pipeline

`slidr` processes Slide-Tag spatial transcriptomics data. Starting from BCL files or pre-demultiplexed FASTQs, it runs demultiplexing → gene expression quantification → ambient RNA removal → spatial barcode assignment → spatial analysis.

---

## Quick reference: launching the pipeline

`./slidr` is the single entry point for local, GCP and Slurm runs — it replaces the older, now-removed `gcp_runner.sh`/`slidr.sh` scripts.

```bash
./slidr --bcl <BCL_ID> [pipeline flags]                        # run locally
./slidr --bcl <BCL_ID> --gcp [GCP options] [pipeline flags]    # run on a GCP VM
./slidr --bcl <BCL_ID> --slurm [Slurm options] [pipeline flags] # submit as a Slurm job
```

Without `--gcp`, `./slidr` installs any missing local dependencies (`uv`, Miniforge/conda) and runs `uv run python workflow/main.py` directly in the current shell.

With `--gcp`, `./slidr` creates a GCE VM, passes your flags through as instance metadata, and streams its logs back via `watch_run.sh`.

**Before running with `--gcp`**, the following must be uploaded manually to the GCS prefix named by `paths.input_path` — the pipeline never puts them there for you:
- `config.yaml` — pipeline configuration
- `auth_key.json` — Google service account key, needed only when `settings.metadata_source` is a Google Sheet
- `<BCL_ID>/` — the sequencing run folder
- any further run folder named by a `Merge RNA From BCL` / `Merge Spatial From BCL` metadata column — see [Split-BCL runs](#split-bcl-runs)

There is no built-in default and no command-line equivalent: set `paths.input_path` in your local `config/config.yaml` to that `gs://` prefix. `./slidr` reads it locally and passes the location to the VM as instance metadata; the VM downloads its own config and key from there, and the pipeline then downloads the run folders themselves into `settings.gcs_download_dest`.

`paths.input_path` is dual-purpose, exactly like `reference_path`/`puck_path`/`raw_barcodes_path`: a local directory of run folders for a local run, and the bucket prefix they are staged out of when staging. Which one it is is decided by `--stage-gcs`/`--gcp`, never by the shape of the value — so a `gs://` value without the flag, or a local path with it, is a clear error rather than a silent misread.

### Minimal example

```bash
./slidr --bcl 20240101_RUNID --gcp --run-all
```

### Full example with GCP options

```bash
./slidr \
  --bcl 20240101_RUNID \
  --gcp \
  --project my-project-id \
  --zone us-central1-a \
  --machine n1-standard-16 \
  --disk 500GB \
  --gpu \
  --run-all
```

### Run only specific stages

```bash
# Cellbender only
./slidr --bcl 20240101_RUNID --gcp --cellbender

# Spatial analysis only (skip ambient removal)
./slidr --bcl 20240101_RUNID --gcp --spatial-count --spatial-analysis --no-cellbender

# Force re-run everything
./slidr --bcl 20240101_RUNID --gcp --run-all --force
```

---

## `./slidr` — all flags

### Required
| Flag | Description |
|---|---|
| `--bcl BCL_ID` | BCL run ID — **required** |

### Informational (exit immediately, need no other flag)
| Flag | Short | Description |
|---|---|---|
| `--version` | `-v` | Print the pipeline version (read from `pyproject.toml`) and exit |
| `--help` | `-h` | Print usage and exit |

### Where the run executes
| Flag | Short | Description |
|---|---|---|
| `--gcp` | `-gc` | Create a GCE VM, run there, and stream its logs back |
| `--slurm` | `-sl` | Submit as a Slurm batch job |

Neither takes a location. The GCS prefix inputs are staged from is `paths.input_path` in the local
`config/config.yaml` — there is no `--input` flag. `--gcp` always stages (a fresh VM has nothing else to
read); `--slurm` stages only when asked with `--stage-gcs`. Where the download lands is
`settings.gcs_download_dest`, which defaults into the run's own output tree.

### GCP VM options (only apply with `--gcp`)
| Flag | Default | Description |
|---|---|---|
| `--project PROJECT` | none — pass explicitly | GCP project ID |
| `--zone ZONE` | `us-central1-a` | GCP compute zone |
| `--machine TYPE` | `n1-standard-16` | Machine type |
| `--disk SIZE` | `200GB` | Boot disk size |
| `--gpu [TYPE]` | off; `nvidia-tesla-t4` when bare | Attach a GPU. Takes an optional GCE accelerator name (e.g. `--gpu nvidia-l4`), validated as lowercase alphanumerics/hyphens before the VM is created |

### Slurm options (only apply with `--slurm`)
| Flag | Default | Description |
|---|---|---|
| `--partition NAME` | cluster default | Partition to submit to |
| `--time LIMIT` | `24:00:00` | Job time limit |
| `--workdir PATH` | derived | Where staged inputs and outputs go. Resolves to `--workdir`, then the parent of `paths.output_path` in the config the job will run with, then `$SCRATCH/slidr`, then `<repo>/slidr-work/<BCL_ID>`. The chosen value and its origin are printed at submit time. Only valid with `--stage-gcs` |

CPUs and memory are not flags: they come from `settings.threads`/`settings.memory` in the config the job
will run with (fetched from the bucket at submit time for a staged run).

### Pipeline stage flags (passed through to `workflow/main.py`)
`./slidr` accepts the same short forms as `workflow/main.py` and always forwards the long form, so VM
metadata and the run summary record one canonical spelling regardless of what was typed.

| Flag | Short | Description |
|---|---|---|
| `--fastqs [PATH]` | `-fq` | Treat the input as already-demultiplexed FASTQs; skips mkfastq. PATH is optional — with no value the FASTQs are read from `paths.input_path`. Mutually exclusive with `--mkfastq` |
| `--metadata PATH` | `-md` | Metadata `.tsv`/`.csv` path or Google Sheet URL for this run, overriding `settings.metadata_source`. Recorded in the run summary as a `Metadata override` line |
| `--stage-gcs` | `-sg` | Fetch this run's inputs from GCS: `reference_path`/`puck_path`/`raw_barcodes_path`, plus `software_path` when it is a `gs://` URI, plus any extra BCL run folder a split-BCL run merges from, plus — under `--slurm` — the sequencing data and the auth key if the bucket has one. A `--slurm` run's config is never staged. Implied by `--gcp` |
| `--run-all` | `-ra` | Run all stages end-to-end (default if no stage flag given) |
| `--mkfastq` | `-mf` | Run Cellranger mkfastq only |
| `--count` | `-ct` | Run Cellranger count only |
| `--cellbender` | `-cb` | Run Cellbender only |
| `--no-cellbender` | `-nb` | Skip Cellbender when running the full pipeline |
| `--spatial-count` | `-sc` | Run spatial barcode counting only |
| `--spatial-analysis` | `-sa` | Run spatial analysis only |
| `--force` | `-f` | Overwrite existing outputs and re-run |

`workflow/main.py` additionally accepts `--config`/`-cf PATH` (a custom config location, used by the
Slurm payload and not exposed on `./slidr`). There is no `--reference` or `--generate-bam` flag anywhere
in the pipeline — the equivalent settings (`reference_genome`, `generate_bam`) are config-only (see
below). There is no `--gcs-key` or `--software` flag either; see the Slurm section for what replaced
them.

### Running from pre-demultiplexed FASTQs (`--fastqs`)

Which kind of input a run has is stated by the caller, never inferred from the contents of `input_path`. Passing `--fastqs` means "these reads are already demultiplexed"; omitting it means "`input_path` holds BCL run folders". Every input check keys off that flag rather than probing.

The directory is optional. A bare `--fastqs` means "the reads are this run's input folder", which is the
usual case; `--fastqs DIR` points at somewhere else. Below, `DIR` means whichever applies:

| | A bare `--fastqs` resolves to |
|---|---|
| local run | `paths.input_path` itself |
| `--stage-gcs`/`--gcp` | `<gcs_download_dest>/<BCL_ID>` — the local copy of `<input_path>/<BCL_ID>`, staged for you. The same folder a BCL run would be staged into, since the bucket prefix holds one folder per run either way |

That makes `--gcp --fastqs` need no path at all; before, you had to know the VM's own staging directory
and pass it by hand.

| Without `--fastqs` | With `--fastqs` |
|---|---|
| mkfastq runs (with `--mkfastq`/`--run-all`), and `input_path/<BCL_ID>` **must** validate as an Illumina run folder — a failure is a hard error | mkfastq is skipped; `DIR` must exist and contain `.fastq.gz` files (checked in `config.py` before any stage runs, or — for a staged bare `--fastqs`, whose directory does not exist yet — in `main.py` right after it is downloaded) |
| `count`/`spatial-count` read FASTQs from `output/mkfastq/` | `count`/`spatial-count` read FASTQs from `DIR`, and `helpers.missing_fastqs` first verifies an R1+R2 exists per sample per library, per declared lane |
| cellranger `--sample` is the sample name (mkfastq names libraries that way) | library names are read back off the filenames in `DIR` |
| `paths.input_path` is required | `paths.input_path` is required for a bare `--fastqs` (it *is* the input); with an explicit `DIR` it may be omitted and defaults to `DIR` |

FASTQs in `DIR` must be named `<library>_S<n>_L<lane>_<R1|R2|I1|I2>_<chunk>.fastq.gz` (the cellranger/bcl2fastq convention), with the library identifying the sample and its type:

```
mysample_GEX_S1_L001_R1_001.fastq.gz    # gene expression
mysample_SB_S2_L003_R1_001.fastq.gz     # spatial barcodes
mysample_S1_L001_R1_001.fastq.gz        # gene expression, bare form (as mkfastq names it)
mysample_sb_S2_L003_R1_001.fastq.gz     # spatial barcodes, as mkfastq names it
```

`-` also works as the separator, the mode token is case-insensitive, and trailing `_`/`-`-delimited extras are allowed. Matching is anchored at the sample name (`helpers.library_matches_sample`), so `mysample` never picks up `mysample_ln`'s FASTQs. Because `spatial_count.jl` selects files by substring-matching the sample ID inside a directory, `helpers.stage_spatial_fastqs` symlinks each sample's spatial library into `tmp/spatial_fastqs/<sample>/` first, so gene-expression reads can't be counted as spatial reads.

```bash
# run everything except demultiplexing, reading FASTQs from `paths.input_path`
./slidr --bcl 20240101_RUNID --fastqs

# the same, from a different directory
./slidr --bcl 20240101_RUNID --fastqs /path/to/fastqs

# count only
./slidr --bcl 20240101_RUNID --count --fastqs /path/to/fastqs

# on GCP, pass no path: the run folder under `paths.input_path` is staged and read automatically
./slidr --bcl 20240101_RUNID --gcp --fastqs

# an explicit DIR under --gcp is for FASTQs already on the VM image, not for staged ones
./slidr --bcl 20240101_RUNID --gcp --fastqs /opt/reference_fastqs
```

---

## Pipeline stages

| Order | Stage | Tool | Flag |
|---|---|---|---|
| 1 | BCL → FASTQ | Cellranger mkfastq + bcl2fastq | `--mkfastq` |
| 2 | Gene expression quantification | Cellranger count | `--count` |
| 3 | Ambient RNA removal | CellBender | `--cellbender` |
| 4 | Spatial barcode counting | Julia (spatial_count.jl) | `--spatial-count` |
| 5 | Spatial analysis | R / Seurat | `--spatial-analysis` |

Stages 3–5 can run in parallel after stage 2. Stage 5 depends on stage 4.

---

## Split-BCL runs

A library sequenced across more than one run is declared per sample in the metadata: `Merge RNA From
BCL` / `Merge Spatial From BCL` name the extra run folder, and `Add RNA Index` / `Add SB Index` the
index that library carries there. mkfastq then demultiplexes each BCL in turn, with a separate
samplesheet per run (`metadata/<BCL_ID>_samplesheet.csv`), and routes the extra libraries —
`<sample>_split_rna` / `<sample>_split_sb` — into the same `mkfastq/<sample>/` and
`mkfastq/<sample>_sb/` folders as the primary reads, so everything downstream sees one merged
library.

The extra run folders are ordinary inputs and are resolved exactly like the primary one, as
`input_path/<BCL>` (`helpers.resolve_bcl_dir`, which handles an `input_path` that points at the
primary run folder itself by looking beside it). Under `--gcp`/`--stage-gcs` they are staged out of
the `paths.input_path` prefix, the same one the primary BCL comes from — so a split-BCL remote run
needs every run folder uploaded there, each named as the metadata spells it.

Staging is the pipeline's job, not the launcher's. `slidr_gcp.sh`/`slidr_slurm.sh` know only the one
BCL ID `./slidr` handed them; which *other* BCLs a run needs is stated only by the metadata, which
nothing outside the pipeline reads. So neither script copies sequencing data at all any more —
`pipeline.stage_input_data` brings down every run folder, from inside `run_mkfastq`, at the first point
the answer is known:

- it is a no-op without `--stage-gcs`/`--gcp` — a missing run folder is then just a missing local
  input, reported by the run-folder check rather than silently downloaded
- it runs only if mkfastq itself runs, so a re-run whose FASTQs already exist never re-downloads
  hundreds of GB
- a run folder already present is left alone; `RESTAGE=1` in the environment forces a re-download
- each staged folder is then put through `validate_bcl_dir` — the primary included, since this is the
  only place all of them are known together

The same function stages the reads for a `--fastqs` run, called from `main.py` instead (see above):
one implementation, one already-staged rule, one destination.

`helpers.merge_bcls` is the single place the two merge columns are turned into a list of BCL IDs, so
the samplesheet writer, the staging step and the demultiplexing loop cannot disagree about which runs
those are. It drops blanks, this run's own BCL (there is no second run folder to fetch, and
`create_samplesheet` warns) and duplicates.

```bash
# every run folder the metadata names is fetched from `paths.input_path` automatically
./slidr --bcl 20240101_RUNID --gcp --run-all
./slidr --bcl 20240101_RUNID --slurm --stage-gcs --run-all
```

---

## config/config.yaml — field reference

Edit this file locally and stage it to GCS (as `config.yaml`) before launching with `--gcp`.

```yaml
paths:
  output_path:          # where outputs are written; the run lands in <output_path>/<BCL_ID>/, and a
                         # path that already ends in the BCL ID is taken as the run directory itself
                         # (so it is never nested twice)
  input_path:            # root dir containing BCL run folders (BCL_ID must be a subdirectory of this
                         # path), or -- under --stage-gcs/--gcp -- the gs:// prefix those folders are
                         # staged out of, which for --gcp also holds config.yaml and auth_key.json.
                         # Optional when running with --fastqs DIR, which then supplies the input instead
  software_path:         # directory scanned for Cellranger/bcl2fastq executables. May instead be a
                         # full gs:// prefix, which is downloaded to output/software/ on first use;
                         # that form requires --stage-gcs, but a local directory stays valid with it
  raw_barcodes_path:     # raw puck barcode files (BeadBarcodes.txt / BeadLocations.txt); a gs:// path when running with --gcp
  puck_path:             # processed puck CSV files (generated from raw_barcodes_path if absent); a gs:// path when running with --gcp
                         # required only by the spatial-count stage, which Flex runs never reach, so an
                         # unset value is reported there rather than at startup
  reference_path:        # directory containing reference genome subdirectories; a gs:// path (or bare bucket name) when running with --gcp
  auth_key_path:         # Google service account JSON key. Required only when metadata_source is a
                         # Google Sheet -- a local .tsv/.csv run never reads it. May be a local path or
                         # a gs:// object, downloaded to the run's tmp/ (mode 600) when needed.
                         # `gcloud auth login` is NOT a substitute: the file is read directly, not
                         # through Application Default Credentials

settings:
  memory: 50             # GB allocated to Cellranger
  threads: 16            # cores for Cellranger and Julia
  metadata_source:       # Google Sheet URL (https://docs.google.com/spreadsheets/d/...#gid=...)
                         # OR absolute path to a .tsv/.csv metadata file
                         # the worksheet tab is taken from the URL's gid — no separate worksheet-name field
  gcs_download_dest:     # local dir a staged run downloads the `paths.input_path` run folders into.
                         # Unset, it defaults to <output_path>/<BCL_ID>/data -- inside the run's own
                         # directory, so a staged run is self-contained and needs no per-machine setting,
                         # but a sibling of output/ rather than inside it, since output/ is what gets
                         # uploaded to output_bucket. Ignored by a local run.
                         # Overridable per host with SLIDR_GCS_DOWNLOAD_DEST
  output_bucket:         # GCS bucket/prefix outputs are uploaded to after a run completes
  reference_genome:      # reference genome directory name under reference_path
  alerts: false          # send Slack alerts on errors/completion. Every alert is tagged with the
                         # run's BCL ID as a bold first line, so concurrent runs are distinguishable
  slack_token:           # where to find the Slack bot token (required if alerts: true): a local file
                         # path, a gs:// object (read with `gcloud storage cat`, never written to disk),
                         # or the literal token. A file or object is preferred -- a literal token
                         # travels with this file as it is copied and staged to buckets

workflow:
  generate_bam: false
  cellbender_total_droplets:       # set null for auto-detection from barcode count
  cellbender_estimated_cells:      # set null to let CellBender estimate
  cellbender_epochs: 160
  cellbender_learn_rate: 0.5
  spatial_downsampling:            # optional float; downsamples spatial reads before spatial_count.jl
  top_n_percent_umi_filter:        # optional float 0-100; bead UMI filtering percentile
  flex_emptydrops_minimum_umis: 100
  flex_probe_set: /path/to/probe_set.csv
  flex_spatial_R1_path: /path/to/flex/fastqs/R1
  flex_spatial_R2_path: /path/to/flex/fastqs/R2
  flex_gex_fastqs:                 # list of GEX FASTQ prefixes
    - fastq_prefix_1
```

There is no separate `cloud:` section — bucket locations for the sequencing data, reference genome, puck files, and raw barcodes are the same `paths.*` fields used locally, just pointed at a `gs://` location when the run stages from GCS (`--gcp` or `--stage-gcs`). Similarly `settings.output_bucket` replaces what used to be `cloud.output_bucket`, etc.

Whether `input_path`, `reference_path`, `puck_path` and `raw_barcodes_path` are read locally or staged from GCS is decided by the flag, never by the shape of the value:

| | Without `--stage-gcs`/`--gcp` | With `--stage-gcs`/`--gcp` |
|---|---|---|
| Expected value | a local directory, validated on disk at startup | a `gs://bucket/prefix` (or bare `bucket/prefix`) string; only the form is checked at startup |
| `input_path` | used in place | `<BCL>/` downloaded to `settings.gcs_download_dest`, for the primary BCL and every merged-in one |
| `reference_path` | used in place | `<genome>` downloaded to `output/reference/` |
| `puck_path` | read and written in place | `<puck>.csv` downloaded to `output/pucks/`, which is also where generated pucks are written |
| `raw_barcodes_path` | read in place | `<puck>/` downloaded to `output/barcodes/` |

Mixing the two forms is caught rather than misread: a `gs://` value without the flag is an error naming
`--stage-gcs`, and an absolute local path *with* it is an error too (`config.looks_local`). The latter
matters because a bucket name cannot begin with `/`, so `gcs_uri` would otherwise quietly turn
`/data/runs` into `gs:///data/runs` and fail much later as an unreadable bucket.

`paths.software_path` is stageable too, but by a deliberately narrower rule — it is the one field where
both forms remain valid under the flag:

| | Value is a local directory | Value is a full `gs://…` URI |
|---|---|---|
| Without `--stage-gcs` | scanned in place | hard error, naming `--stage-gcs` |
| With `--stage-gcs`/`--gcp` | scanned in place — *not* an error | contents downloaded to `output/software/`, then scanned there |

The asymmetry is intentional: Cellranger is usually preinstalled on the machine that runs the pipeline
(the `--gcp` VM image is built that way), and since `--gcp` implies `--stage-gcs`, requiring a bucket
would break every existing VM run. A bare `bucket/prefix` is *not* accepted here, unlike the three
fields above — it is indistinguishable from a relative local directory, and guessing wrong costs a
multi-GB download rather than an error message. Staging is lazy (a warm `software_cache.txt` skips it),
happens at most once per run, is reused by later runs in the same output tree, and ends by restoring
execute bits, which `gcloud storage cp` does not preserve.

A `gs://` value without the flag is a clear error ("not a directory... pass `--stage-gcs`"), not a silent misread.

---

## Metadata format

Metadata must be provided as a Google Sheet or `.tsv`/`.csv` file with these columns (enforced by `workflow/pipeline/config.py`'s `METADATA_SPECS` / `workflow/pipeline/helpers.py`'s `METADATA_OPT`):

| Column | Description |
|---|---|
| `Run` | `YES` to include this sample, anything else to skip |
| `Email` | User email (used for Slack notifications) |
| `Sample Name` | Unique sample identifier |
| `BCL` | BCL run ID (matched to `--bcl`); the BCL directory is resolved as `input_path/<BCL>` — there is no separate "BCL Path" column |
| `Species` | `Human` or `Mouse` (auto-selects reference genome) |
| `Chemistry` | e.g. `3Pv3`, `5P`, `Flex` |
| `RNA Index` | Cellranger index for the RNA library |
| `Lane` | Sequencing lane(s), comma-separated |
| `SB Index` | Cellranger index for the spatial barcode library |
| `SB Lane` | Lane(s) for spatial barcode library, comma-separated |
| `Puck ID` | Identifier matching a puck CSV in `puck_path` |
| `Merge RNA From BCL` | BCL ID to merge RNA FASTQs from (optional). Its run folder is demultiplexed alongside the primary one, and staged from the `paths.input_path` prefix on a `--gcp`/`--stage-gcs` run — see [Split-BCL runs](#split-bcl-runs) |
| `Merge Spatial From BCL` | BCL ID to merge spatial FASTQs from (optional), staged and demultiplexed the same way |
| `Add RNA Index` | Index for merged RNA library (optional) |
| `Add SB Index` | Index for merged spatial library (optional) |
| `Add Puck ID` | Puck override for a merged sample (optional) |
| `Cellranger` | Per-sample Cellranger version override (optional). Accepts `V8`/`v8`/`8` shorthand, expanded to a full release via `helpers.CELLRANGER_VERSIONS`, or a full version (`8.0.1`) used as written. Defaults to `8.0.1` (`9.0.1` for Flex, which needs `cellranger multi`). Honoured per sample by `count`; mkfastq is one invocation for the whole run, so it uses the single declared version and warns if samples disagree |
| `Flex Probe Barcode IDs` | Probe barcode IDs for Flex chemistry (required for `Flex`). Separated by `,` or `\|`, or a mixture — parsed by `helpers.split_probe_barcodes` |

`example_metadata.tsv` in the repository root carries this exact schema — all 11 required columns followed by the 7 optional ones — with rows demonstrating a two-sample 3' run, a skipped row (`Run: NO`), a cross-BCL RNA merge, a `Cellranger` version override and a Flex pool. Optional columns left empty are ignored, so it can be edited down in place.

---

## GCP architecture

```
Local machine                  GCP
─────────────────              ──────────────────────────────────────────
./slidr --gcp          ──→    creates VM from image slidr-vm-image6, passing
  (writes .last_run)           input_path / bcl-id / flags as instance metadata
watch_run.sh           ←──    VM runs workflow/bash/slidr_gcp.sh as startup script:
  (reads .last_run,             1. sets up the ops-agent so runtime.log reaches Cloud Logging
   streams Cloud Logging)       2. drops privileges to the `runner` user
                                3. clones the `stable` branch of the repo
                                4. downloads config.yaml + auth_key.json from input_path
                                5. uv sync
                                6. uv run python workflow/main.py --bcl ... <flags>
                                     main.py stages the run folders themselves out of
                                     paths.input_path, once the metadata says which
                                7. gcloud storage cp results → settings.output_bucket
                                8. self-deletes VM once idle (or on failure)
```

Key GCP resources:
- **Service account**: `slidr-runner@<project>.iam.gserviceaccount.com`
- **VM image**: `slidr-vm-image6`
- **Network**: `vpc1` / subnet `my-subnet-central`
- **VM scopes**: `storage-rw,compute-rw,logging-write,monitoring-write` (narrowed from `cloud-platform`; actual permissions are still governed by the service account's IAM roles)
- **VM paths**: outputs at `/pipeline/out/<BCL_ID>/output`, staged sequencing data at `/pipeline/out/<BCL_ID>/data/<BCL_ID>` (there is no separate `/pipeline/data` any more — size `--disk` for data **and** outputs together)

There is no gcsfuse-mounted data bucket anymore — BCL data, reference genomes, and puck files are each copied on demand from the `gs://` paths in `config.yaml` via `gcloud storage cp`.

---

## Running on Slurm

```bash
# on a cluster that already sees the lab filesystem — everything read locally
./slidr --bcl 20240101_RUNID --slurm --run-all

# on an independent cluster — the job downloads its data from GCS, but still
# runs with this checkout's config/config.yaml
./slidr --bcl 20240101_RUNID --slurm --stage-gcs --run-all
```

Passing `--stage-gcs` is what makes a Slurm run self-staging; without it nothing is downloaded and the behaviour is unchanged. The *location* always comes from `paths.input_path`; the flag states only whether to read it as a bucket prefix or as a local directory, so a configured bucket alone never silently turns a run into a staging one.

### What a staged Slurm run does

`./slidr --slurm --stage-gcs` submits `workflow/bash/slidr_slurm.sh` (the Slurm counterpart of `slidr_gcp.sh`) instead of invoking `main.py` directly:

```
Submit node                       Compute node (inside the batch job)
─────────────────                 ──────────────────────────────────────────────
./slidr --slurm --stage-gcs ─→   workflow/bash/slidr_slurm.sh:
  reads config/config.yaml          1. checks gcloud is present and authenticated
  to size the allocation and        2. downloads auth_key.json → workdir, if the
  derive the workdir                   bucket has one (absence is not an error)
  ensure_uv / ensure_conda          3. exports SLIDR_OUTPUT_PATH=<workdir>/outs
  uv sync                           4. uv sync
  sbatch slidr_slurm.sh             5. main.py --stage-gcs …  (no --config: the
                                       checkout's config/config.yaml is used)
                                    6. main.py stages every run folder the metadata
                                       names out of paths.input_path, plus the
                                       reference genome, puck maps, raw barcodes and
                                       any gs:// software_path
                                    7. main.py uploads outputs → settings.output_bucket
```

Differences from the GCP path, all because there is no VM to own: the repo is not cloned (the job runs from the checkout it was submitted from), privileges are not dropped (Slurm already runs as the submitting user), and nothing self-deletes. Staging happens inside the job, not at submit time, so hundreds of GB never move through a shared login node. Sequencing data already staged is not re-downloaded — set `RESTAGE=1` to force it.

The script itself no longer downloads the sequencing data: `paths.input_path` names the bucket prefix, and `main.py` fetches from it into `<workdir>/outs/<BCL_ID>/data`. That is what lets a split-BCL job fetch the run folders its metadata merges from — the script only ever knew the one BCL ID it was handed.

**A Slurm run never stages `config.yaml`.** It runs with `config/config.yaml` from the checkout it was submitted from, and the bucket needs no copy. The GCP path still downloads one because a fresh VM has no checkout to read; a cluster job does, and that file is the one the submitter just edited — staging a copy written for another machine on top of it only let the two disagree. Two follow-on effects: a local config change takes effect on the next `sbatch` with nothing to re-upload, and the CPUs/memory requested from Slurm are read from the same file the job runs with, so the allocation cannot be sized off a different config than the run uses.

| Flag | Applies to | Description |
|---|---|---|
| `--stage-gcs` | `--slurm` | Opt a Slurm run into reading `input_path` (and its reference/puck/barcode/software paths) as GCS locations and staging them. The key is staged if the bucket holds one; the config never is. Implied by `--gcp` |
| `--workdir PATH` | staged `--slurm` | Where outputs are written, and — by default — where inputs are staged beneath them. Rarely needed: defaults to the parent of `paths.output_path` in the config, then `$SCRATCH/slidr`, then `<repo>/slidr-work/<BCL_ID>`. The resolved value and its origin are printed at submit time |

### Pointing a config at the job's workdir

A staged Slurm job writes into a workdir that `config/config.yaml` cannot know at the time it was written, so those locations are overridden from the environment (`config.py`'s `PATH_ENV_OVERRIDES` and `SETTINGS_ENV_OVERRIDES`) rather than by rewriting the YAML:

| Config field | Environment variable | Set by `slidr_slurm.sh` to |
|---|---|---|
| `paths.output_path` | `SLIDR_OUTPUT_PATH` | `<workdir>/outs` |
| `settings.gcs_download_dest` | `SLIDR_GCS_DOWNLOAD_DEST` | not set by the script — the default puts staged data at `<workdir>/outs/<BCL_ID>/data`. Export it to stage onto a different filesystem |
| `paths.auth_key_path` | `SLIDR_AUTH_KEY_PATH` | `<workdir>/auth_key.json` — **only if** a key was staged; otherwise the config's own value stands |
| `paths.software_path` | `SLIDR_SOFTWARE_PATH` | not set by the script — export it yourself, or point the config at a `gs://` prefix |
| `paths.input_path` | `SLIDR_INPUT_PATH` | **not set by the script, and must not be**: under `--stage-gcs` this field is the `gs://` prefix to stage *from*, so pointing it at a local directory would leave the job with nothing to download. Where data lands is `gcs_download_dest`, above |

An override takes precedence over the file, is still validated normally, and is recorded in the run's summary log. Everything else — threads, memory, buckets, reference/puck/barcode locations, metadata source — comes from `config/config.yaml` as committed. The overrides are exported into the job's environment rather than written into the config, so concurrent jobs from one checkout never fight over that file.

Prerequisites on the cluster: `gcloud` on `PATH` (load the module first if yours provides one) and working credentials — either `gcloud auth login`, or `gcloud auth activate-service-account --key-file=...` for a headless cluster. Both persist to `~/.config/gcloud`; if that is not shared with compute nodes, set `CLOUDSDK_CONFIG` to somewhere that is. There is no `--gcs-key` flag: slidr does not re-implement what `gcloud auth` already does persistently.

---

## Monitoring a running job

After `./slidr --gcp` launches the VM, it automatically streams logs via `watch_run.sh`. To reattach later:

```bash
# Reattach to the most recent run
./watch_run.sh

# Reattach to a specific VM
./watch_run.sh slidr-run-1719123456-12345

# ... one in a different project than the last run's
./watch_run.sh slidr-run-1719123456-12345 my-other-project
```

`.last_run` records the most recent VM name on line 1 and the project it was created in on line 2,
and `watch_run.sh` reads both. No project ID is hardcoded anywhere in the repo: the project is
resolved from the second argument, then `$SLIDR_PROJECT`, then `.last_run`, then
`gcloud config get-value project`, and the script errors with those four options listed if all are
empty. A `.last_run` written before this (VM name only) still works — it just supplies no project.

---

## Output structure

```
<output_path>/<BCL_ID>/
├── log/
│   ├── runtime.log                 # main execution log
│   ├── <timestamp>_summary.log     # per-run summary of settings and result
│   └── <stage>.log                 # mkfastq / count / cellbender / spatial_barcodes /
│                                   #   spatial_analysis / takara_pipeline (Flex)
├── metadata/
│   ├── metadata_summary.csv        # the metadata rows selected for this run
│   ├── samplesheet.csv             # cellranger mkfastq sample sheet
│   └── multi_samplesheet.csv       # Flex only: cellranger multi sample sheet
├── output/
│   ├── mkfastq/<library>/          # demultiplexed FASTQs (<sample> and <sample>_sb)
│   ├── count/<sample>/             # Cellranger count outputs (Flex: count/flex/<sample>_<BC>/)
│   ├── cellbender/<sample>/        # CellBender filtered matrices
│   ├── spatial_barcodes/<sample>/  # SBcounts.h5
│   ├── spatial_analysis/<sample>/  # Seurat objects, coords, summary PDFs
│   ├── flex/                       # Flex only: demux/ pucks/ samplesheets/ trekker/
│   ├── reference/                  # staged runs only: local copy of the reference genome
│   ├── pucks/                      # staged runs only: puck CSVs, downloaded or generated
│   ├── barcodes/                   # staged runs only: raw bead barcodes
│   └── software/                   # only when software_path is a gs:// prefix: Cellranger/bcl2fastq
├── data/<BCL>/                     # staged runs only: run folders downloaded from paths.input_path,
│                                   #   one per BCL the run demultiplexes. NOT under output/ — see below
└── tmp/                            # scratch; also holds auth/auth_key.json for a gs:// key.
                                    #   Safe to delete after a run
```

No single run contains every entry: the `flex/`, `multi_samplesheet.csv` and `takara_pipeline.log`
entries appear only for Flex chemistry, `data/`/`reference/`/`pucks/`/`barcodes/` only when staging from
GCS, and `software/` only when `paths.software_path` is a `gs://` prefix.

`data/` is the only staged directory that sits outside `output/`, and the reason is the upload: on
success `main.py` copies the whole of `output/` to `settings.output_bucket`, and `data/` holds raw inputs
that came *out* of a bucket — hundreds of GB per run. Staging them under `output/` would send the entire
sequencing run straight back to GCS every time a run finished. Keeping them a sibling means the run is
still self-contained (nothing is written outside its own directory, which is all a fresh VM or a cluster
job can count on) without that cost. `settings.gcs_download_dest` moves them off this filesystem
entirely, e.g. onto cluster scratch.

### Key output files per sample

| File | Stage | Description |
|---|---|---|
| `count/<sample>/filtered_feature_bc_matrix.h5` | count | Cellranger expression matrix |
| `cellbender/<sample>/cellbender_output_filtered.h5` | cellbender | Ambient-corrected matrix |
| `spatial_barcodes/<sample>/SBcounts.h5` | spatial-count | Barcode → coordinate mapping |
| `spatial_analysis/<sample>/seurat.qs` | spatial-analysis | Seurat object with spatial embeddings |
| `spatial_analysis/<sample>/coords.csv` | spatial-analysis | Cell coordinates |
| `spatial_analysis/<sample>/summary.pdf` | spatial-analysis | QC summary plots |

---

## Key files

| File | Purpose |
|---|---|
| `slidr` | Unified entry point — runs locally, creates a GCP VM and streams its logs with `--gcp`, or submits a Slurm job with `--slurm` |
| `workflow/bash/slidr_gcp.sh` | VM startup script (runs inside the VM) |
| `workflow/bash/slidr_slurm.sh` | Slurm batch payload: stages inputs from GCS, then runs the pipeline (runs on the compute node) |
| `watch_run.sh` | Re-attaches to a running VM's log stream |
| `config/config.yaml` | Pipeline configuration (must be staged to GCS as `config.yaml` before GCP runs) |
| `workflow/main.py` | Top-level pipeline orchestration |
| `workflow/pipeline/config.py` | Argument parsing + config loading |
| `workflow/pipeline/pipeline.py` | Stage implementations (mkfastq, count, cellbender, spatial) |
| `workflow/pipeline/helpers.py` | Utilities: metadata loading, software detection, logging |
| `pyproject.toml` | Python dependencies (managed by uv) |
| `envs/conda.yml` | Conda environment for R/CellBender |
| `envs/julia/` | Julia project/manifest for spatial barcode counting |
| `.last_run` | Most recent run's VM name (line 1) and GCP project (line 2), auto-written by `./slidr --gcp`; gitignored |
| `software_cache.txt` | Cached paths to Cellranger/bcl2fastq/Julia executables (repo root) |

---

## Software cache

On first run, slidr scans `software_path` for Cellranger, bcl2fastq, and Julia, and caches paths in `software_cache.txt` (repo root). To pre-populate or pin specific versions:

```
/path/to/cellranger-8.0.1/bin/cellranger
/path/to/bcl2fastq2/bin/bcl2fastq
/home/runner/.juliaup/bin/julia
```

---

## Troubleshooting

- **Cellranger refuses to run on re-run**: manually delete `output/mkfastq/` and/or `output/count/` — Cellranger will not overwrite existing outputs. Use `--force` to have slidr handle this.
- **`... is not a usable Illumina BCL run directory`**: a run folder failed the run-folder check; the preceding log lines list exactly what's missing. If the reads are in fact already demultiplexed, pass them with `--fastqs` instead — the pipeline will not guess this for you. The check covers every BCL the run demultiplexes, so for a merged-in one the fix is usually that the `Merge ... From BCL` value does not match the folder name under `paths.input_path`. On a staged run, a partial download is re-fetched with `RESTAGE=1`.
- **`input_path must be a GCS location when staging from GCS`**: `--stage-gcs`/`--gcp` reads `paths.input_path` as a bucket prefix, and this one is a local path. Point it at the `gs://` prefix holding the run folders; where they are downloaded *to* is `settings.gcs_download_dest`, not this field. The same error names `reference_path`/`puck_path`/`raw_barcodes_path` for the same mistake.
- **`input_path is a GCS location, but this run was not asked to stage from GCS`**: the converse — add `--stage-gcs` (or `--gcp`), or point the field at a local directory.
- **`no .fastq.gz files were staged to the FASTQ directory for --fastqs`**: a staged bare `--fastqs` downloaded `<input_path>/<BCL_ID>` and found no FASTQs in it. Usually that folder holds BCLs, in which case drop `--fastqs`.
- **`the following FASTQ files are missing from ...`**: a `--fastqs` directory is incomplete or its filenames don't follow the convention above. The message names each sample, library and read that couldn't be found; check the `Lane`/`SB Lane` metadata columns against the `_L00N_` tokens in the filenames.
- **Missing reference genome when staging**: set `paths.reference_path` to a `gs://` location (or a bare bucket path — `gs://` is prepended automatically) — the pipeline will `gcloud storage cp` the `settings.reference_genome` subdirectory from there if it isn't already staged under `output/reference/`.
- **Puck file not found**: when staging, `paths.puck_path` must be a `gs://` location — the pipeline downloads `<puck_id>.csv` from it; if that object is absent it falls back to copying the corresponding raw barcodes from `paths.raw_barcodes_path` and generating the puck locally, and only errors if `raw_barcodes_path` is unset too.
- **`could not stage ... from GCS`**: `gcloud`'s own stderr is echoed under the message, along with the exact command that failed. Check the object exists (`gcloud storage ls <uri>`) and that the active account can read the bucket (`gcloud auth print-access-token`).
- **`gcloud is not on PATH` on a cluster**: many clusters ship the CLI as a module — `module load google-cloud-sdk` (or equivalent) before submitting. `./slidr --slurm --stage-gcs` also checks this at submit time so the job doesn't fail minutes later.
- **Google Sheets auth failure**: verify `auth_key_path` points to a valid *service account* JSON key (an OAuth client-secret file will not work), and that the sheet is shared with the key's `client_email`. `gcloud auth login` is not a substitute — the file is read directly. A `gs://` value is downloaded to the run's `tmp/auth/` first. If the metadata is a local `.tsv`/`.csv`, the key is not needed at all.
- **Slack alerts not arriving**: `settings.slack_token` may be a local file path, a `gs://` object or the literal token; a path-shaped value that does not resolve is a hard error rather than being tried as a token. The bot needs `users:read.email` (to map the `Email` column to a user) and `chat:write`. A token file that is group- or world-readable produces a warning.
- **Cellranger version not found**: `test_and_install_software` matches an install directory named `cellranger[-_v]<version>` exactly, so a `V8` in the `Cellranger` column is expanded via `helpers.CELLRANGER_VERSIONS`. If your site has a different patch release, edit that table — an unmapped value is passed through as written. The separator between name and version may be any run of `-`, `_` or `v` (`cellranger-8.0.1`, `cellranger_8.0.1`, `cellranger-v8.0.1`), but there must be at least one: a bare `cellranger8.0.1` does not match. The match is exact at both ends, so `cellranger-8.0.11` and `cellranger-8.0.1-beta` are not picked up for `8.0.1`.
- **Software staged from GCS but nothing found**: the scan looks under `output/software/` once `paths.software_path` is a `gs://` prefix. Check the objects exist (`gcloud storage ls <software_path>`), and that the prefix holds the install *directories* themselves — the contents of the prefix are copied in, so `<software_path>/cellranger-8.0.1/bin/cellranger` becomes `output/software/cellranger-8.0.1/bin/cellranger`. Execute bits are restored automatically after the download, since `gcloud storage cp` does not preserve them.
