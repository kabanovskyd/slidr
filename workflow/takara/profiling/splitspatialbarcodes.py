#!/usr/bin/env python

import pandas as pd
import sys
import os.path
import math


def die(message: str, *hints: str) -> None:
    """
    Report a fatal error with actionable troubleshooting bullets, then exit non-zero.

    This script is one link in the Trekker/Flex chain, and its stdout is captured into
    takara_pipeline.log; the caller only ever sees the exit status. So a failure has to explain
    itself here, and it has to exit *non-zero* -- a bare sys.exit() returns 0 and would let the
    pipeline carry on as though the barcode files had been written.
    """

    print(f"[ERROR]: {message}", flush=True)
    if hints:
        print("Troubleshooting:", flush=True)
        for hint in hints:
            print(f" • {hint}", flush=True)
    sys.exit(1)


# every failure below is some form of "this puck whitelist is not shaped the way we need", so the
# same two escape hatches apply throughout
WHITELIST_HINTS = (
    "The whitelist is the puck's bead barcode file, downloaded from the vendor by the Flex path",
    "Check the `Puck ID` metadata column names the puck this sample was actually run on",
    "Delete the cached copy under the run's flex/pucks/ directory to force a fresh download",
)

print(sys.argv)
print(len(sys.argv))

if len(sys.argv) != 3:
    die(
        f"expected 2 arguments (spatial barcode whitelist, output directory), got {len(sys.argv) - 1}",
        "This script is normally invoked by the Flex spatial-profiling stage, not by hand",
        "Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`",
        f"Received: {sys.argv[1:]}",
    )
else:
    barcode_path=sys.argv[1]
    barcode_outdir=sys.argv[2]

try:
    whitelist=pd.read_csv(barcode_path, header=None, index_col=0, dtype=str).apply(lambda col: col.str.strip())
except pd.errors.EmptyDataError:
    die(
        f"{barcode_path} is empty",
        "An empty file usually means the vendor download failed or was interrupted",
        *WHITELIST_HINTS,
    )

# read with index_col=0, so the raw file's first column becomes the (discarded) index and the
# barcode/x/y are iloc columns 0/1/2 below — i.e. the file needs 4 raw columns (id, barcode, x, y).
# whitelist.shape[1] is the post-index count, so raw column count = whitelist.shape[1] + 1.
if whitelist.shape[1] < 3:
    die(
        f"{barcode_path} must have at least 4 columns (id, barcode, x, y), found {whitelist.shape[1] + 1}",
        "Columns are positional, so a header row or a different column order also fails this check",
        f"Inspect the file: `head -3 {barcode_path}`",
        *WHITELIST_HINTS,
    )

#BARCODES
barcodes=whitelist.iloc[:,0]
barcode_lengths = barcodes.str.len()
if barcode_lengths.nunique() > 1:
    die(
        f"{barcode_path} contains barcodes of differing lengths: {sorted(barcode_lengths.unique())}",
        "All bead barcodes must be the same length; they are split base-by-base into BeadBarcodes.txt below",
        "A whitelist concatenated from two pucks is the usual cause",
        *WHITELIST_HINTS,
    )
barcodes_bybase=pd.DataFrame([list(x) for x in list(barcodes)])
barcodes_bybase.to_csv(os.path.join(barcode_outdir,"BeadBarcodes.txt"), header=None, index=None)

#BARCODE LOCATIONS
x=whitelist.iloc[:,1]
y=whitelist.iloc[:,2]
if x.isna().any() or y.isna().any():
    die(
        f"{barcode_path} contains missing x or y coordinate values",
        "Every bead needs both coordinates, or it cannot be placed on the puck",
        f"Find the offending rows: `awk -F, 'NF<4 || $3==\"\" || $4==\"\"' {barcode_path} | head`",
        *WHITELIST_HINTS,
    )
try:
    x_float = [float(i) for i in x]
    y_float = [float(i) for i in y]
except ValueError as exc:
    die(
        f"{barcode_path} contains non-numeric x or y coordinate values: {exc}",
        "Columns 3 and 4 must be plain numbers; a stray unit suffix or thousands separator fails here",
        f"Inspect the file: `head -3 {barcode_path}`",
        *WHITELIST_HINTS,
    )
# float() accepts the literal text "nan"/"inf"/"-inf" without raising, and isna() above is False
# for those strings (the column is dtype=str), so guard finiteness explicitly — otherwise a
# garbled coordinate cell writes a NaN/Inf bead position straight through to BeadLocations.txt
if not all(math.isfinite(v) for v in x_float) or not all(math.isfinite(v) for v in y_float):
    die(
        f"{barcode_path} contains non-finite (NaN/Inf) x or y coordinate values",
        "The literal text 'nan' or 'inf' in a coordinate cell reaches this point without raising, so it is rejected here",
        f"Find them: `grep -in 'nan\\|inf' {barcode_path} | head`",
        *WHITELIST_HINTS,
    )
x=','.join([str(i) for i in x])
y=','.join([str(i) for i in y])

location_out = open(os.path.join(barcode_outdir,"BeadLocations.txt"), "w")
location_out.write(x+"\n")
location_out.write(y+"\n")
location_out.close()
