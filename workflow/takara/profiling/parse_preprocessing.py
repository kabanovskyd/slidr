#!/usr/bin/env python


import pandas as pd
import subprocess
import io
from scipy.io import mmread, mmwrite
import gzip
import os
import sys


#==================== Function to transpose matrix ====================
def transpose_matrix(mtx_in, mtx_out):
   """Tanspose the input matrix in R"""
   print("Transposing matrix...")
   if mtx_in.endswith(".gz"):
      with gzip.open(mtx_in, "rt") as infile:
         matrix = mmread(infile)
   else:
      matrix = mmread(mtx_in)
   transposed = matrix.transpose()
   with gzip.open(mtx_out, "wb") as f:
      mmwrite(f, transposed)
   return mtx_out


#============== Function to translate cell barcodes ===================
def translate_cellbarcodes(barcodes_in, barcodes_out, parseCBDict):
   """Translate the input Cell barcodes"""
   print("Translating cell barcodes...")
   dict_R1dT, dict_R1hex, dict_R3, dict_R2 = dict(zip(parseCBDict.Well, parseCBDict.R1_dT)), dict(zip(parseCBDict.Well, parseCBDict.R1_hex)), dict(zip(parseCBDict.Well, parseCBDict.R3)), dict(zip(parseCBDict.Well, parseCBDict.R2))
   barcodes_in["bc1_seqdT"]=barcodes_in['bc1_well'].map(dict_R1dT)
   barcodes_in["bc2_seq"]=barcodes_in['bc2_well'].map(dict_R2)
   barcodes_in["bc3_seq"]=barcodes_in['bc3_well'].map(dict_R3)
   for col, well_col in [("bc1_seqdT", "bc1_well"), ("bc2_seq", "bc2_well"), ("bc3_seq", "bc3_well")]:
      unmapped = barcodes_in.loc[barcodes_in[col].isna(), well_col].unique()
      if len(unmapped) > 0:
         raise ValueError(
            f"translate_cellbarcodes: well ID(s) {sorted(unmapped)} in column '{well_col}' not found in parseCBDict\n"
            "Troubleshooting:\n"
            " \u2022 Well IDs in the single-cell cell_metadata.csv are translated to barcode sequences via parseCBDict.csv\n"
            " \u2022 An unknown well ID means the two came from different Parse kit versions\n"
            " \u2022 Check parseCBDict.csv in the Trekker profiling script directory matches the kit used for this library"
         )
   barcodes_in['bc']=barcodes_in["bc3_seq"]+barcodes_in["bc2_seq"]+barcodes_in["bc1_seqdT"]
   barcodes_in['bc'].to_csv(barcodes_out, header=False, index=False, compression='gzip')
   return barcodes_out


#============== Function to generate features ===================
def generate_features(features_in, features_out):
   """Upadate the input features"""
   print("Generating features...")
   open_fun = gzip.open if features_in.endswith(".gz") else open
   with open_fun(features_in, 'rt', encoding="utf-8") as infile:
      features = infile.readlines()
   with gzip.open(features_out, "wt", encoding="utf-8") as outfile:
      for feature in features[1:]:
         outfile.write(feature.replace(',', '\t'))
   return features_out


def main():
   if len(sys.argv) < 3:
      print(f"[ERROR]: expected 2 arguments (matrix folder, script folder), got {len(sys.argv) - 1}", flush=True)
      print("Troubleshooting:", flush=True)
      print(" \u2022 This script is normally invoked by the Flex spatial-profiling stage, not by hand", flush=True)
      print(" \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`", flush=True)
      print(f" \u2022 Received: {sys.argv[1:]}", flush=True)
      sys.exit(1)
   sc_outdir = sys.argv[1]
   script_dir = sys.argv[2]
   #========================== Validate Input Arguments ==========================  
   parseCBDict = pd.read_csv(os.path.join(script_dir, "parseCBDict.csv"))
   non_gz_files = ['count_matrix.mtx', 'cell_metadata.csv', 'all_genes.csv']
   gz_files = ['count_matrix.mtx.gz', 'cell_metadata.csv.gz', 'all_genes.csv.gz']
   if all(os.path.isfile(os.path.join(sc_outdir, file)) for file in non_gz_files):
      print("Found all non-gzipped files")
      barcodes_file = pd.read_csv(os.path.join(sc_outdir, "cell_metadata.csv"))
      features_file = os.path.join(sc_outdir, "all_genes.csv")
      mtx_file = os.path.join(sc_outdir, "count_matrix.mtx")
   elif all(os.path.isfile(os.path.join(sc_outdir, file)) for file in gz_files):
      print("Found all gzipped files")
      barcodes_file = pd.read_csv(os.path.join(sc_outdir, "cell_metadata.csv.gz"))
      features_file = os.path.join(sc_outdir, "all_genes.csv.gz")
      mtx_file = os.path.join(sc_outdir, "count_matrix.mtx.gz")
   else:
      missing_non_gz = [f for f in non_gz_files if not os.path.isfile(os.path.join(sc_outdir, f))]
      missing_gz = [f for f in gz_files if not os.path.isfile(os.path.join(sc_outdir, f))]
      print(f"[ERROR]: required single-cell matrix files are missing in '{sc_outdir}'", flush=True)
      print(f"  Missing non-gzipped files: {missing_non_gz}", flush=True)
      print(f"  Missing gzipped files: {missing_gz}", flush=True)
      print("Troubleshooting:", flush=True)
      print("  \u2022 All three of count_matrix.mtx, cell_metadata.csv and all_genes.csv must be present, either all gzipped or all plain", flush=True)
      print("  \u2022 A mixture of the two forms fails this check even when every file exists", flush=True)
      print("  \u2022 These are the Parse pipeline's outputs; check that upstream run completed and wrote them here", flush=True)
      print(f"  \u2022 List what is actually there: `ls -l {sc_outdir}`", flush=True)
      sys.exit(1)
   
   #======================= Define Output Files =======================
   matrix_transposed_file = os.path.join(sc_outdir, "matrix.mtx.gz")
   barcodes_translated_file = os.path.join(sc_outdir, "barcodes.tsv.gz")
   features_updated_file = os.path.join(sc_outdir, "features.tsv.gz")

   transpose_matrix(mtx_file, matrix_transposed_file)
   translate_cellbarcodes(barcodes_file, barcodes_translated_file, parseCBDict)
   generate_features(features_file, features_updated_file)

   print("Single Cell Output Conversion Done")

if __name__ == "__main__":
    main()

