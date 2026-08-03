import pandas as pd
import os
import sys

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

# Every failure in this file is a malformed barcode reaching the mismatch analysis, which points at
# the puck whitelist or the matching step upstream rather than at anything the user typed. Shared so
# the two guards below give the same advice.
BARCODE_HINTS = (
    "Troubleshooting:\n"
    " • Barcodes reaching this point come from the puck whitelist and the bead-matching step\n"
    " • Check the `Puck ID` metadata column names the puck this sample was actually run on\n"
    " • Check the puck's BeadBarcodes file holds uniform-length 14bp barcodes\n"
    " • Delete the cached puck under the run's flex/pucks/ directory to force a fresh download, then re-run --spatial-analysis --force"
)


def find_mismatch(chaine1, chaine2):
    if len(chaine1) != len(chaine2):
        raise ValueError(
            f"find_mismatch: barcode length mismatch ({len(chaine1)} vs {len(chaine2)}): "
            f"{chaine1} vs {chaine2}\n{BARCODE_HINTS}"
        )
    return [i for i in range(len(chaine1)) if chaine1[i] != chaine2[i]]

def increment_mismatch_reads(chaine1, chaine2, nreads):
    BARCODE_LENGTH = 14
    if len(chaine1) != BARCODE_LENGTH or len(chaine2) != BARCODE_LENGTH:
        raise ValueError(
            f"increment_mismatch_reads: expected {BARCODE_LENGTH}bp barcodes, got lengths "
            f"{len(chaine1)} and {len(chaine2)}: {chaine1} vs {chaine2}\n{BARCODE_HINTS}"
        )
    list_mismatches = [0]*BARCODE_LENGTH
    pos = find_mismatch(chaine1, chaine2)
    target = [nreads]*len(pos)
    for x,y in zip(pos, target):
        list_mismatches[x]=y
    return list_mismatches

print(sys.argv)
print(len(sys.argv))

if len(sys.argv) != 3:
    # sys.exit(1), not sys.exit(): a bare exit() returns 0 and the calling stage would treat this
    # failure as a success and carry on with no mismatch analysis written
    print(f"[ERROR]: expected 2 arguments (sample id, misc directory), got {len(sys.argv) - 1}", flush=True)
    print("Troubleshooting:", flush=True)
    print(" • This script is normally invoked by the Flex spatial-profiling stage, not by hand", flush=True)
    print(" • Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`", flush=True)
    print(f" • Received: {sys.argv[1:]}", flush=True)
    sys.exit(1)
else:
    sample_id = sys.argv[1]
    dir_misc = sys.argv[2]

bcmatching_results = pd.read_csv(os.path.join(dir_misc,"matching_result_"+sample_id+".csv"), index_col=0)
sequenced_barcodes = pd.read_csv(os.path.join(dir_misc,"reads_per_SB_"+sample_id+".txt"), index_col=0)
sequenced_barcodes.columns = ["num_reads"]


bcmatching_results = bcmatching_results.join(sequenced_barcodes, how="left")
missing_reads = bcmatching_results["num_reads"].isna().sum()
if missing_reads > 0:
    print(f"WARNING: {missing_reads} matched barcodes had no corresponding entry in reads_per_SB file; treating as 0 reads")
    bcmatching_results["num_reads"] = bcmatching_results["num_reads"].fillna(0)
bcmatching_results["illumina_barcodes"] = bcmatching_results.index

total_reads = bcmatching_results["num_reads"].sum()
print("Total useful reads based on SB")
print(total_reads)

# calculate the number of reads that have a mismatch at each spatial barcode location
bcmatching_results["mismatches"] = bcmatching_results.apply(lambda x: increment_mismatch_reads(x.illumina_barcodes, x.matched_beadbarcode, x.num_reads), axis=1)
pd_mismatches = pd.DataFrame(bcmatching_results["mismatches"].tolist())
print(pd_mismatches)

if total_reads == 0:
    print("WARNING: total_reads is 0; reporting 0% mismatch frequency instead of dividing by zero")
    list_mismatches = pd_mismatches.sum() * 0.0
else:
    list_mismatches = round(100*pd_mismatches.sum()/total_reads,2)

list_mismatches.to_csv(os.path.join(dir_misc, sample_id+"_mismatchfreq.csv"), index=None)

#useful reads per nucei
# reads_perMatchedCB_SB.txt has no header row; columns are (n_reads, cb, sb) in that order
nucleimatching_results = pd.read_csv(
    os.path.join(dir_misc, sample_id+"_reads_perMatchedCB_SB.txt"),
    sep=",", header=None, names=["n_reads", "cb", "sb"],
)
nucleimatching_results.set_index("sb", inplace=True)

nucleimatching_results = nucleimatching_results.join(bcmatching_results["matched_beadbarcode"], how="inner")
print(nucleimatching_results.head())
print("Done")

df_useful_reads_perNuclei = nucleimatching_results.groupby("cb")["n_reads"].sum()
#print(df_useful_reads_perNuclei)

print("Total useful reads based on CB")
print(df_useful_reads_perNuclei.sum())

if df_useful_reads_perNuclei.empty:
    print("WARNING: no nuclei had any useful reads matched to the spatial whitelist; reporting 0 for median/mean")
    Median_useful_reads_per_nuclei = 0
    Mean_useful_reads_per_nuclei = 0
else:
    Median_useful_reads_per_nuclei = round(df_useful_reads_perNuclei.median(),0)
    Mean_useful_reads_per_nuclei = round(df_useful_reads_perNuclei.mean(),0)

summary_reads = pd.DataFrame({'proper_reads_matched_to_spatial_whitelist': [total_reads],
                              'Median_useful_reads_per_nuclei': [Median_useful_reads_per_nuclei],
                              'Mean_useful_reads_per_nuclei': [Mean_useful_reads_per_nuclei]})
print(summary_reads)

summary_reads.to_csv(os.path.join(dir_misc, sample_id+"_properreads_matched_to_spatial_whitelist.csv"))
