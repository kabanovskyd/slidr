#!/usr/bin/env python

import pandas as pd
import sys
import os
import numpy as np

def log(msg):
    print(msg)

def extract_info(summary, key_word, next_row, info_index):
    keyword_row_index = summary.index[summary[0].str.contains(key_word)]
    if len(keyword_row_index) == 0:
        raise ValueError(
            f"extract_info: could not find line containing '{key_word}' in bead_matching.log\n"
            "Troubleshooting:\n"
            " • This metric is scraped out of the bead-matching step's log, so a missing line means that step did not run to completion\n"
            " • Check the tail of bead_matching.log in the run's log directory for its own error\n"
            " • Re-run the Flex spatial analysis with `--spatial-analysis --force` to regenerate it"
        )
    info_row_index = keyword_row_index[0]+next_row
    info_row = summary.iloc[info_row_index][0]
    info_split = info_row.strip().split(" ")
    info = int(info_split[info_index])
    return info

def create_df(name_lst, value_list):
    df = pd.DataFrame(list(zip(name_lst, value_list)),
               columns =['Metrics', 'Value'])
    return df

def safe_div(numerator, denominator, default=0.0):
    """Division that reports `default` instead of inf/nan when denominator is zero."""
    if denominator == 0:
        log(f"WARNING: division by zero avoided (numerator={numerator}); reporting {default}")
        return default
    return numerator / denominator

def pct(numerator, denominator, default=0.0):
    """Percentage (numerator/denominator*100), reporting `default` instead of inf/nan when denominator is zero."""
    return 100 * safe_div(numerator, denominator, default)

print(sys.argv)
print(len(sys.argv))

if len(sys.argv) != 7:
    # sys.exit(1), not sys.exit(): a bare exit() returns 0 and the calling stage would treat this
    # failure as a success and carry on with no metrics written
    print(f"[ERROR]: expected 6 arguments (sample id, tile id, single cell platform, log dir, "
          f"misc dir, output dir), got {len(sys.argv) - 1}", flush=True)
    print("Troubleshooting:", flush=True)
    print(" • This script is normally invoked by the Flex spatial-profiling stage, not by hand", flush=True)
    print(" • Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`", flush=True)
    print(f" • Received: {sys.argv[1:]}", flush=True)
    sys.exit(1)
else:
    sample_id = sys.argv[1]
    tile_id = sys.argv[2] 
    scRNAseq_assay = sys.argv[3]
    logdir = sys.argv[4]
    trekkerout_misc = sys.argv[5]
    trekkerout_main = sys.argv[6]

path_nuclei_bc_match_summary = os.path.join(trekkerout_misc, "matcher_summary_"+sample_id+".txt") 
path_nuclei_zero_rep_summary = os.path.join(trekkerout_misc, sample_id+"_rep_summary.csv")
path_spatial_bc_match_summary = os.path.join(logdir, "bead_matching.log")
path_spatial_positioning_summary = os.path.join(trekkerout_main, "plots", "summary_metrics.csv") 
path_spatial_bc_read_match_summary = os.path.join(trekkerout_misc, sample_id+"_properreads_matched_to_spatial_whitelist.csv")
path_position_confidence_summary = os.path.join(trekkerout_misc, sample_id+"_summary_position_conf_consolidated.txt")
path_bead_removal_summary = os.path.join(trekkerout_misc, sample_id+"_beadRemoval_summary.csv")

spatial_positioning_summary = pd.read_csv(path_spatial_positioning_summary, sep=",", index_col=0)
spatial_nuclei_zero_rep_summary = pd.read_csv(path_nuclei_zero_rep_summary, sep=",")
spatial_bc_match_summary = pd.read_csv(path_spatial_bc_match_summary, sep="\t", header=None)
nuclei_bc_match_summary = pd.read_csv(path_nuclei_bc_match_summary, sep=",")
spatial_bc_read_match_summary = pd.read_csv(path_spatial_bc_read_match_summary, sep=",")
position_confidence_summary = pd.read_csv(path_position_confidence_summary, sep="\t")
bead_removal_summary = pd.read_csv(path_bead_removal_summary, index_col=0, header=0)

# every metric below is read from row 0 of these single-row summaries via .iloc[0]; fail with a
# clear message if any of them is header-only (no data rows) instead of an opaque IndexError
for _label, _df in [
    ("spatial positioning summary (summary_metrics.csv)", spatial_positioning_summary),
    ("nuclei barcode match summary", nuclei_bc_match_summary),
    ("spatial barcode read-match summary", spatial_bc_read_match_summary),
    ("position confidence summary", position_confidence_summary),
]:
    if len(_df) == 0:
        raise ValueError(
            f"genmetrics: {_label} has no data rows; cannot compute metrics for {sample_id}\n"
            "Troubleshooting:\n"
            " • A header-only summary means the step that produces it found nothing to report\n"
            f" • This usually means no nuclei were positioned for {sample_id}; check the spatial-profiling log for that sample\n"
            " • Check the `Flex Probe Barcode IDs` metadata column matches the barcodes actually used, so reads are not all discarded\n"
            " • Re-run the Flex spatial analysis with `--spatial-analysis --force` once the cause is fixed"
        )

#df_PARAMETERS
log("df_PARAMETERS")
name_lst = [
        "eps",
        "Min_spatial_barcodes_used_to_locate_a_nucleus_centroid",
        "Maximum_UMI_cutoff"
        ]

value_list = [
        '%.0f'%(spatial_positioning_summary.iloc[0]["eps"]),
        '%.0f'%(spatial_positioning_summary.iloc[0]["minPts"]),
        '%.0f'%(spatial_positioning_summary.iloc[0]["nUMI_cutoff_combined"])
        ]

df_PARAMETERS = create_df(name_lst, value_list)
print(df_PARAMETERS.shape)

#df_METADATA
log("df_METADATA")
name_lst = [
        "Sample_ID",
        "Single_cell_assay",
        "Tile_ID"
        ]

value_list=[
        sample_id,
        scRNAseq_assay,
        tile_id
        ]

df_METADATA = create_df(name_lst, value_list)
print(df_METADATA.shape)

#df_NUCLEI_RECOVERY_QC
log("df_NUCLEI_RECOVERY_QC")
Total_nuclei_from_WTA_library = nuclei_bc_match_summary.iloc[0]["cb_snRNAseq"]
Total_nuclei_from_WTA_library_found_in_tags_library = nuclei_bc_match_summary.iloc[0]["cb_matched"]
Pct_nuclei_from_WTA_library_found_in_tags_library = pct(Total_nuclei_from_WTA_library_found_in_tags_library, Total_nuclei_from_WTA_library)

name_lst = [
        "Total_nuclei_from_single-nuclei_sequencing_library",
        "Nuclei_from_single-nuclei_sequencing_library_found_in_Trekker_library",
        "Pct_nuclei_in_Trekker_library",
        "Nuclei_from_single-nuclei_sequencing_library_found_in_Trekker_library_with_valid_spatial_barcodes",
        "Pct_nuclei_in_Trekker_library_with_valid_spatial_barcodes",
        "Nuclei_positioned",
        "Pct_nuclei_positioned",
        "Nuclei_confidently_positioned",
        "Pct_nuclei_confidently_positioned"
        #"Pct_nuclei_positioned_with_2+_spatial_locations"
        ]

value_list = [
        '%.0f'%(Total_nuclei_from_WTA_library),
        '%.0f'%(Total_nuclei_from_WTA_library_found_in_tags_library),
        '%.2f'%(Pct_nuclei_from_WTA_library_found_in_tags_library), 
        '%.0f'%(spatial_positioning_summary.iloc[0]["total_cells"]),
        '%.2f'%(pct(spatial_positioning_summary.iloc[0]["total_cells"], Total_nuclei_from_WTA_library)),
        '%.0f'%(spatial_positioning_summary.iloc[0]['spatially_mapped_any_cluster_number']),
        '%.2f'%(pct(spatial_positioning_summary.iloc[0]['spatially_mapped_any_cluster_number'], Total_nuclei_from_WTA_library)),
        '%.0f'%(spatial_positioning_summary.iloc[0]['spatially_mapped_top_cluster']),
        '%.2f'%(pct(spatial_positioning_summary.iloc[0]['spatially_mapped_top_cluster'], Total_nuclei_from_WTA_library))
        #'%.2f'%(100*((spatial_positioning_summary.iloc[0]['spatially_mapped_any_cluster_number']-spatial_positioning_summary.iloc[0]['spatially_mapped_top_cluster'])/Total_nuclei_from_WTA_library))
        ]

df_NUCLEI_RECOVERY_QC = create_df(name_lst, value_list)
print(df_NUCLEI_RECOVERY_QC.shape)

#df_READ_RECOVERY_QC_PART1
log("df_READ_RECOVERY_QC_PART1")

Total_readpairs_in_tags_library = nuclei_bc_match_summary.iloc[0]["Total_reads"]
Readpairs_with_exact_UP = nuclei_bc_match_summary.iloc[0]["UP_site_matches"]
Pct_readpairs_with_exact_UP = pct(Readpairs_with_exact_UP, Total_readpairs_in_tags_library)

name_lst = [
        "Total_readpairs_in_Trekker_library",
        "Readpairs_with_proper_structure",
        "Pct_readpairs_with_proper_structure"
        ]

value_list = [
        '%.0f'%(Total_readpairs_in_tags_library),
        '%.0f'%(Readpairs_with_exact_UP),
        '%.2f'%(Pct_readpairs_with_exact_UP) 
        ]

df_READ_RECOVERY_QC_PART1 = create_df(name_lst, value_list)
print(df_READ_RECOVERY_QC_PART1.shape)

#df_READ_RECOVERY_QC_PART2
log("df_READ_RECOVERY_QC_PART2")
nuclei_bc_match_summary = pd.read_csv(path_nuclei_bc_match_summary, sep=",")
Readpairs_matched_to_WTA_nuclei_barcodes = nuclei_bc_match_summary.iloc[0]["cb_matched_reads"]
Fraction_properstructure_readpairs_matched_to_WTA_nuclei_barcodes = safe_div(Readpairs_matched_to_WTA_nuclei_barcodes, nuclei_bc_match_summary.iloc[0]["UP_site_matches"])
Pct_readpairs_matched_to_WTA_nuclei_barcodes = Pct_readpairs_with_exact_UP*Fraction_properstructure_readpairs_matched_to_WTA_nuclei_barcodes
Pct_readpairs_matched_to_WTA_nuclei_barcodes_and_spatial_whitelist = pct(spatial_bc_read_match_summary.iloc[0]['proper_reads_matched_to_spatial_whitelist'], Readpairs_matched_to_WTA_nuclei_barcodes)
Pct_useful_reads = Pct_readpairs_matched_to_WTA_nuclei_barcodes*Pct_readpairs_matched_to_WTA_nuclei_barcodes_and_spatial_whitelist/100

name_lst = [
        "Readpairs_used_for_matching_to_single_nuclei_barcodes",
        "Readpairs_matched_single_nuclei_barcodes",
        "Pct_readpairs_matched_to_single_nuclei_barcodes",
        "Readpairs_matched_to_single_nuclei_barcodes_with_valid_spatial_barcodes",
        "Pct_readpairs_matched_to_single_nuclei_barcodes_with_valid_spatial_barcodes",
        "Pct_useful_reads"
        ]

value_list = [
        '%.0f'%(nuclei_bc_match_summary.iloc[0]["UP_site_matches"]),
        '%.0f'%(Readpairs_matched_to_WTA_nuclei_barcodes),
        '%.2f'%(Pct_readpairs_matched_to_WTA_nuclei_barcodes),
        '%.0f'%(spatial_bc_read_match_summary.iloc[0]['proper_reads_matched_to_spatial_whitelist']),
        '%.2f'%(Pct_readpairs_matched_to_WTA_nuclei_barcodes_and_spatial_whitelist),
        '%.2f'%(Pct_useful_reads)
        ]

df_READ_RECOVERY_QC_PART2 = create_df(name_lst, value_list)
print(df_READ_RECOVERY_QC_PART2.shape)

#df_BARCODE_MATCHING_QC
log("df_BARCODE_MATCHING_QC")

Total_seqenced_spatial_barcodes_in_tags_data = extract_info(spatial_bc_match_summary, "Number of barcode from Illumina reads", 0, -1)
Total_seqenced_spatial_barcodes_perfectly_matched_to_whitelist = extract_info(spatial_bc_match_summary, "exact match", 0, -1)
Total_seqenced_spatial_barcodes_approximately_matched_to_whitelist = extract_info(spatial_bc_match_summary, "fuzzy match", 0, -1)

name_lst = [
        "Total_sequenced_spatial_barcodes", 
        "Sequenced_spatial_barcodes_perfectly_matched_to_whitelist", 
        "Sequenced_spatial_barcodes_approximately_matched_to_whitelist",
        "Pct_valid_spatial_barcodes"
        ]

value_list = [
        '%.0f'%(Total_seqenced_spatial_barcodes_in_tags_data),
        '%.0f'%(Total_seqenced_spatial_barcodes_perfectly_matched_to_whitelist),
        '%.0f'%(Total_seqenced_spatial_barcodes_approximately_matched_to_whitelist),
        '%.2f'%(pct(Total_seqenced_spatial_barcodes_perfectly_matched_to_whitelist+Total_seqenced_spatial_barcodes_approximately_matched_to_whitelist, Total_seqenced_spatial_barcodes_in_tags_data))
        ]

df_BARCODE_MATCHING_QC = create_df(name_lst, value_list)
print(df_BARCODE_MATCHING_QC.shape)

## df_DEPTH
log("df_DEPTH")

name_lst = [
        "Median_spatial_barcodes_per_nuclei",
        "Mean_spatial_barcodes_per_nuclei",
        "Median_spatial_barcodes_UMI_per_nuclei",
        "Mean_spatial_barcodes_UMI_per_nuclei",
        "Median_reads_per_nuclei",
        "Mean_reads_per_nuclei",
        "Median_useful_reads_per_nuclei",
        "Mean_useful_reads_per_nuclei"
        ]
value_list = [
        '%.0f'%(spatial_positioning_summary.iloc[0]["median_unique_SB_per_cell"]),
        '%.0f'%(spatial_positioning_summary.iloc[0]["mean_unique_SB_per_cell"]),
        '%.0f'%(spatial_positioning_summary.iloc[0]["median_SB_UMI_per_cell"]),
        '%.0f'%(spatial_positioning_summary.iloc[0]["mean_SB_UMI_per_cell"]),
        '%.0f'%(nuclei_bc_match_summary.iloc[0]["Median_reads_per_nuclei"]),
        '%.0f'%(nuclei_bc_match_summary.iloc[0]["Mean_reads_per_nuclei"]),
        '%.0f'%(spatial_bc_read_match_summary.iloc[0]['Median_useful_reads_per_nuclei']),
        '%.0f'%(spatial_bc_read_match_summary.iloc[0]['Mean_useful_reads_per_nuclei'])
        ]

df_DEPTH = create_df(name_lst, value_list)
print(df_DEPTH.shape)

## df_SIGNAL_VS_NOISE
log("df_SIGNAL_VS_NOISE")
name_lst = [
        "median_proportion_unique_SB_top_cluster_total",
        "mean_proportion_unique_SB_top_cluster_total",
        "median_proportion_SB_UMI_top_cluster_total",
        "mean_proportion_SB_UMI_top_cluster_total",
        "median_mapped_cells_proportion_unique_SB_top_cluster",
        "mean_mapped_cells_proportion_unique_SB_top_cluster",
        "median_mapped_cells_proportion_SB_UMI_top_cluster",
        "mean_mapped_cells_proportion_SB_UMI_top_cluster"
        ]

# name_lst here is a hardcoded set of expected stat columns; a reduced summary_metrics.csv that is
# non-empty but missing one of them would otherwise raise an opaque KeyError on the .iloc[0] subset
_missing_cols = [c for c in name_lst if c not in spatial_positioning_summary.columns]
if _missing_cols:
    raise ValueError(
        f"genmetrics: summary_metrics.csv is missing expected column(s): {_missing_cols}\n"
        "Troubleshooting:\n"
        f" • File: {path_spatial_positioning_summary}\n"
        " • These columns are written by the spatial-positioning step, so a partial file means that step did not finish\n"
        " • Check the spatial-profiling log for that step's own error\n"
        " • A file left over from an older version of the Trekker scripts is the other common cause; delete it and re-run --spatial-analysis --force"
    )
value_list = ['%.4f'% m for m in spatial_positioning_summary.iloc[0][name_lst]]
name_lst = [m.capitalize() for m in name_lst]
df_SIGNAL_VS_NOISE = create_df(name_lst, value_list)
print(df_SIGNAL_VS_NOISE.shape)

## df_position_confidence_summary
log("df_position_confidence_summary")

name_lst = list(position_confidence_summary.columns)
value_list = ['%.0f'% m for m in position_confidence_summary.iloc[0][name_lst]]
name_lst = ["Nuclei_" + s for s in list(position_confidence_summary.columns)]
df_position_confidence_summary = create_df(name_lst, value_list)
print(df_position_confidence_summary.shape)

## df_pct_position_confidence_summary
log("df_pct_position_confidence_summary")

position_confidence_total = int(position_confidence_summary.sum(axis=1).iloc[0])
if position_confidence_total == 0:
    log("WARNING: position_confidence_summary row sums to zero; reporting 0.0 for all pct_position_confidence_summary values")
    pct_position_confidence_summary = position_confidence_summary * 0.0
else:
    pct_position_confidence_summary = position_confidence_summary.div(position_confidence_total)*100
name_lst = ["pct_" + s for s in list(pct_position_confidence_summary.columns)]
value_list = ['%.2f'% m for m in list(pct_position_confidence_summary.iloc[0])]
df_pct_position_confidence_summary = create_df(name_lst, value_list)
print(df_pct_position_confidence_summary.shape)

## bead_removal_summary (% salvaged)
bead_removal_summary = pd.read_csv(path_bead_removal_summary, index_col=0, header=0)

## df_BEAD_REMOVAL_SUMMARY
log("df_bead_removal_summary")

# catch the col name for the ">=4" column
if ">=4" in bead_removal_summary.columns:
    col_4plus = ">=4"
else:
    # Step 2: fallback to any column name that contains "4"
    candidates = [c for c in bead_removal_summary.columns if "4" in str(c)]
    if candidates:
        col_4plus = candidates[0]  # take the first match
    else:
        raise ValueError(
            "no column for the '>=4 clusters' bead-removal bin (nor any variant containing '4') in "
            f"{path_bead_removal_summary}\n"
            f" Columns present: {list(bead_removal_summary.columns)}\n"
            "Troubleshooting:\n"
            " • This file is written by the bead-removal step, so a missing bin means that step wrote an unexpected layout\n"
            " • A file left over from an older version of the Trekker scripts is the usual cause; delete it and re-run --spatial-analysis --force\n"
            " • If the layout is genuinely new, this scraper needs updating -- please report it at https://github.com/kabanovskyd/slidr/issues"
        )

#print(f"Using column: '{col_4plus}' as the 4plus column name")

def bead_count(row_label, col_label):
    """
    Nuclei count from bead_removal_summary, defaulting to 0 when the cluster-count column is
    absent -- a small sample may have no nuclei in a given spatial-location bin, in which case the
    column is missing rather than present-with-zero, and a direct .loc would raise KeyError.
    """
    if col_label not in bead_removal_summary.columns:
        log(f"WARNING: bead_removal_summary has no '{col_label}' column; using 0")
        return 0
    return bead_removal_summary.loc[row_label, col_label]

name_lst = [
        "Nuclei_o_1",
        "Nuclei_o_2",
        "Nuclei_o_3",
        "Nuclei_o_>=4",
        "pct_o_1",
        "pct_o_2",
        "pct_o_3",
        "pct_o_4",
        "Nuclei_salvaged_2",
        "Nuclei_salvaged_3",
        "Nuclei_salvaged_>=4",
        "Pct_nuclei_salvaged_from_2_spatial_locations",
        "Pct_nuclei_salvaged_from_3_spatial_locations",
        "Pct_nuclei_salvaged_from_4+_spatial_locations",
        "Pct_nuclei_salvaged"
            ]

Total_nuclei_from_WTA_library = nuclei_bc_match_summary.iloc[0]["cb_snRNAseq"]
Nuclei_o_1 = spatial_positioning_summary.iloc[0]["spatially_mapped_top_cluster_o"]
Nuclei_salvaged = bead_removal_summary.loc["Nuclei_with_beads"].sum()

value_list = [
        '%.0f'%(Nuclei_o_1),
        '%.0f'%(bead_count("Nuclei","2")),
        '%.0f'%(bead_count("Nuclei","3")),
        '%.0f'%(bead_removal_summary.loc["Nuclei",col_4plus]),
        '%.2f'%(pct(Nuclei_o_1, Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(bead_count("Nuclei","2"), Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(bead_count("Nuclei","3"), Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(bead_removal_summary.loc["Nuclei",col_4plus], Total_nuclei_from_WTA_library)),
        '%.0f'%(bead_count("Nuclei_with_beads","2")),
        '%.0f'%(bead_count("Nuclei_with_beads","3")),
        '%.0f'%(bead_removal_summary.loc["Nuclei_with_beads",col_4plus]),
        '%.2f'%(pct(bead_count("Nuclei_with_beads","2"), Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(bead_count("Nuclei_with_beads","3"), Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(bead_removal_summary.loc["Nuclei_with_beads",col_4plus], Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(Nuclei_salvaged, Total_nuclei_from_WTA_library))
        ]

df_bead_removal_summary = create_df(name_lst, value_list)
print(df_bead_removal_summary.shape)


#df_ZERO_REPARAMETER_SUMMARY
log("df_ZERO_REPARAMETER_SUMMARY")
if 1 in spatial_nuclei_zero_rep_summary.index:
    Nuclei_1_rep = spatial_nuclei_zero_rep_summary.loc[1, "2"]
else:
    log("WARNING: spatial_nuclei_zero_rep_summary has no row for rep index 1; reporting 0")
    Nuclei_1_rep = 0
Nuclei_2_rep = spatial_nuclei_zero_rep_summary.loc[2:, "2"].sum()

name_lst = [
        "Pct_nuclei_1_rep",
        "Pct_nuclei_2plus_rep",
        "nuclei_1_rep",
        "nuclei_2plus_rep"
        ]

value_list = [
        '%.2f'%(pct(Nuclei_1_rep, Total_nuclei_from_WTA_library)),
        '%.2f'%(pct(Nuclei_2_rep, Total_nuclei_from_WTA_library)),
        '%.0f'%(Nuclei_1_rep),
        '%.0f'%(Nuclei_2_rep)
        ]

df_zero_reparameter_summary = create_df(name_lst, value_list)
print(df_zero_reparameter_summary.shape)

#Concatenation
log("Concatenating all metrics")
df_vertical_concat = pd.concat([df_METADATA,
                                df_PARAMETERS, 
                                df_NUCLEI_RECOVERY_QC,
                                df_READ_RECOVERY_QC_PART1,
                                df_READ_RECOVERY_QC_PART2,
                                df_BARCODE_MATCHING_QC,
                                df_DEPTH,
                                df_SIGNAL_VS_NOISE,
                                df_position_confidence_summary,
                                df_pct_position_confidence_summary, 
                                df_bead_removal_summary,
                                df_zero_reparameter_summary],axis=0)

print(df_vertical_concat.shape)
df_vertical_concat.set_index('Metrics', inplace=True)

log("Output compiled metrics")
df_vertical_concat.to_csv(os.path.join(trekkerout_main, sample_id+"_summary_metrics.csv"), index=True, header=True)

