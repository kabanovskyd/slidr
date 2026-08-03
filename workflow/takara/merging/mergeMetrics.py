import pandas as pd
import sys
import os

def safe_pct(numerator, denominator, default=0.0):
    """Percentage (numerator/denominator*100), reporting `default` instead of inf/nan when denominator is zero."""
    if denominator == 0:
        print(f"WARNING: division by zero avoided (numerator={numerator}); reporting {default}")
        return default
    return 100 * numerator / denominator


def load_and_merge(input_dir, file_suffix='summary_metrics.csv'):
    """
    This function loads CSV files from the specified directory, concatenates them by columns, and returns the merged DataFrame.
    
    Parameters:
    input_dir (str): The directory containing the CSV files.
    file_suffix (str): The suffix of the CSV files to load (default is 'summary_metrics.csv').
    
    Returns:
    DataFrame: The merged DataFrame with Sample_IDs as column names.
    """

    # List all CSV files that end with the specified suffix
    csv_files = [file for file in os.listdir(input_dir) if file.endswith(file_suffix)]

    dfs = []

    # Load each CSV file into a dataframe and store it in the list
    for file in csv_files:
        file_path = os.path.join(input_dir, file)
        print(f"Loading file: {file_path}")
        df = pd.read_csv(file_path, index_col="Metrics")
        dfs.append(df)

    if not dfs:
        raise ValueError(
            f"No files matching suffix '{file_suffix}' found in {input_dir}\n"
            "Troubleshooting:\n"
            " \u2022 This step merges the per-partition metrics files produced by the Flex spatial-profiling step\n"
            " \u2022 An empty directory means none of this sample's partitions completed profiling\n"
            " \u2022 Check the per-partition logs in the run's takara_pipeline.log for the underlying error\n"
            f" \u2022 List what is actually there: `ls -l {input_dir}`"
        )

    # pd.concat(axis=1) silently fills any row mismatch between inputs with NaN, so
    # flag any file whose set of metrics differs from the first before merging
    reference_index = dfs[0].index
    for file, df in zip(csv_files, dfs):
        if not df.index.equals(reference_index):
            missing = set(reference_index) - set(df.index)
            extra = set(df.index) - set(reference_index)
            print(f"WARNING: {file} has a different set of metrics than {csv_files[0]} "
                  f"(missing: {sorted(missing)}, extra: {sorted(extra)}); merged values may contain NaN")

    # Merge the dataframes by concatenating columns
    merged_df = pd.concat(dfs, axis=1)

    # Set column names based on "Sample_ID" row and remove the "Sample_ID" row
    merged_df.columns = merged_df.loc["Sample_ID"]
    merged_df = merged_df.drop("Sample_ID")
    
    return merged_df


def process_df(merged_df):
    """
    This function selects specific rows from the input DataFrame, converts the values to float, 
    sums them across the columns, and returns a DataFrame with the summed values in a 'Value' column.
    
    Parameters:
    merged_df (DataFrame): The input DataFrame.
    
    Returns:
    DataFrame: A processed DataFrame containing only the 'Value' column with summed values.
    """
    
    selected_rows = [
        "Total_nuclei_from_single-nuclei_sequencing_library",
        "Nuclei_from_single-nuclei_sequencing_library_found_in_Trekker_library",
        "Nuclei_from_single-nuclei_sequencing_library_found_in_Trekker_library_with_valid_spatial_barcodes",
        #"Total_nuclei_positioned",
        #"Total_nuclei_positioned_with_1_spatial_location",
        "Total_readpairs_in_Trekker_library",
        "Readpairs_with_proper_structure",
        "Readpairs_used_for_matching_to_single_nuclei_barcodes",
        "Readpairs_matched_single_nuclei_barcodes",
        "Readpairs_matched_to_single_nuclei_barcodes_with_valid_spatial_barcodes"
    ]
    
    # a selected metric absent from the merged frame (i.e. present in no input file) would make
    # .loc raise an opaque KeyError; surface it clearly instead
    missing_rows = [r for r in selected_rows if r not in merged_df.index]
    if missing_rows:
        print(f"[ERROR]: the following metrics are absent from all input summary files: {missing_rows}", flush=True)
        print("Troubleshooting:", flush=True)
        print(" \u2022 These metrics are written per partition by genmetrics.py, so a missing one means that step wrote an unexpected layout", flush=True)
        print(f" \u2022 Metrics that were found: {sorted(merged_df.index)}", flush=True)
        print(" \u2022 Files left over from an older version of the Trekker scripts are the usual cause; delete them and re-run --spatial-analysis --force", flush=True)
        sys.exit(1)

    merged_df = merged_df.loc[selected_rows]
    merged_df = merged_df.astype(float)

    # sum(axis=1) defaults to skipna=True, which would silently under-count the total
    # for any metric missing from one or more samples instead of surfacing the gap
    nan_counts = merged_df.isna().sum(axis=1)
    rows_with_missing = nan_counts[nan_counts > 0]
    if not rows_with_missing.empty:
        print(f"WARNING: the following metrics are missing for one or more samples; "
              f"their summed totals below will under-count the true total:\n{rows_with_missing}")

    merged_df['Value'] = merged_df.sum(axis=1)
    merged_df = merged_df[['Value']]
    
    return merged_df


def calculate_qc_metrics(merged_df):
    """
    This function calculates QC metrics for nuclei and reads in the input DataFrame,
    and updates the DataFrame with the calculated percentage values.

    Parameters:
    merged_df (DataFrame): The input DataFrame with the necessary metrics.

    Returns:
    DataFrame: The updated DataFrame with calculated QC metrics.
    """

    # QC Nuclei
    merged_df.loc["Pct_nuclei_in_Trekker_library", 'Value'] = safe_pct(
        merged_df.loc["Nuclei_from_single-nuclei_sequencing_library_found_in_Trekker_library", 'Value'],
        merged_df.loc["Total_nuclei_from_single-nuclei_sequencing_library", 'Value']
    )

    merged_df.loc["Pct_nuclei_in_Trekker_library_with_valid_spatial_barcodes", 'Value'] = safe_pct(
        merged_df.loc["Nuclei_from_single-nuclei_sequencing_library_found_in_Trekker_library_with_valid_spatial_barcodes", 'Value'],
        merged_df.loc["Total_nuclei_from_single-nuclei_sequencing_library", 'Value']
    )

    """
    merged_df.loc["Pct_nuclei_positioned", 'Value'] = (
        100 * merged_df.loc["Total_nuclei_positioned", 'Value'] /
        merged_df.loc["Total_nuclei_from_single-nuclei_sequencing_library", 'Value']
    )

    merged_df.loc["Pct_nuclei_positioned_with_1_spatial_location", 'Value'] = (
        100 * merged_df.loc["Total_nuclei_positioned_with_1_spatial_location", 'Value'] /
        merged_df.loc["Total_nuclei_from_single-nuclei_sequencing_library", 'Value']
    )

    merged_df.loc["Pct_nuclei_positioned_with_2+_spatial_locations", 'Value'] = (
        100 * (merged_df.loc["Total_nuclei_positioned", 'Value'] -
        merged_df.loc["Total_nuclei_positioned_with_1_spatial_location", 'Value']) /
        merged_df.loc["Total_nuclei_from_single-nuclei_sequencing_library", 'Value']
    )
    """

    # QC Reads
    merged_df.loc["Pct_readpairs_with_proper_structure", 'Value'] = safe_pct(
        merged_df.loc["Readpairs_with_proper_structure", 'Value'],
        merged_df.loc["Total_readpairs_in_Trekker_library", 'Value']
    )

    merged_df.loc["Pct_readpairs_matched_to_single_nuclei_barcodes", 'Value'] = safe_pct(
        merged_df.loc["Readpairs_matched_single_nuclei_barcodes", 'Value'],
        merged_df.loc["Total_readpairs_in_Trekker_library", 'Value']
    )

    merged_df.loc["Pct_readpairs_matched_to_single_nuclei_barcodes_with_valid_spatial_barcodes", 'Value'] = safe_pct(
        merged_df.loc["Readpairs_matched_to_single_nuclei_barcodes_with_valid_spatial_barcodes", 'Value'],
        merged_df.loc["Readpairs_matched_single_nuclei_barcodes", 'Value']
    )

    merged_df.loc["Pct_useful_reads", 'Value'] = safe_pct(
        merged_df.loc["Readpairs_matched_to_single_nuclei_barcodes_with_valid_spatial_barcodes", 'Value'],
        merged_df.loc["Total_readpairs_in_Trekker_library", 'Value']
    )

    return merged_df


def format_dataframe(df):
    """
    This function formats the 'Value' column in the DataFrame.
    If a row contains 'Pct' in its name, values are formatted to two decimal places.
    Otherwise, values are formatted with zero decimal places. The result is saved to a CSV file.
    
    Parameters:
    df (DataFrame): The input DataFrame.

    Returns:
    DataFrame: The reformatted DataFrame.
    """ 
    export_df = df.copy()

    # Convert all columns to strings to avoid errors with non-numeric types
    export_df = export_df.apply(lambda col: col.map(str))

    # Iterate through rows and format based on row names
    for row in export_df.index:
        value = export_df.loc[row, 'Value']
        if "Pct" in row:
            # Format to two decimal places
            export_df.loc[row, 'Value'] = '%.2f' % float(value) if isinstance(value, (int, float, str)) else value
        else:
            # Format to zero decimal places
            export_df.loc[row, 'Value'] = '%.0f' % float(value) if isinstance(value, (int, float, str)) else value
    
    return export_df

                                 
def main():
    """
    Main function to execute the pipeline for loading, processing, and exporting CSV data.
    """
    
    if len(sys.argv) < 4:
        print(f"[ERROR]: expected 3 arguments (input dir, output prefix, output dir), got {len(sys.argv) - 1}", flush=True)
        print("Troubleshooting:", flush=True)
        print(" \u2022 This script is normally invoked by the Flex merge stage, not by hand", flush=True)
        print(" \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`", flush=True)
        print(f" \u2022 Received: {sys.argv[1:]}", flush=True)
        sys.exit(1)
    
    # Get input arguments
    input_dir = sys.argv[1]
    output_prefix = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Load and merge data
    merged_df = load_and_merge(input_dir)
    print("Merged DataFrame:\n", merged_df.shape)
    print(merged_df.head())
    
    # Process DataFrame
    merged_df = process_df(merged_df)
    print("Processed DataFrame:\n", merged_df.shape)
    
    # Calculate QC metrics
    merged_df = calculate_qc_metrics(merged_df)
    print("DataFrame with QC metrics:\n", merged_df.shape)
    
    # Format the DataFrame and save it
    formatted_df = format_dataframe(merged_df)
    formatted_df.to_csv(os.path.join(output_dir, output_prefix + "_summary_metrics_merged.csv"))
    print(f"Formatted DataFrame saved to {output_dir}")

if __name__ == "__main__":
    main()                                 
                                 
                                 
                                 
                                 
