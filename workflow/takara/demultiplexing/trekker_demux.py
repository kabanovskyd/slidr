#!/usr/bin/env python3

# Description:
#   Parses paried-end sequencing data for 10x Flex and demultiplexes the reads based on barcode sequences
#   Fastq.gz files are written to separate files based on the barcode sequence
# Input:
#   - Paired-end multiplexed sequencing data in fastq.gz format
#      - R1.fastq.gz
#      - R2.fastq.gz
#   - A csv file containing Sample_ID and Barcode_Label (optional) OR 'NA'
# Output:
#   - Paired demultiplexed fastq.gz files with the Sample_ID or Barcode_Label in the filename
#   - Log file with detailed progress and errors
# Usage:
#   python trekkerFX_demux.py <R1.fastq.gz> <R2.fastq.gz> <sample_file.csv>
#   python trekkerFX_demux.py <R1.fastq.gz> <R2.fastq.gz> NA

# Import modules
import re
import os
import sys
import csv
import gzip
import psutil
import concurrent.futures
from datetime import datetime

# Global variables (populated by validate_input() once argument count is confirmed)
output_dir = None
LOG = None

# Assorted functions
def format_time(td):
    """Format a timedelta into a human-readable string."""
    hours, remainder = divmod(td.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

def memory_usage():
    # Return total memory usage of script in megabytes
    process = psutil.Process()
    memory_mb = process.memory_info().rss / (1024 * 1024)
    return memory_mb

def print_valid_usage_and_exit():
    """
    Print the valid usage of the script, then exit
    """
    print('Usage: python trekkerFX_demux.py <R1.fastq.gz> <R2.fastq.gz> <sample_file.csv> <output_dir>')
    print("Usage: if you want to run the script without a sample file, enter 'NA' as the third argument")
    print('Fatal Error: Script terminated')
    if LOG is not None:
        LOG.write('Usage: python trekkerFX_demux.py <R1.fastq.gz> <R2.fastq.gz> <sample_file.csv> <output_dir>\n')
        LOG.write("Usage: if you want to run the script without a sample file, enter 'NA' as the third argument\n")
        LOG.write('Fatal Error: trekkerFX_demux.py script terminated\n')
        LOG.flush()
    exit(1)

def validate_input():
    """
    Validates input args and if the input files are present in the current directory
    Also creates the output directory and opens the log file once the output_dir arg is known
    Returns:
        Tuple containing the paths to the input fastq.gz and sample translation file
    """
    global output_dir, LOG
    file_gz1 = ''
    file_gz2 = ''
    file_smp = ''
    if len(os.sys.argv) == 5:
        file_gz1 = os.sys.argv[1]
        file_gz2 = os.sys.argv[2]
        file_smp = os.sys.argv[3]
        output_dir = os.sys.argv[4]
    else:
        print_valid_usage_and_exit()

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f'Created output directory: {output_dir}')

    # Open log file
    log_file = os.path.join(output_dir, 'trekkerFX_demux.log')
    LOG = open(log_file, 'w')
    LOG.write(f'Created output directory: {output_dir}\n')
    print(f'Created log file: {log_file}')
    LOG.write(f'Created log file: {log_file}\n')
    LOG.flush()

    # Check if the input fastq.gz files are present in the current directory
    if not os.path.isfile(file_gz1):
        print(f'Invalid input file: {file_gz1}')
        LOG.write(f'Invalid input file: {file_gz1}')
        print_valid_usage_and_exit()
    if not os.path.isfile(file_gz2):
        print(f'Invalid input file: {file_gz2}')
        LOG.write(f'Invalid input file: {file_gz2}')
        print_valid_usage_and_exit()
    # Verify 3rd input argument
    if file_smp != 'NA' and not os.path.isfile(file_smp):
        print(f'Invalid input file: {file_smp}')
        LOG.write(f'Invalid input file: {file_smp}')
        print_valid_usage_and_exit()
    # Retrieve file names from full paths
    file_gz1_name = os.path.basename(file_gz1)
    file_gz2_name = os.path.basename(file_gz2)
    # Verify that the first file name contains 'R1' (case agnostic)
    if 'R1' not in file_gz1_name.upper():
        print('First fastq.gz file must contain "R1" or "r1" the name')
        LOG.write('First fastq.gz file must contain "R1" or "r1" the name')
        print_valid_usage_and_exit()
    # Verify that the second file name contains 'R2' (case agnostic)
    if 'R2' not in file_gz2_name.upper():
        print('Second fastq.gz file must contain "R2" or "r2" the name')
        LOG.write('Second fastq.gz file must contain "R2" or "r2" the name')
        print_valid_usage_and_exit()
    # Check if the input files end with "fastq.gz"
    if not file_gz1.endswith('fastq.gz') or not file_gz2.endswith('fastq.gz'):
        print('Input fastq files must end with "fastq.gz"')
        LOG.write('Input fastq.gz files must end with "fastq.gz"')
        print_valid_usage_and_exit()
    return file_gz1, file_gz2, file_smp

def initialize_samples(file_smp, s_labels):
    """
    Initialize a dictionary mapping barcode labels to user supplied sample names
    Also, verify that all barcode labels are present in the barcode dictionary
    Args:
        file_smp (string): path to the sample translation file or 'NA' if not provided
        This is a 2 column CSV file with user supplied sample names and barcode labels
        Column 1: Sample name
        Column 2: Barcode label
        Example: 
            Sample1,AB001
            Sample2,AB002
            Sample3,AB003
    Returns:
        Dictionary mapping barcode labels to user supplied sample names
        Example: d_smp[AB001] = Sample1
    """
    
    def sanitize_string(input_str):
        """Sanitize a string to contain only alphanumeric, underscore and hyphen characters"""
        original = input_str
        # Replace any character that's not alphanumeric, underscore or hyphen with underscore
        sanitized = re.sub(r'[^A-Za-z0-9_-]', '_', input_str)
        if sanitized != original:
            LOG.write(f"initialize_samples|WARNING: String '{original}' sanitized to '{sanitized}'\n")
            print(f"initialize_samples|WARNING: String '{original}' sanitized to '{sanitized}'")
        return sanitized
        
    d_smp = {}
    s_sample_names = set()
    c_records = 0
    c_captured = 0
    
    if file_smp != 'NA':
        with open(file_smp, 'r', encoding='utf-8-sig') as file:  # utf-8-sig handles BOM characters
            # Try different delimiters that Excel might use
            sample = file.read(1024)
            file.seek(0)
            
            delimiter = ','  # Default delimiter
            if ';' in sample:
                delimiter = ';'  # European Excel sometimes uses semicolons
            elif '\t' in sample:
                delimiter = '\t'  # Excel can also export as TSV
            
            try:
                # Use csv module to properly handle quoted fields
                reader = csv.reader(file, delimiter=delimiter)
                for row in reader:
                    if not row or row[0].startswith('#'):  # Skip empty rows and comments
                        continue
                    
                    c_records += 1
                    # Clean fields of extra whitespace
                    fields = [field.strip() for field in row]
                    
                    if len(fields) >= 2:  # Allow extra columns but require at least 2
                        barcode_label = fields[1]
                        sample_name = sanitize_string(fields[0])  # Sanitize sample name
                        
                        # Verify that barcode_label and sample_name are not empty
                        if len(barcode_label) == 0 and len(sample_name) == 0:
                            print(f"initialize_samples|WARNING: Empty barcode label and sample name in row: {','.join(fields)}")
                            LOG.write(f"initialize_samples|WARNING: Empty barcode label and sample name in row: {','.join(fields)}\n")
                            LOG.flush()
                        elif len(barcode_label) == 0 or len(sample_name) == 0:
                            print(f"initialize_samples|ERROR: Empty barcode label or sample name in row: {','.join(fields)}")
                            LOG.write(f"initialize_samples|ERROR: Empty barcode label or sample name in row: {','.join(fields)}\n")
                            LOG.flush()
                            exit(1)
                        elif barcode_label in d_smp:
                            print(f"initialize_samples|ERROR: Duplicate barcode label {barcode_label} found in sample file: {file_smp}\nScript terminated.")
                            LOG.write(f"initialize_samples|ERROR: Duplicate barcode label {barcode_label} found in sample file: {file_smp}\nScript terminated.\n")
                            LOG.flush()
                            exit(1)
                        elif barcode_label not in s_labels:
                            print(f"initialize_samples|ERROR: Barcode label {barcode_label} not expected in sample file: {file_smp}\nScript terminated.")
                            LOG.write(f"initialize_samples|ERROR: Barcode label {barcode_label} not expected in sample file: {file_smp}\nScript terminated.\n")
                            LOG.flush()
                            exit(1)
                        elif sample_name in s_sample_names:
                            print(f"initialize_samples|ERROR: Duplicate sample name {sample_name} found in sample file: {file_smp}\nScript terminated.")
                            LOG.write(f"initialize_samples|ERROR: Duplicate sample name {sample_name} found in sample file: {file_smp}\nScript terminated.\n")
                            LOG.flush()
                            exit(1)
                        else:
                            d_smp[barcode_label] = sample_name
                            s_sample_names.add(sample_name)
                            print(f"initialize_samples|INFO: Captured barcode label {barcode_label} for sample name: {sample_name}")
                            LOG.write(f"initialize_samples|INFO: Captured barcode label {barcode_label} for sample name: {sample_name}\n")
                            LOG.flush()
                            c_captured += 1
                    else:
                        print(f"initialize_samples|WARNING: Incorrect formatting in row: {row}")
                        LOG.write(f"initialize_samples|WARNING: Incorrect formatting in row: {row}\n")
                        LOG.flush()
            except csv.Error as e:
                print(f"initialize_samples|ERROR: CSV parsing error: {e}\nScript terminated.")
                LOG.write(f"initialize_samples|ERROR: CSV parsing error: {e}\nScript terminated.\n")
                LOG.flush()
                exit(1)
                
        if len(d_smp) > 0:
            print(f"initialize_samples|INFO: Captured {c_captured} records from {c_records} total in sample file: {file_smp}")
            LOG.write(f"initialize_samples|INFO: Captured {c_captured} records from {c_records} total in sample file: {file_smp}\n")
            LOG.flush()
    
    if file_smp == 'NA' or len(d_smp) == 0:
        print("initialize_samples|WARNING: No sample file provided or no valid records found. Using default barcode labels")
        LOG.write("initialize_samples|WARNING: No sample file provided or no valid records found. Using default barcode labels\n")
        LOG.flush()
    
    return d_smp

def hamming_distance(base_seq, test_seq, n):
    """determine if test_seq is within n Hamming distance of base_seq
    Requirements:
        test_seq and base_seq must be same length
    Args:
        base_seq (string): base nucleotide sequence
        test_seq (string): test nucleotide sequence
        n (int): hamming distance allowed
    Returns:
        b_match (bool): True if test_seq is within n Hamming distance of base_seq, False otherwise
    """
    b_match = False
    if len(base_seq)!= len(test_seq):
        LOG.write(f"hamming_distance|WARNING: base_seq ({base_seq}) length ({len(base_seq)}) and test_seq ({test_seq}) length ({len(test_seq)}) are not equal\n")
        LOG.flush()
    else:
        n_match = 0
        for i in range(len(base_seq)):
            if base_seq[i] == test_seq[i]:
                n_match += 1
        if n_match >= len(base_seq) - n:
            b_match = True
    return b_match

def initialize_barcodes():
    """
    Initialize barcode to label dictionary
    Returns:
        Dictionary mapping barcode sequences to their corresponding labels
    """
    d_codes = {}
    s_labels = set()
    # Add barcode sequences and their corresponding labels
    # Example: d_codes['AGGCTCCT'] = 'AB001'
    str_codes = "AGGCTCCT,AB001|GGCTCCTC,AB001|CCTTGTAG,AB001|CTTGTAGC,AB001|GTCAAGGA,AB001|TCAAGGAC,AB001|TAAGCATC,AB001|AAGCATCC,AB001|ACGACTGC,AB002|CGACTGCC,AB002|CAACACAA,AB002|AACACAAC,AB002|GTTGTACG,AB002|TTGTACGC,AB002|TGCTGGTT,AB002|GCTGGTTC,AB002|AGTCTGGA,AB003|GTCTGGAC,AB003|CAAACCTG,AB003|AAACCTGC,AB003|GTGGGAAC,AB003|TGGGAACC,AB003|TCCTATCT,AB003|CCTATCTC,AB003|AACGGAGG,AB004|ACGGAGGC,AB004|CCATTGAT,AB004|CATTGATC,AB004|GTGCATCA,AB004|TGCATCAC,AB004|TGTACCTC,AB004|GTACCTCC,AB004|AGACCTAG,AB005|GACCTAGC,AB005|CAGTAGTT,AB005|AGTAGTTC,AB005|GTCATAGC,AB005|TCATAGCC,AB005|TCTGGCCA,AB005|CTGGCCAC,AB005|AGTTCTAC,AB006|GTTCTACC,AB006|CCAAGCGA,AB006|CAAGCGAC,AB006|GACCAATT,AB006|ACCAATTC,AB006|TTGGTGCG,AB006|TGGTGCGC,AB006|ACTGAAAC,AB007|CTGAAACC,AB007|CTATGTGA,AB007|TATGTGAC,AB007|GACACCCT,AB007|ACACCCTC,AB007|TGGCTGTG,AB007|GGCTGTGC,AB007|AGCACTGG,AB008|GCACTGGC,AB008|CTTCTAAC,AB008|TTCTAACC,AB008|GAGGAGCT,AB008|AGGAGCTC,AB008|TCATGCTA,AB008|CATGCTAC,AB008|ATGCGTAG,AB009|TGCGTAGC,AB009|CCCGTGTT,AB009|CCGTGTTC,AB009|GAATCAGC,AB009|AATCAGCC,AB009|TGTAACCA,AB009|GTAACCAC,AB009|ACATTGGC,AB010|CATTGGCC,AB010|CTCACTTT,AB010|TCACTTTC,AB010|GATGACCA,AB010|ATGACCAC,AB010|TGGCGAAG,AB010|GGCGAAGC,AB010|ATCTGTTC,AB011|TCTGTTCC,AB011|CGAGTCGA,AB011|GAGTCGAC,AB011|GAGCCAAT,AB011|AGCCAATC,AB011|TCTAAGCG,AB011|CTAAGCGC,AB011|ATAGTGTC,AB012|TAGTGTCC,AB012|CCGTGAGA,AB012|CGTGAGAC,AB012|GATCCTCT,AB012|ATCCTCTC,AB012|TGCAACAG,AB012|GCAACAGC,AB012|AGTCTCAC,AB013|GTCTCACC,AB013|CAGGATTA,AB013|AGGATTAC,AB013|GCCAGACT,AB013|CCAGACTC,AB013|TTATCGGG,AB013|TATCGGGC,AB013|AACGAACA,AB014|ACGAACAC,AB014|CGTTTGGG,AB014|GTTTGGGC,AB014|GCAAGCTT,AB014|CAAGCTTC,AB014|TTGCCTAC,AB014|TGCCTACC,AB014|ATGCAGGT,AB015|TGCAGGTC,AB015|CATGCCAG,AB015|ATGCCAGC,AB015|GCATTACA,AB015|CATTACAC,AB015|TGCAGTTC,AB015|GCAGTTCC,AB015|AGTCCATC,AB016|GTCCATCC,AB016|CTAGATGG,AB016|TAGATGGC,AB016|GACTTCCA,AB016|ACTTCCAC,AB016|TCGAGGAT,AB016|CGAGGATC,AB016"
    l_codes = str_codes.split('|')
    for code in l_codes:
        split_code = code.split(',')
        d_codes[split_code[0]] = split_code[1]
        s_labels.add(split_code[1])
    print(f'Initialized barcode to label dictionary with {len(d_codes)} entries')
    LOG.write(f'Initialized barcode to label dictionary with {len(d_codes)} entries\n')
    print(f'barcode_init|INFO: Initialized set of unique labels: {len(s_labels)}')
    LOG.write(f'barcode_init|INFO: Initialized set of unique labels: {len(s_labels)}\n')
    LOG.flush()
    return d_codes, s_labels

# Class to store total sequences, exact matches, percent exact matches, no matches, and percent no matches
class ParseStats:
    def __init__(self):
        self.total_seqs = 0
        self.exact_matches = 0
        self.percent_exact_matches = 0.0
        self.no_matches = 0
        self.percent_no_matches = 0.0

    def update(self, total_seqs, exact_matches, no_matches):
        self.total_seqs = total_seqs
        self.exact_matches = exact_matches
        self.no_matches = no_matches
        if total_seqs > 0:
            self.percent_exact_matches = (exact_matches / total_seqs) * 100
            self.percent_no_matches = (no_matches / total_seqs) * 100
        else:
            self.percent_exact_matches = 0.0
            self.percent_no_matches = 0.0

def find_r2_barcode_matches(file_r2_fastq, d_codes):
    """
    Parse and extract sequence IDs and barcodes from R2.fastq.gz and match to sample labels
    Args:
        input_file (str): The path to the input file
        d_codes (dict): Dictionary mapping barcode sequences to their corresponding labels
        Example: d_codes[AGGCTCCT] = Sample1 OR AB001
    Returns:
        Dictionary mapping sequence IDs to their corresponding sample labels and set of unique labels
        d_r2[seq_id] = label
    """
    d_r2 = {}
    count_seqs = 0
    count_no_match = 0
    count_too_short = 0
    count_exact_match = 0
    #count_hamming_match = 0

    print(f'find_r2_barcode_matches|INFO: Parsing R2 input file: {file_r2_fastq}')
    LOG.write(f'find_r2_barcode_matches|INFO: Parsing R2 input file: {file_r2_fastq}\n')
    LOG.flush()
    with gzip.open(file_r2_fastq, 'rt') as f_in:
        for line in f_in:
            count_seqs += 1
            # Parse input file (4 lines per record)
            seq_id = line.strip().split(' ')[0]
            seq = f_in.readline().strip()
            plus = f_in.readline().strip()
            qual = f_in.readline().strip()

            if len(seq) < 83:
                LOG.write(f"find_r2_barcode_matches|WARNING: Sequence ID: {seq_id} has length ({len(seq)}) less than 83\n")
                count_no_match += 1
                count_too_short += 1
            else:
                # Extract barcode sequence at position 76-83
                seq_barcode = seq[75:83].upper()
                
                # Match barcodes and assign sample label
                if seq_barcode in d_codes:
                    # Exact match
                    label = d_codes[seq_barcode]
                    d_r2[seq_id] = label
                    count_exact_match += 1
                else:
                    count_no_match += 1

                    '''
                    # Search for barcode with Hamming distance of 1
                    s_ham_labels = set()
                    for barcode, ham_label in d_codes.items():
                        if hamming_distance(barcode, seq_barcode, 1):
                            s_ham_labels.add(ham_label)
                    # Hamming distance match, ONLY if there is exactly one match
                    if len(s_ham_labels) == 1:
                        label = s_ham_labels.pop()
                        d_r2[seq_id] = label
                        count_hamming_match += 1
                    # No matches OR multiple Hamming distance matches
                    else:
                    '''
                    
    parse_stats = ParseStats()
    parse_stats.update(count_seqs, count_exact_match, count_no_match)
    percent_too_short = (count_too_short / count_seqs * 100) if count_seqs > 0 else 0.0

    print(f'R2 input sequences total: {count_seqs}')
    print(f'R2 input too short: {count_too_short} and percentage: {percent_too_short:.2f}%')
    print(f'R2 input exact matches: {count_exact_match} and percentage: {parse_stats.percent_exact_matches:.2f}%')
    #print(f'R2 input Hamming matches: {count_hamming_match} and percentage: {pct_hamming_match:.2f}%')
    print(f'R2 input no matches: {count_no_match} and percentage: {parse_stats.percent_no_matches:.2f}%')
    LOG.write(f'R2 input sequences total: {count_seqs}\n')
    LOG.write(f'R2 input too short: {count_too_short} and percentage: {percent_too_short:.2f}%\n')
    LOG.write(f'R2 input exact matches: {count_exact_match} and percentage: {parse_stats.percent_exact_matches:.2f}%\n')
    #LOG.write(f'R2 input Hamming matches: {count_hamming_match} and percentage: {pct_hamming_match:.2f}%\n')
    LOG.write(f'R2 input no matches: {count_no_match} and percentage: {parse_stats.percent_no_matches:.2f}%\n')
    LOG.flush()

    if len(d_r2) == 0:
        print(f'find_r2_barcode_matches|ERROR: No barcode matches found in R2 input file: {file_r2_fastq}\nScript terminated.')
        LOG.write(f'find_r2_barcode_matches|ERROR: No barcode matches found in R2 input file: {file_r2_fastq}\nScript terminated.\n')
        LOG.flush()
        exit(1)

    return d_r2, parse_stats

def write_demux_fastqs(file_fastq_gz, d_r2, s_labels, read_pair, d_smp):
    """
    Write demultiplexed fastq.gz file based on barcode matches in d_r2
    Args:
        file_fastq_gz (str): Path to the input fastq.gz file
        d_r2 (dict): Dictionary mapping barcode sequence to their corresponding sample labels, d_r2[seq_id] = label
        s_labels (set): Set of unique sample labels
        read_pair (str): Read pair identifier, either 'R1' or 'R2
        d_smp (dict): Dictionary mapping barcode labels to user supplied sample names (d_smp[AB001] = custom_sample_name)
    Returns:
        Dictionary mapping sample labels to their counts of records
        Example: d_counts[AB001] = count
    """
    
    # Create output files for each sample label
    d_file_handles = {}
    d_buffers = {}
    for label in s_labels:
        if label in d_smp:
            file_out = os.path.join(output_dir, f'{d_smp[label]}_{read_pair}.fastq.gz')
        else:
            file_out = os.path.join(output_dir, f'{label}_{read_pair}.fastq.gz')
        d_file_handles[label] = gzip.open(file_out, 'wt', compresslevel=1)  # Use gzip with compression level 1 for faster writing
        d_buffers[label] = []
        print(f'Writing demultiplexed fastq.gz file: {file_out}')
        LOG.write(f'Writing demultiplexed fastq.gz file: {file_out}\n')
        LOG.flush()
    
    count_seqs = 0  # Count of sequences processed
    d_counts = {}   # Count of sequences per label
    for label in s_labels:
        d_counts[label] = 0
    
    # Set buffer size (number of sequence records before flushing)
    BUFFER_SIZE = 8000000

    # Function to flush a buffer to the corresponding file
    def flush_buffer(label):
        if d_buffers[label]:
            d_file_handles[label].write(''.join(d_buffers[label]))
            d_buffers[label] = []
    
    # Read input file and write demultiplexed fastq.gz files
    with gzip.open(file_fastq_gz, 'rt') as f_in:
        for line in f_in:
            count_seqs += 1
            # Periodically report progress
            if count_seqs % 5000000 == 0:
                LOG.write(f"Processed {count_seqs} sequences | Memory usage: {memory_usage():.2f} MB\n")
            # Parse input file (4 lines per sequence)
            seq_id = line.strip()
            seq_id_comp = seq_id.split(' ')[0]
            seq = f_in.readline().strip()
            plus = f_in.readline().strip()
            qual = f_in.readline().strip()
            if seq_id_comp in d_r2:
                label = d_r2[seq_id_comp]
                #d_file_handles[label].write(f'{seq_id}\n{seq}\n{plus}\n{qual}\n')
                d_buffers[label].extend([seq_id, '\n', seq, '\n', plus, '\n', qual, '\n'])
                d_counts[label] += 1

                # If buffer is full, flush it to disk
                if len(d_buffers[label]) >= BUFFER_SIZE:
                    flush_buffer(label)
    
    # Flush any remaining buffers to disk and close the files
    for label in s_labels:
        flush_buffer(label)
        d_file_handles[label].close()
    
    print(f'Total fastq records processed: {count_seqs} from parent file: {file_fastq_gz}')
    print(f'Counts of records per label:')
    for label, count in d_counts.items():
        if label in d_smp:
            print(f'{label}: {count} to file: {d_smp[label]}_{read_pair}.fastq.gz')
        else:
            print(f'{label}: {count} to file: {label}_{read_pair}.fastq.gz')
    
    LOG.write(f'Total fastq records processed: {count_seqs} from parent file: {file_fastq_gz}\n')
    LOG.write(f'Counts of records per label:\n')
    for label, count in d_counts.items():
        if label in d_smp:
            LOG.write(f'{label}: {count} to file: {d_smp[label]}_{read_pair}.fastq.gz\n')
        else:
            LOG.write(f'{label}: {count} to file: {label}_{read_pair}.fastq.gz\n')
    LOG.flush()

    return d_counts

def write_summary_file(d_counts_r1, d_counts_r2, d_smp, parse_stats):
    """
    Write summary file with stats of demultiplexing for each label
    Args:
        - d_counts_r1 (dict): Dictionary mapping sample labels to their counts of records in R1, d_counts_r1[AB001] = count
        - d_counts_r2 (dict): Dictionary mapping sample labels to their counts of records in R2, d_counts_r2[AB001] = count
        - d_smp (dict): Dictionary mapping barcode labels to user supplied sample names, d_smp[AB001] = custom_sample_name
        - parse_stats (ParseStats): Object containing parse statistics
    Output file columns: demux_summary.csv:
        1. Sample_ID = custom sample label (if provided)
        2. Barcode_ID = generic sample label
        3. Reads = Number of reads for each label
        4. Pct_reads = Percent reads representing each label
    """
    print(f'write_summary_file|INFO: Writing summary file')
    LOG.write(f'write_summary_file|INFO: Writing summary file\n')
    LOG.flush()

    # Verify R1 and R2 counts are equal
    total_count_r1 = sum(d_counts_r1.values())
    total_count_r2 = sum(d_counts_r2.values())
    if total_count_r1 != total_count_r2:
        print(f'write_summary_file|ERROR: R1 count ({total_count_r1}) and R2 count ({total_count_r2}) do not match\nScript terminated.')
        LOG.write(f'write_summary_file|ERROR: R1 count ({total_count_r1}) and R2 count ({total_count_r2}) do not match\nScript terminated.\n')
        LOG.flush()
        for label_r1, count_r1 in d_counts_r1.items():
            if d_counts_r2[label_r1] != count_r1:
                print(f'write_summary_file|ERROR: Label {label_r1} R1 count ({count_r1}) and R2 count ({d_counts_r2[label_r1]}) do not match\nScript terminated.')
                LOG.write(f'write_summary_file|ERROR: Label {label_r1} R1 count ({count_r1}) and R2 count ({d_counts_r2[label_r1]}) do not match\nScript terminated.\n')
                LOG.flush()
        exit(1)
    
    # Create and write summary file
    file_summary = os.path.join(output_dir, "metrics.csv")
    with open(file_summary, 'wt') as f_out:
        # Write parse statistics
        f_out.write(f"Total_sequences:{parse_stats.total_seqs}\n")
        f_out.write(f"Exact_matches:{parse_stats.exact_matches}\n")
        f_out.write(f"Percent_exact_matches:{parse_stats.percent_exact_matches:.2f}\n")
        f_out.write(f"No_matches:{parse_stats.no_matches}\n")
        f_out.write(f"Percent_no_matches:{parse_stats.percent_no_matches:.2f}\n")
        # Write header
        f_out.write("Sample_ID,Barcode_ID,Reads,Pct_reads\n")
        for label, count in d_counts_r1.items():
            percent = (count / total_count_r1 * 100) if total_count_r1 > 0 else 0.0
            if label in d_smp:
                f_out.write(f"{d_smp[label]},{label},{count},{percent:.2f}\n")
            else:
                f_out.write(f"{label},{label},{count},{percent:.2f}\n")

                # For labels not in d_smp (not custom) AND percent < 1, delete label_R1/2.fastq.gz files
                if percent < 1:
                    file_out = os.path.join(output_dir, f'{label}_R1.fastq.gz')
                    if os.path.isfile(file_out):
                        os.remove(file_out)
                        print(f'write_summary_file|INFO: Deleted file: {file_out}')
                        LOG.write(f'write_summary_file|INFO: Deleted file: {file_out}\n')
                        LOG.flush()
                    file_out = os.path.join(output_dir, f'{label}_R2.fastq.gz')
                    if os.path.isfile(file_out):
                        os.remove(file_out)
                        print(f'write_summary_file|INFO: Deleted file: {file_out}')
                        LOG.write(f'write_summary_file|INFO: Deleted file: {file_out}\n')
                        LOG.flush()

    print(f'Summary file created: {file_summary}')
    LOG.write(f'Summary file created: {file_summary}\n')
    LOG.flush()

def process_files_parallel(file_gz1, file_gz2, d_r2, s_labels, d_smp):
    """Process R1 and R2 files in parallel using threads"""
    # Create output directory before starting threads to avoid race conditions
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f'process_files_parallel()|created output directory: {output_dir}')
        LOG.write(f'process_files_parallel()|created output directory: {output_dir}\n')
    
    # Use ThreadPoolExecutor for I/O bound tasks
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        # Submit both tasks
        r1_future = executor.submit(write_demux_fastqs, file_gz1, d_r2, s_labels, "R1", d_smp)
        r2_future = executor.submit(write_demux_fastqs, file_gz2, d_r2, s_labels, "R2", d_smp)
        
        # Wait for both tasks to complete and get results
        d_counts_r1 = r1_future.result()
        d_counts_r2 = r2_future.result()
    
    return d_counts_r1, d_counts_r2

# Section rule and timestamped lines matching the rest of the pipeline's logs (helpers.log_ts /
# log_section on the Python side, and the equivalents in the Julia and R scripts). Only the entry and
# exit of this script are restyled: the per-function `name|LEVEL: message` lines below are a coherent
# convention of their own and are far more useful for debugging than a cosmetic rewrite would be.
def _log_section(title, width=74):
    prefix = f"── {title} "
    print("", flush=True)
    print(prefix + "─" * max(0, width - len(prefix)), flush=True)


def _log_ts(message):
    print(f"{datetime.now().strftime('%H:%M:%S')}  {message}", flush=True)


def _log_detail(message):
    print(f"{'':<10}{message}", flush=True)


if __name__ == '__main__':
    start_time_beginning = datetime.now()
    _log_section("Trekker demultiplexing")
    # Validate input args and check input files
    file_gz1, file_gz2, file_smp = validate_input()
    
    # Initialize barcode to label dictionary
    # d_codes[AGGCTCCT] = AB001
    d_codes, s_labels = initialize_barcodes()

    # Initialize label to sample dictionary
    # d_smp[AB001] = user_supplied_sample_name
    d_smp = initialize_samples(file_smp, s_labels)

    # Parse and extract sequence IDs and barcodes from R2.fastq.gz and match to sample labels
    start_time = datetime.now()
    d_r2, parse_stats = find_r2_barcode_matches(file_gz2, d_codes)
    LOG.write(f'trekkerFX_demux.py find_r2_barcode_matches() completed in time: {format_time(datetime.now() - start_time)}\n')
    LOG.flush()

    # Write demultiplexed fastq.gz files
    start_time = datetime.now()
    d_counts_r1, d_counts_r2 = process_files_parallel(file_gz1, file_gz2, d_r2, s_labels, d_smp)
    LOG.write(f'trekkerFX_demux.py write_demux_fastqs (parallel) completed in time: {format_time(datetime.now() - start_time)}\n')
    LOG.flush()

    # Write summary file
    write_summary_file(d_counts_r1, d_counts_r2, d_smp, parse_stats)

    elapsed = format_time(datetime.now() - start_time_beginning)
    LOG.write(f'trekkerFX_demux.py completed in time: {elapsed}\n')
    LOG.flush()

    # the LOG file records this too, but the pipeline captures stdout into takara_pipeline.log, which
    # is where anyone looking for "did the demux finish, and how long did it take" will look first
    _log_ts(f"✓ demultiplexing complete — {elapsed}")

    LOG.close()
    sys.exit(0)
