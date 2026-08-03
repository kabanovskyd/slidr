
# Import modules
import sys
import os
import csv
from Bio import SeqIO
import gzip
import psutil
import gc
import time
import statistics

# Used for testing only
s_test_cb = set()
s_cb1 = set()
f_test_cb = ""
samplesheet_override = ""
whitelist_override = ""

print(f'START|fastq_parser.py', flush=True)
def get_memory():
    """calculate current memory usage
    Returns:
        (string): memory usage in B, KB, MB, GB
    """
    process = psutil.Process(os.getpid())
    memory = process.memory_info().rss
    if memory < 1000:
        return f"{memory:.2f}B"
    elif memory < 1_000_000:
        return f"{(memory/1000):.2f}KB"
    elif memory < 1_000_000_000:
        return f"{(memory/1000/1000):.2f}MB"
    else:
        return f"{(memory/1000/1000/1000):.2f}GB"

def swap_fastqs(fastq_R1, fastq_R2):
    """Returns swapped R1 and R2 fastqs if sc_platform is 'TrekkerQ_P' and 'TrekkerQ_S'."""
    return fastq_R2, fastq_R1

def cleanup_barcodes_TrekkerFX(input_file, output_file):
    with gzip.open(input_file, "rt") as infile, gzip.open(output_file, "wt") as outfile:
        for line in infile:
            barcode = line[:16]
            if len(line.rstrip("\n")) < 16:
                print(f"cleanup_barcodes_TrekkerFX|[ERROR]: line in {input_file} is shorter than the expected 16bp barcode: {line!r}", flush=True)
                print("Troubleshooting:", flush=True)
                print(" \u2022 This file holds the cell barcodes read back from the single-cell matrix, one 16bp barcode per line", flush=True)
                print(" \u2022 A malformed line means the upstream count matrix is not the 10x Flex layout this step expects", flush=True)
                print(f" \u2022 Inspect it: `zcat {input_file} | head -3`", flush=True)
                print(" \u2022 Re-run the count stage with --count --force, then --spatial-analysis --force", flush=True)
                sys.exit(1)
            if "-" in barcode:
                print(f"cleanup_barcodes_TrekkerFX|[ERROR]: unexpected '-' within first 16bp of barcode in {input_file}: {line!r}", flush=True)
                print("Troubleshooting:", flush=True)
                print(" \u2022 A '-' inside the first 16 bases means the lane suffix (e.g. '-1') has already been stripped or the barcodes are shorter than 16bp", flush=True)
                print(f" \u2022 Inspect it: `zcat {input_file} | head -3`", flush=True)
                print(" \u2022 Re-run the count stage with --count --force, then --spatial-analysis --force", flush=True)
                sys.exit(1)
            if "1" in barcode:
                print(f"cleanup_barcodes_TrekkerFX|[ERROR]: unexpected '1' within first 16bp of barcode in {input_file}: {line!r}", flush=True)
                print("Troubleshooting:", flush=True)
                print(" \u2022 A digit inside the first 16 bases means this is not a raw ACGT barcode -- the lane suffix has bled into the barcode field", flush=True)
                print(f" \u2022 Inspect it: `zcat {input_file} | head -3`", flush=True)
                print(" \u2022 Re-run the count stage with --count --force, then --spatial-analysis --force", flush=True)
                sys.exit(1)
            outfile.write(barcode + "-1\n")
    return output_file

def get_rhapsody_keys():
    """initialize keys for Rhapsody
    Returns:
        Enh_linker1: Rhapsody enhancer linker 1
        Enh_linker2: Rhapsody enhancer linker 2
        B384_cell_key1: Rhapsody cell barcode key 1
        B384_cell_key2: Rhapsody cell barcode key 2
        B384_cell_key3: Rhapsody cell barcode key 3
    """
    Enh_linker1 = 'GTGA'
    Enh_linker2 = 'GACA'
    B384_cell_key1 = ("TGTGTTCGC","TGTGGCGCC","TGTCTAGCG","TGGTTGTCC","TGGTTCCTC","TGGTGTGCT","TGGCGACCG","TGCTGTGGC","TGCTGGCAC","TGCTCTTCC",
                    "TGCCTCACC","TGCCATTAT","TGATGTCTC","TGATGGCCT","TGATGCTTG","TGAAGGACC","TCTGTCTCC","TCTGATTAT","TCTGAGGTT","TCTCGTTCT",
                    "TCTCATCCG","TCCTGGATT","TCAGCATTC","TCACGCCTT","TATGTGCAC","TATGCGGCC","TATGACGAG","TATCTCGTG","TATATGACC","TAGGCTGTG",
                    "TACTGCGTT","TACGTGTCC","TAATCACAT","GTTGTGTTG","GTTGTGGCT","GTTGTCTGT","GTTGTCGAG","GTTGTCCTC","GTTGTATCC","GTTGGTTCT",
                    "GTTGGCGTT","GTTGGAGCG","GTTGCTGCC","GTTGCGCAT","GTTGCAGGT","GTTGCACTG","GTTGATGAT","GTTGATACG","GTTGAAGTC","GTTCTGTGC",
                    "GTTCTCTCG","GTTCTATAT","GTTCGTATG","GTTCGGCCT","GTTCGCGGC","GTTCGATTC","GTTCCGGTT","GTTCCGACG","GTTCACGCT","GTTATCACC",
                    "GTTAGTCCG","GTTAGGTGT","GTTAGAGAC","GTTAGACTT","GTTACCTCT","GTTAATTCC","GTTAAGCGC","GTGTTGCTT","GTGTTCGGT","GTGTTCCAG",
                    "GTGTTCATC","GTGTCACAC","GTGTCAAGT","GTGTACTGC","GTGGTTAGT","GTGGTACCG","GTGGCGATC","GTGCTTCTG","GTGCGTTCC","GTGCGGTAT",
                    "GTGCGCCTT","GTGCGAACT","GTGCAGCCG","GTGCAATTG","GTGCAAGGC","GTCTTGCGC","GTCTGGCCG","GTCTGAGGC","GTCTCAGAT","GTCTCAACC",
                    "GTCTATCGT","GTCGGTGTG","GTCGGAATC","GTCGCTCCG","GTCCTCGCC","GTCCTACCT","GTCCGCTTG","GTCCATTCT","GTCCAATAC","GTCATGTAT",

                    "GTCAGTGGT","GTCAGATAG","GTATTAACT","GTATCAGTC","GTATAGCCT","GTATACTTG","GTATAAGGT","GTAGCATCG","GTACCGTCC","GTACACCTC",
                    "GTAAGTGCC","GTAACAGAG","GGTTGTGTC","GGTTGGCTG","GGTTGACGC","GGTTCGTCG","GGTTCAGTT","GGTTATATT","GGTTAATAC","GGTGTACGT",
                    "GGTGCCGCT","GGTGCATGC","GGTCGTTGC","GGTCGAGGT","GGTAGGCAC","GGTAGCTTG","GGTACATAG","GGTAATCTG","GGCTTGGCC","GGCTTCACG",
                    "GGCTTATGT","GGCTTACTC","GGCTGTCTT","GGCTCTGTG","GGCTCCGGT","GGCTCACCT","GGCGTTGAG","GGCGTGTAC","GGCGTGCTG","GGCGTATCG",
                    "GGCGCTCGT","GGCGCTACC","GGCGAGCCT","GGCGAGATC","GGCGACTTG","GGCCTCTTC","GGCCTACAG","GGCCAGCGC","GGCCAACTT","GGCATTCCT",
                    "GGCATCCGC","GGCATAACC","GGCAACGAT","GGATGTCCG","GGATGAGAG","GGATCTGGC","GGATCCATG","GGATAGGTT","GGAGTCGTG","GGAGAAGGC",
                    "GGACTCCTT","GGACTAGTC","GGACCGTTG","GGAATTAGT","GGAATCTCT","GGAATCGAC","GGAAGCCTC","GCTTGTAGC","GCTTGACCG","GCTTCGGAC",
                    "GCTTCACAT","GCTTAGTCT","GCTGGATAT","GCTGGAACC","GCTGCGATG","GCTGATCAG","GCTGAGCGT","GCTCTTGTC","GCTCTCCTG","GCTCGGTCC",
                    "GCTCCAATT","GCTATTCGC","GCTATGAGT","GCTAGTGTT","GCTAGGATC","GCTAGCACT","GCTACGTAT","GCTAACCTT","GCGTTCCGC","GCGTGTGCC",
                    "GCGTGCATT","GCGTCGGTT","GCGTATGTG","GCGTATACT","GCGGTTCAC","GCGGTCTTG","GCGGCGTCG","GCGGCACCT","GCGCTGGAC","GCGCTCTCC",

                    "GCGCGGCAG","GCGCGATAC","GCGCCGACC","GCGAGCGAG","GCGAGAGGT","GCGAATTAC","GCCTTGCAT","GCCTGCGCT","GCCTAACTG","GCCGTCCGT",
                    "GCCGCTGTC","GCCATGCCG","GCCAGCTAT","GCCAACCAG","GCATGGTTG","GCATCGACG","GCAGGCTAG","GCAGGACGC","GCAGCCATC","GCAGATACC",
                    "GCAGACGTT","GCACTATGT","GCACACGAG","GATTGTCAT","GATTGGTAG","GATTGCACC","GATTCTACT","GATTCGCTT","GATTAGGCC","GATTACGGT",
                    "GATGTTGGC","GATGTTATG","GATGGCCAG","GATCGTTCG","GATCGGAGC","GATCGCCTC","GATCCTCTG","GATCCAGCG","GATACACGC","GAGTTACCT",
                    "GAGTCGTAT","GAGTCGCCG","GAGGTGTAG","GAGGCATTG","GAGCGGACG","GAGCCTGAG","GAGATCTGT","GAGATAATT","GAGACGGCT","GACTTCGTG",
                    "GACTGTTCT","GACTCTTAG","GACCGCATT","GAATTGAGC","GAATATTGC","GAAGGCTCT","GAAGAGACT","GAACTGCCG","GAACGCGTG","CTTGTGTAT",
                    "CTTGTGCGC","CTTGTCATG","CTTGGTCTT","CTTGGTACC","CTTGGATGT","CTTGCTCAC","CTTGCAATC","CTTGAGGCC","CTTGACGGT","CTTCTGATC",
                    "CTTCTCGTT","CTTCTAGGC","CTTCGTTAG","CTTATGTCC","CTTATGCTT","CTTATATAG","CTTAGGTTG","CTTAGGAGC","CTTACTTAT","CTGTTCTCG",
                    "CTGTGCCTC","CTGTCGCAT","CTGTCGAGC","CTGTAGCTG","CTGTACGTT","CTGCTTGCC","CTGCGTAGT","CTGCACACC","CTGATGGAT","CTGAGTCAT",
                    "CTGACGCCG","CTGAACGAG","CTCTTGTAG","CTCTTAGTT","CTCTTACCG","CTCTGCACC","CTCTCGTCC","CTCGTATTG","CTCGACTAT","CTCCTGACG",

                    "CTCACTAGC","CTATACGGC","CGTTCGCTC","CGTTCACCG","CGTATAGTT","CGGTGTTCC","CGGTGTCAG","CGGTCCTGC","CGGCGACTC","CGGCACGGT",
                    "CGGATAGCC","CGGAGAGAT","CGCTAATAG","CGCGTTGGC","CGCGCAGAG","CGCACTGCC","CCTTGTCTC","CCTTGGCGT","CCTTCTGAG","CCTTCTCCT",
                    "CCTTCGACC","CCTTACTTG","CCTGTTCGT","CCTGTATGC","CCTCGGCCG","CCGTTAATT","CCATGTGCG","CCAGTGGTT","CCAGGCATT","CCAGGATCC",
                    "CCAGCGTTG","CATTCCGAT","CATTATACC","CATGTTGAG","ATTGCGTGT","ATTGCGGAC","ATTGCGCCG","ATTGACTTG","ATTCGGCTG","ATTCGCGAG",
                    "ATTCCAAGT","ATTATCTTC","ATTACTGTT","ATTACACTC","ATGTTCTAT","ATGTTACGC","ATGTGTATC","ATGTGGCAG","ATGTCTGTG","ATGGTGCAT",
                    "ATGCTTACT","ATGCTGTCC","ATGCTCGGC","ATGAGGTTC","ATGAGAGTG","ATCTTGGCT","ATCTGTGCG","ATCGGTTCC","ATCATGCTC","ATCATCACT",
                    "ATATCTTAT","ATAGGCGCC","AGTTGGTAT","AGTTGAGCC","AGTGCGACC","AGGTGCTAC","AGGCTTGCG","AGGCCTTCC","AGGCACCTT","AGGAATATG",
                    "AGCGGCCAG","AGCCTGGTC","AGCCTGACT","AGCAATCCG","AGAGATGTT","AGAGAATTC","ACTCGCTTG","ACTCGACCT","ACGTACACC","ACGGATGGT",
                    "ACCAGTCTG","ACATTCGGC","ACATGAGGT","ACACTAATT"
                    )

    B384_cell_key2 = ("TTGTGTTGT","TTGTGGTAG","TTGTGCGGA","TTGTCTGTT","TTGTCTAAG","TTGTCATAT","TTGTCACGA","TTGTATGAA","TTGTACAGT","TTGGTTAAT",
                    "TTGGTGCAA","TTGGTCGAG","TTGGTATTA","TTGGCACAG","TTGGATACA","TTGGAAGTG","TTGCGGTTA","TTGCCATTG","TTGCACGCG","TTGCAAGGT",
                    "TTGATGTAT","TTGATAATT","TTGAGACGT","TTGACTACT","TTGACCGAA","TTCTGGTCT","TTCTGCACA","TTCTCCTTA","TTCTCCGCT","TTCTAGGTA",
                    "TTCTAATCG","TTCGTCGTA","TTCGTAGAT","TTCGGCTTG","TTCGGAATA","TTCGCCAGA","TTCGATTGT","TTCGATCAG","TTCCTCGGT","TTCCGGCAG",
                    "TTCCGCATT","TTCCAATTA","TTCATTGAA","TTCATGCTG","TTCAGGAGT","TTCACTATA","TTCAACTCT","TTCAACGTG","TTATGCGTT","TTATGATTG",
                    "TTATCCTGT","TTATCCGAG","TTATATTAT","TTAGGCGCG","TTACTGGAA","TTACTAGTT","TTACGTGGT","TTACGATAT","TTACCTAGA","TTACATGAG",
                    "TTACAGCGT","TTACACGGA","TTACACACT","TTAATCAGT","TTAATAGGA","TTAAGTGTG","TTAACCTTG","TTAACACAA","TGTTCACTT","TGTTCAAGA",
                    "TGTTAAGTG","TGTGTTATG","TGTGTCCAA","TGTGGAGCG","TGTCAGTTA","TGTCAGAAG","TGGTTAGTT","TGGTTACAA","TGGCGTTAT","TGGCGCCAA",
                    "TGGAGTCTT","TGCGTATTG","TGATAGAGA","TGAGGTATT","TGAGAATCT","TCTTGGTAA","TCTTCATAG","TCTGTCCTT","TCTGGAATT","TCTACCGCG",
                    "TCGTTCGAA","TCGTCAGTG","TCGACGAGA","TCATGGCTT","TCACACTTA","TATTCCGAA","TATTATGGT","TATGCTATT","TATCAAGGA","TAGTTCAAT",

                    "TAGCTGCTT","TAGAGGAAG","TACCTGTTA","TACACCTGT","GTTGTGCGT","GTTGGCTAT","GTTGCCAAG","GTTGACCTT","GTTCTGCTA","GTTCTGAAT",
                    "GTTCTATCA","GTTCGCGTG","GTTCCTTAT","GTTAGCAGT","GTTACTGTG","GTTACTCAA","GTTAAGAGA","GTTAACTTA","GTGTCGGCA","GTGTCCATT",
                    "GTGCTTGAG","GTGCTCGTT","GTGCTCACA","GTGCCTGGA","GTCTTGTCG","GTCTTGATT","GTCTTCCGT","GTCTTAAGA","GTCTCATCT","GTCTACGAG",
                    "GTCGTTGCT","GTCGTGTTA","GTCGGTAAT","GTCGGATGT","GTCGAGCTG","GTCCGGACT","GTCCAACAT","GTCAGACGA","GTCAGAATT","GTCACTCTT",
                    "GTCAAGGAA","GTATGTCTT","GTATGTACA","GTATCGGTT","GTATATGTA","GTATACAAT","GTAGTTAAG","GTAGTCGAT","GTAGCCTTA","GTAGATACT",
                    "GTACGATTA","GTACAGTCT","GTAATTCGT","GCTTGGCAG","GCTTGCTTG","GCTTGAGGA","GCTTCATTA","GCTTATGCG","GCTGTGTAG","GCTGTCATG",
                    "GCTGGTTGT","GCTGGACTG","GCTGCCTAA","GCTGATATT","GCTCTTAGT","GCTCTATTG","GCTCGCCGT","GCTCCGCTG","GCTATTCTG","GCTATACGA",
                    "GCTACTAAG","GCTACATGT","GCTAACTCT","GCGTTGTAA","GCGTTCTCT","GCGTGCGTA","GCGTCTTGA","GCGTCCGAT","GCGTAAGAG","GCGCTTACG",
                    "GCGCGGATT","GCGCCATAT","GCGCATGAA","GCGATCAAT","GCGAGCCTT","GCGAGATTG","GCGAGAACA","GCCTTGGTA","GCCTTCTAG","GCCTTCACA",
                    "GCCTGAGTG","GCCTCACGT","GCCGGCGAA","GCCGCACAA","GCCATGCTT","GCCATATAT","GCCAATTCG","GCATTCGTT","GCATGATGT","GCAGTTGGA",

                    "GCAGTGTCT","GCACTTGTG","GCAATCTGT","GCAACACTT","GATTGTATT","GATTGCGAG","GATTCCAGT","GATTCATAT","GATTATCAG","GATTAGGTT",
                    "GATGTTGCG","GATGGATCT","GATGCTGAT","GATGCCTTG","GATCTCCTT","GATCGCTTA","GATATTGAA","GATATTACT","GAGTGTTAT","GAGCTCAGT",
                    "GAGCGTGCT","GAGCGTCGA","GAGCGGTTG","GAGCGACTT","GAGCCGAAT","GAGATAGAT","GAGACCTAT","GACGGTCGT","GACGCAGGT","GACGATATG",
                    "GACCTATCT","GAATTAGGA","GAATCAGCT","GAAGTTCAT","GAAGTGGTT","GAAGTATTG","GAAGGCATT","GAACGCTGT","CTTGTCCAG","CTTGGATTG",
                    "CTTGCTGAA","CTTGCCGTG","CTTGATTCT","CTTCTGTCG","CTTCGGCGT","CTTATGAGT","CTTACCGAT","CTGTTAGGT","CTGTCGTCT","CTGTATAAT",
                    "CTGGCTCAT","CTGGATGCG","CTGCGTGTG","CTGCGCGGT","CTGCCGATT","CTGCATTGT","CTGATTAAG","CTGAGATAT","CTGACCTGT","CTCGTATCT",
                    "CTCGGCAAG","CTCGCAATT","CTCCTGCTT","CTCCTAAGT","CTCCGGATG","CTCCGAGCG","CTCACAGGT","CTATTCTAT","CTATTAGTG","CTATGAATT",
                    "CTACATATT","CGTGGCATT","CGTCTTAAT","CGTCTGGTT","CGTCACTGT","CGTAGGTCT","CGGTTCGAG","CGGTTCATT","CGGTGCTCT","CGGTAATTG",
                    "CGGCCTGAT","CGGATATAG","CGGAATATT","CGCTCCAAT","CGCGTTCGT","CGCAGGTTG","CGAGGATGT","CGAGCTGTT","CGACGGCTT","CCTTGTGTG",
                    "CCTGTCTCA","CCTGACTAT","CCTACCTTG","CCGTAGATT","CCGGCTGGT","CATCGGACG","CATCGATAA","CATCCTTCT","CAGTTCTGT","CAGTGCCAG",

                    "CAGGCACTG","CAGCCTCTT","CACTTATAT","CACTGGTCG","CACTGCATG","CACGCGTTG","CACGATGTT","CACCATCTG","CACAGGCGT","ATTGTACAA",
                    "ATTGGTATG","ATTGCTAAT","ATTGCATAG","ATTGCAGTT","ATTCTGCAG","ATTCTACGT","ATTCGGATT","ATTCCGTTG","ATTCATCAA","ATTCAAGAG",
                    "ATTAGCCTT","ATTAATATT","ATGTTAGAG","ATGTTAACT","ATGTAGTCG","ATGGTGTAG","ATGGATTAT","ATCTTGAAG","ATCTGATAT","ATCTCAGAA",
                    "ATCGCTCAA","ATCGCGTCG","ATCCATGGT","ATCATGAGA","ATCATAGTT","ATCAGCGAG","ATCACCATT","ATAGTAATT","ATAGCTGTG","ATACTCTCG",
                    "ATACCTCAT","AGTTGCGCG","AGTTGAATT","AGTTATGAT","AGTGTCCGT","AGTGGCTTG","AGTGCTTCT","AGTATCATT","AGTACACAA","AGGTATGCG",
                    "AGGTATAGT","AGGCTACTT","AGGCCAGGT","AGGAGCGAT","AGCTTATAG","AGCTCTAGA","AGCGTGTAT","AGCGTCACA","AGCCTTCAT","AGCCTGTCG",
                    "AGCCTCGAG","AGCACTGAA","AGATGTACG","AGAGTTAAT","AGACCTCTG","ACTTCTATA","ACTGTCGAG","ACTGTATGT","ACTCTGTAA","ACTCGCGAA",
                    "ACTAGATCT","ACTAACGTT","ACGTTACTG","ACGTGGAAT","ACGGACTCT","ACGCCTAAT","ACGCCGTTA","ACGACGTGT","ACCTCGCAT","ACCATCATA",
                    "ACATATATT","ACAGGCACA","ACACCTGAG","ACACATTCT"
                    )

    B384_cell_key3 = ("TTGTGGCTG","TTGTGGAGT","TTGTGCGAC","TTGTCTTCA","TTGTAAGAT","TTGGTTCTG","TTGGTGCGT","TTGGTCTAC","TTGGTAACT","TTGGCGTGC",
                    "TTGGATTAG","TTGGAGACG","TTGGAATCA","TTGCGGCGA","TTGCGCTCG","TTGCCTTAC","TTGCCGGAT","TTGCATGCT","TTGCACGTC","TTGCACCAT",
                    "TTGAACCTG","TTCTCGCGT","TTCTCAACT","TTCTACTCA","TTCGTCCAT","TTCGGATAC","TTCGGACGT","TTCGCAATC","TTCCGGTGC","TTCCGACTG",
                    "TTCATTATG","TTCATGGAT","TTCAGCGCA","TTCACCTCG","TTCAAGCAG","TTCAACTAC","TTATGCCAG","TTATGCATC","TTATCGTAC","TTATACCTA",
                    "TTATAATAG","TTATAAGTC","TTAGTTAGC","TTAGCTCAT","TTAGCACTA","TTAGATATG","TTACTACGA","TTACCGTCA","TTACAGAGC","TTAATTGCA",
                    "TTAACAGAT","TGTTGGCTA","TGTTGATGA","TGTTAAGCT","TGTGGCCGA","TGTGCTAGC","TGTGCGTCA","TGTCGCAGT","TGTCGAGCA","TGTACAACG",
                    "TGGTTCCGA","TGGTTCACT","TGGTCAAGT","TGGCTTGTA","TGGCTGTCG","TGGCGTATG","TGGCGCGCT","TGGATGTAC","TGGACTTGC","TGGAATACT",
                    "TGCTAGCGA","TGCGTTGCT","TGCGGTCTG","TGCGCTTAG","TGCGCGACG","TGCCTGCAT","TGCCTAGAC","TGCACGAGT","TGAGTGTGC","TGAGGCTCG",
                    "TCTTCCGTC","TCTTATAGT","TCTTACCAT","TCTGTTGTC","TCTGTTACT","TCTGGCTAG","TCTCAGATC","TCTAGTTGA","TCTAGTACG","TCGTACTAC",
                    "TCGGTGTAG","TCGGCTGCT","TCGCTACTG","TCGATCACG","TCGAGGCAT","TCCGGCGTC","TCCGGAGCT","TCCGCTCGT","TCCGAGTAC","TCCATTCAT",

                    "TCCATGGTC","TCCAAGTCG","TCATTACGT","TCATGCACT","TCAGGTTGC","TCAGACCGT","TCACTCAGT","TCAAGCTCA","TATTGCGCA","TATTCGGCT",
                    "TATTCCAGC","TATTCATCA","TATGTTCAG","TATGGTATG","TATGCAAGT","TATCTGGTC","TATCTGACT","TATCCAGAT","TATCAGTCG","TATCACGCT",
                    "TAGGCGCGA","TAGGCACAT","TAGGATCGT","TAGCATTGC","TAGAGTTAC","TAGACTGAT","TACTTGTCG","TACGTCCGA","TACCGTACT","TACCGCGAT",
                    "TACCAGGAC","TACAGAAGT","TAAGTGCAT","TAAGCTACT","GTTGACCGA","GTTCTCGAC","GTTCCTGCT","GTTATGATG","GTGCTTGCA","GTGCCGCGT",
                    "GTATTGCTG","GTATTCCGA","GTATTAAGC","GTATGACGT","GTAGTTGTC","GTAGTACAT","GTAGCTCGA","GGTTGCTCA","GGTTGAGTA","GGTTAACGT",
                    "GGTGTGGCA","GGTCTTCAG","GGTCGTCTA","GGTCGGCGT","GGTCCGACT","GGTCATGTC","GGTCACATG","GGTAGTGCT","GGTAGCGTC","GGTACCAGT",
                    "GGTAAGGAT","GGCTTGTGC","GGCTTGACT","GGCTTACGA","GGCTGTAGT","GGCTGGCAG","GGCTCCATC","GGCGTGGAT","GGCGTAATC","GGCGCAAGT",
                    "GGCGAGTAG","GGCGACCGT","GGCCTGTCA","GGCCATTGC","GGCACTCTG","GGATGTCAT","GGAGTAACT","GGAGAACGA","GGACTGGCT","GGACGTTCA",
                    "GGAACGTGC","GCTGTCCAT","GCTGGTTCA","GCTGCAACT","GCTCGTTAC","GCTATAGAT","GCTAGTCGT","GCTACCATG","GCGTTCTGA","GCGTGTTAG",
                    "GCGGTATCG","GCGGAGCAT","GCGCGGTGC","GCGCCTAGT","GCGCCGGCT","GCCTTCATG","GCCATACTG","GCATGTTGA","GCATGCTAC","GCAGTATAC",

                    "GCAGGTACT","GCAGCGCGT","GCACCTCAT","GCAATTCGA","GATTGCCGT","GATGAACAT","GATCTTCGA","GATCTGCAT","GAGTGGCAT","GAGTCGGAC",
                    "GAGTATGAT","GAGGCGAGT","GAGGCAACG","GAGCGCACT","GAATAGGCT","ATTGTCACT","ATTGTATCA","ATTGGTCAG","ATTGGCGAT","ATTGATCGT",
                    "ATTCGTAGT","ATTCATACG","ATTCAGGAC","ATTACTTCA","ATTAATTAG","ATTAAGCAT","ATGTCTCTA","ATGTAGCGT","ATGGCATAC","ATGGAGATC",
                    "ATGGACTCG","ATGGAACGA","ATGCTTCAT","ATGCTCGCT","ATGCGACGT","ATGCCGTAG","ATGAGTTCG","ATGACTATC","ATGACCGAC","ATCTTATGC",
                    "ATCTTACTA","ATCTATCAG","ATCGTGTAC","ATCGTCTGA","ATCGGCATG","ATCGCGAGC","ATCGCAACG","ATCGATGCT","ATCGAATAG","ATCCTTCTG",
                    "ATCCTGCGT","ATCCGCACT","ATCCATTAC","ATCCAAGCA","ATCAGATCA","ATCACACAT","ATCAACGTC","ATCAACCGA","ATATTGAGT","ATATTCGTC",
                    "ATATTACAG","ATATCTTGA","ATATCGCAT","ATATCAATC","ATAGTCCTG","ATAGGTCTA","ATAGCTGAC","ATAGCGGTA","AGTTCGCTG","AGTTACAGC",
                    "AGTTAACTA","AGTGCAATC","AGTCTGGTA","AGTCTGAGC","AGTCTACAT","AGTCGAACT","AGTCCATCG","AGTCATTCA","AGTATCCAG","AGTAGACTG",
                    "AGTAATCGA","AGTAAGTGC","AGGTTGGCT","AGGTTCTAG","AGGTGTTCA","AGGTGCCAT","AGGTCTGAT","AGGTCGTAC","AGGTCAGCA","AGGCTTATC",
                    "AGGCTATGA","AGGCCGACG","AGGCCAAGC","AGGCAGGTC","AGGCAAGAT","AGGAGCAGT","AGGACCGCT","AGGAATTAC","AGCTTGGAC","AGCTTAAGT",

                    "AGCTACACG","AGCGTTACG","AGCGGTGCA","AGCGGAGTC","AGCGGACGA","AGCGCGCTA","AGCGATAGC","AGCGACTCA","AGCCTCTAC","AGCCGTCGT",
                    "AGCATGATC","AGCACTTCG","AGCACGGCA","AGATTCTGA","AGATTAGAT","AGATGATAG","AGATATGTA","AGATACCGT","AGAGTGCGT","AGAGCCGAT",
                    "AGACTCACT","ACTTGCCTA","ACTTGAGCA","ACTTCTAGC","ACTTCGACT","ACTTAGTAC","ACTGTTGAT","ACTGTAACG","ACTGGTATC","ACTGACGTC",
                    "ACTGAAGCT","ACTCTGATG","ACTCCTGAC","ACTCCGCTA","ACTCAACTG","ACTATTGCA","ACTAGGCAG","ACTACGCGT","ACTAATACT","ACGTTCGTA",
                    "ACGTGTGCT","ACGTGTATG","ACGTGGAGC","ACGTCTTCG","ACGTCAGTC","ACGGTCTCA","ACGGTCCGT","ACGGTACAG","ACGGCGCTG","ACGCTGCGA",
                    "ACGCGTGTA","ACGCGCCAG","ACGATGTCG","ACGATGGAT","ACGATCTAC","ACGAGCTGA","ACGAGCATC","ACGAATCGT","ACGAACGCA","ACCTTGTAG",
                    "ACCTGTTGC","ACCTGTCAT","ACCTCGATC","ACCTAGGTA","ACCTACTGA","ACCTAATCG","ACCGTAGCA","ACCGGTAGT","ACCGGCTAC","ACCGCTTCA",
                    "ACATTGTGC","ACATTCTCG","ACATGGCTG","ACATGACGA","ACATATGAT","ACATATACG","ACAGCGTAC","ACACTTGCT","ACACTATCA","ACACGCATG",
                    "ACACCAGTA","ACACCAACT","ACACATAGT","ACACACCTA"
                    )
    return Enh_linker1, Enh_linker2, B384_cell_key1, B384_cell_key2, B384_cell_key3

def index_to_sequence(index, Enh_linker1, Enh_linker2, B384_cell_key1, B384_cell_key2, B384_cell_key3):
    """convert numerical index to cell barcode sequence for bead version EnhV2
    Args:
        index (int): Rhapsody cell barcode index
    Returns:
        seq (string): full sequence including linkers and cell barcode segments
        cb (string): 27nt cell barcode sequence
    """
    zerobased = int(index) - 1
    cl1 = (int((zerobased) / 384 / 384) % 384) + 1
    cl2 = (int((zerobased) / 384) % 384) + 1
    cl3 = (zerobased % 384) + 1
    diversityInsert = ''
    subIndex = ((cl1-1) % 96) + 1

    if 1 <= subIndex <= 24:
        diversityInsert = ''
    elif 25 <= subIndex <= 48:
        diversityInsert = 'A'
    elif 49 <= subIndex <= 72:
        diversityInsert = 'GT'
    else: # 73 <= subIndex <= 96:
        diversityInsert = 'TCA'

    cls1_sequence = B384_cell_key1[cl1-1]
    cls2_sequence = B384_cell_key2[cl2-1]
    cls3_sequence = B384_cell_key3[cl3-1]

    seq = diversityInsert + cls1_sequence + Enh_linker1 + cls2_sequence + Enh_linker2 + cls3_sequence
    cb = cls1_sequence + cls2_sequence + cls3_sequence
    return seq,cb


def get_cb_umi(seq, sc_platform):
    """get cell barcode (cb) and umi from read 1 sequence
    Args:
        seq (string): sequence of read 1
        sc_platform (string): One of ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_C', 'TrekkerU_CX', 'TrekkerU_M','TrekkerR', 'Rhapsody', 'TrekkerU_R', 'TrekkerFX_FLEX', 'TrekkerQ_S', 'TrekkerQ_P', 'TrekkerU_IL', 'TrekkerU_PIP', 'TrekkerSHA_WGA', 'TrekkerU_F_TEST', 'TrekkerU_RATAC', 'TrekkerU_RVDJ')
    Returns:
        cb (string): cell barcode sequence
        umi (string): umi sequence
    """
    cb = ""
    umi = ""
    if sc_platform in ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_C', 'TrekkerU_CX', 'TrekkerU_M', 'TrekkerFX_FLEX', 'Trekker5C_CX', 'Trekker5C_C'):
        if len(seq) >= 28:
            cb = seq[0:16]
            umi = seq[16:28]
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short {seq}", flush=True)
    elif sc_platform in ('TrekkerU_PIP', 'TrekkerSHA_WGA'):
        if len(seq) >= 16:
            cb = seq[0:16]
            # UMI is on R2 for TrekkerU_PIP (Fluent)
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short {seq}", flush=True)
    elif sc_platform in ('TrekkerU_IL',):
        if len(seq) >= 39:
            linker = seq[26:34]
            idx = linker.find("TCGAG")
            if idx >= 0:
                if 39+idx > len(seq):
                    print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short (<{39+idx}) for seq: {seq}", flush=True)
                else:
                    t1 = seq[0+idx:8+idx]
                    t2 = seq[11+idx:17+idx]
                    t3 = seq[20+idx:26+idx]
                    t4 = seq[31+idx:39+idx]
                    cb = t1 + t2 + t3 + t4
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short {seq}", flush=True)
    elif sc_platform in ('TrekkerQ_P',):
        if len(seq) >= 58:
            cb = seq[10:18] + seq[30:38] + seq[50:58]     # 8+8+8 = 24
            umi = seq[0:10]
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short {seq}", flush=True)
    elif sc_platform in ('TrekkerQ_S',):
        if len(seq) >= 48:
            cb = seq[0:40]
            umi = seq[40:48]
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short {seq}", flush=True)
    elif sc_platform in ('TrekkerR', 'Rhapsody'):
        if len(seq) >= 63:
            cb = seq[20:29] + seq[33:42] + seq[46:55]     # 9+9+9 = 27
            umi = seq[55:63]
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short {seq}", flush=True)
    elif sc_platform in ('TrekkerU_R', 'TrekkerU_RATAC', 'TrekkerU_RVDJ'):
        if len(seq) >= 43:
            linker = seq[9:16]
            idx = linker.find("GTGA")
            if idx >= 0:
                if 43+idx > len(seq):
                    print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short (<{43+idx}) for seq: {seq}", flush=True)
                else:
                    cb1 = seq[0+idx:9+idx]
                    cb2 = seq[13+idx:22+idx]
                    cb3 = seq[26+idx:35+idx]
                    umi = seq[35+idx:43+idx]
                    cb = cb1 + cb2 + cb3                # 9+9+9=27
        else:
            print(f"get_cb_umi|R1 {sc_platform} sequence length ({len(seq)}) is too short (<43) for seq: {seq}", flush=True)
    else:
        print(f"get_cb_umi|[ERROR]: invalid sc_platform ({sc_platform}) for seq: {seq}", flush=True)
        print("Troubleshooting:", flush=True)
        print(" \u2022 The platform string decides how each read is sliced into cell barcode and UMI, so an unknown value cannot be handled", flush=True)
        print(" \u2022 It is set by the pipeline to TrekkerFX_FLEX for Flex runs; a different value means the samplesheet was hand-edited", flush=True)
        print(" \u2022 Re-run the stage through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis --force`", flush=True)
        sys.exit(1)

    return cb, umi

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
        print(f"hamming_distance|WARNING: base_seq ({base_seq}) length ({len(base_seq)}) and test_seq ({test_seq}) length ({len(test_seq)}) are not equal", flush=True)
    else:
        n_match = 0
        for i in range(len(base_seq)):
            if base_seq[i] == test_seq[i]:
                n_match += 1
        if n_match >= len(base_seq) - n:
            b_match = True
    return b_match

def validate_up(seq, b_up):
    """confirms that 'UP' sequence is present and that sequence has minimum length of 32
    Args:
        seq (string): sequence
        b_up (bool): use original UP sequence verification
    Returns:
        b_up (bool): True if UP is present, False otherwise
    """
    b_vup = False
    if len(seq) >= 32:
        if b_up:
            if "TCTTCAGCGTTCCCGAGA" in seq:
                b_vup = True
        else:
            up_seq = seq[8:26]
            if up_seq == "TCTTCAGCGTTCCCGAGA":
                b_vup = True
            elif hamming_distance(up_seq, "TCTTCAGCGTTCCCGAGA", 1):
                b_vup = True
    else:
        print(f"validate_up|ERROR: Read 2 sequence length < 32: {seq}", flush=True)
    return b_vup

def get_sb(seq):
    """get spatial barcode (sb) from sequence
    Args:
        seq (string): sequence
    Returns:
        sb (string): spatial barcode sequence
    """
    sb = ""
    sb = seq[0:8] + seq[26:32]
    return sb

def parse_whitelist(sc_platform, f_whitelist):
    """parse cell barcode translation from f_whitelist.txt file
    Args:
        sc_platform (string): type of sequencing platform
        f_whitelist (string): path/to whitelist.txt file
    Returns:
        d_CB (dict): dictionary of cell barcodes and their translation
        d_CB is returned empty if sc_platform not in ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_M', 'TrekkerFX_FLEX', 'TrekkerQ_S', 'TrekkerQ_P')
    """
    if sc_platform in ('tenx', 'TrekkerC', 'TrekkerCX'):
        print(f"parse_whitelist|Parsing {f_whitelist}...")
    else:
        print(f"parse_whitelist|No whitelist needed...")
    d_CBx = {}
    if sc_platform in ('tenx', 'TrekkerC', 'TrekkerCX'):
        with gzip.open(f_whitelist,"rt") as FILE_WHITELIST:
            for line in FILE_WHITELIST:
                sline = line.rstrip().split('\t')
                if len(sline) == 2:
                    d_CBx[sline[0]] = sline[1]
                else:
                    print(f"parse_whitelist|ERROR: Invalid line in file: {f_whitelist}, line: {line}", flush=True)
        print(f"parse_whitelist|SUMMARY: {len(d_CBx)} cell barcodes parsed from {f_whitelist}", flush=True)
        print(f"parse_whitelist|SUMMARY: memory usage: {get_memory()}", flush=True)
    return d_CBx

def parse_sc_cb(sample_name, sc_platform, f_barcodes, misc, f_whitelist=None):
    """Parse cell barcodes from f_barcodes file barcodes.tsv.gz
    Args:
        sample_name (string): sample name
        sc_platform (string): Trekker sequencing platform
        f_barcodes (string): path/to/ barcodes.tsv.gz file
        f_whitelist (string): path/to/ whitelist.txt file
    Returns:
        d_CB (dict): barcodes from barcodes.tsv.gz file
    """
    # For 'tenx', 'TrekkerC', 'TrekkerCX', parse cell barcode translation from f_whitelist.txt
    d_CBx = parse_whitelist(sc_platform, f_whitelist)
    d_CB = {}
    c_invalid_cb = 0
    c_all = 0
    print(f"parse_sc_cb|Parsing {f_barcodes}...", flush=True)
    # Get all cell barcodes from f_barcodes file AND translate CB using whitelist.txt file
    if len(d_CBx) > 0:
        if sc_platform in ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_M', 'TrekkerFX_FLEX', 'Trekker5C_CX', 'Trekker5C_C'):
            with gzip.open(f_barcodes,"rt") as FILE_SC_BC:
                for line in FILE_SC_BC:
                    if '-' in line:
                        c_all += 1
                        sline = line.rstrip().split('-')
                        if sline[0] in d_CBx:
                            cbx = d_CBx[sline[0]]
                            if cbx not in d_CB:
                                d_CB[cbx] = 'd'
                        else:
                            c_invalid_cb += 1
                    else:
                        c_invalid_cb += 1
        if sc_platform in ('TrekkerQ_S', 'TrekkerQ_P'):
            with gzip.open(f_barcodes,"rt") as FILE_SC_BC:
                for line in FILE_SC_BC:
                    c_all += 1
                    if any(char.isspace() for char in line):
                        sline = line.rstrip().split()
                        if sline[0] not in d_CB:
                            d_CB[sline[0]] = 'd'
                        else:
                            c_invalid_cb += 1
                    else:
                        c_invalid_cb += 1
    # Get all cell barcodes from f_barcodes file (R index translation)
    elif sc_platform in ('TrekkerR','TrekkerU_R','Rhapsody', 'TrekkerU_RATAC', 'TrekkerU_RVDJ'):
        Enh_linker1, Enh_linker2, B384_cell_key1, B384_cell_key2, B384_cell_key3 = get_rhapsody_keys()
        with gzip.open(f_barcodes,"rt") as FILE_SC_BC:
            for line in FILE_SC_BC:
                c_all += 1
                index = int(line.rstrip())
                seq,cb = index_to_sequence(index, Enh_linker1, Enh_linker2, B384_cell_key1, B384_cell_key2, B384_cell_key3)
                #seq_nolinker=seq[0:9]+seq[13:22]+seq[26:35]
                if len(cb) >= 27:
                    if cb not in d_CB:
                        d_CB[cb] = 'd'
                else:
                    c_invalid_cb += 1

        # Write all cell barcodes to file
        cout = 0
        file_out = sample_name + "_barcodes_perPublicFxn.tsv"
        print(f"parse_sc_cb|writing valid cell barcodes to file: {file_out}", flush=True)
        with open(os.path.join(misc, file_out),"wt") as FILE_OUT:
            for cbi in d_CB:
                cout += 1
                FILE_OUT.write(f"{cbi}\n")
        print(f"parse_sc_cb|SUMMARY|wrote {cout} valid cell barcodes written to file: {file_out}", flush=True)
    # Get all cell barcodes from file f_barcodes file (NO translation)
    else:
        with gzip.open(f_barcodes,"rt") as FILE_SC_BC:
            for line in FILE_SC_BC:
                c_all += 1
                if '-' in line:
                    sline = line.rstrip().split('-')
                    if sline[0] not in d_CB:
                        d_CB[sline[0]] = 'd'
                elif any(char.isspace() for char in line):
                    sline = line.rstrip().split()
                    if sline[0] not in d_CB:
                        d_CB[sline[0]] = 'd'
                else:
                    c_invalid_cb += 1
    
    print(f"parse_sc_cb|SUMMARY: {len(d_CB)}/{c_all} valid cell barcodes parsed from {f_barcodes}", flush=True)
    print(f"parse_sc_cb|SUMMARY: {c_invalid_cb} invalid cell barcodes removed", flush=True)
    print(f"parse_sc_cb|SUMMARY: memory usage: {get_memory()}", flush=True)
    # clean up d_CB to save memory
    del d_CBx
    gc.collect()
    return d_CB

class R1:
    def __init__(self, cb, umi):
        self.cb = cb
        self.umi = umi

def parse_r1(sample_name, fastq_R1, sc_platform, f_barcodes, misc, f_whitelist=None):
    """parse_r1 parses reads from fastq_R1 file line by line and stores read_id, CB, UMI in dictionary d_R1
    Args:
        sample_name (string): name of the sample
        fastq_R1 (string): file path to R1 fastq file
        sc_platform (string): type of trekker ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_C', 'TrekkerU_CX', 'TrekkerU_M', 'TrekkerR', 'Rhapsody', 'TrekkerU_R', 'TrekkerFX_FLEX', 'TrekkerQ_S', 'TrekkerQ_P', 'TrekkerU_IL', 'TrekkerU_PIP', 'TrekkerU_F_TEST', 'TrekkerSHA_WGA', 'TrekkerU_RATAC', 'TrekkerU_RVDJ')
        f_barcodes (string): file path to cell barcodes file (barcodes.tsv.gz)
        f_whitelist (string): file path to whitelist file (whitelist.txt)
    Returns:
        dictionary: d_R1[read_id] = R1 object (R1.CB, R1.UMI)
    """
    # Get all cell barcodes from f_barcodes file
    # if sc_platform == 'tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_M', 'TrekkerFX_FLEX', 'TrekkerQ_S', 'TrekkerQ_P' then use whitelist.txt file to translate cell barcodes
    
    d_CB = parse_sc_cb(sample_name, sc_platform, f_barcodes, misc, f_whitelist)

    # Count unique CB in R1
    s_cbr1 = set()

    # Dictionary to store R1, d_R1[read_id] = R1.CB, R1.UMI
    d_R1 = {}
    # Count invalid CB
    c_icb = 0
    # Count all records in R1 file
    c_all = 0
    # Store all read_ids in R1 file
    s_R1 = set()

    print(f"parse_r1|Parsing R1 file: {fastq_R1}...", flush=True)
    # Read fastq_R1 file line by line
    with gzip.open(fastq_R1,"rt") as FILE_R1:
        recsR1 = SeqIO.parse(FILE_R1,"fastq")
        for rec in recsR1:
            c_all += 1
            id1 = rec.id
            s_R1.add(id1)
            cb,umi = get_cb_umi(rec.seq, sc_platform)

            if sc_platform in ('TrekkerU_IL', 'TrekkerU_PIP', 'TrekkerSHA_WGA') : # Use read id for Fluent
                umi = id1

            if len(cb) > 15 and len(umi) > 7:
                # Filter for valid cell barcodes (matching CB in barcodes.tsv.gz)
                if cb in d_CB:
                    s_cbr1.add(cb)
                    # Store read_id, CB, UMI for all valid matched CB
                    if id1 not in d_R1:
                        d_R1[id1] = R1(cb, umi)
                    else:
                        print(f"parse_r1|ERROR|Duplicate read_id: {id1} in R1 file: {fastq_R1}", flush=True)
                else:
                    c_icb += 1
    print(f"parse_r1|SUMMARY|successfully parsed {len(d_R1)}/{c_all} records from R1 file: {fastq_R1}", flush=True)
    print(f"parse_r1|SUMMARY|number of invalid cell barcodes excluded: {c_icb}", flush=True)
    print(f"parse_r1|SUMMARY|memory usage: {get_memory()}", flush=True)
    
    # Capture count of unique cell barcodes in barcodes.tsv.gz
    c_cb = len(d_CB)
    # Count of unique CB with match in R1 file
    c_cbr1 = len(s_cbr1)
    # Delete s_CB set to save memory
    del d_CB
    del s_cbr1
    gc.collect()
    return s_R1, d_R1, c_cb, c_cbr1, c_all

def trekker_parser(sample_name, fastq_R1, fastq_R2, f_barcodes, sc_platform, misc, b_up, f_whitelist=None):
    """Parses R1 and R2 fastq.gz Trekker files and generates whitelist file for sample
    Args:
        sample_name (string): name of Trekker sample
        fastq_R1 (string): /path/to/R1.fastq.gz
        fastq_R2 (string): /path/to/R2.fastq.gz
        f_barcodes (string): /path/to/barcodes.tsv.gz
        sc_platform (string): type of platform ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_C', 'TrekkerU_CX', 'TrekkerU_M', 'TrekkerR', 'Rhapsody', 'TrekkerU_R', 'TrekkerFX_FLEX', 'TrekkerU_IL', 'TrekkerU_PIP', 'TrekkerU_F_TEST', 'TrekkerSHA_WGA', 'TrekkerU_RATAC', 'TrekkerU_RVDJ')
        dir_whitelist (string): /path/to/whitelist/directory
        misc (string): /path/to/misc/directory
    Output:
        df_whitelist_<sample_name>.txt (file): output file containing CB, SB, nUMI
        <sample_name>_barcodes_perPublicFxn.tsv (file): output file containing all CB in barcodes.tsv.gz (ONLY FOR R-based samples)
        reads_per_SB_<sample_name>.txt (file): output file containing unique SB (after UP filter), number of reads associated with each SB
        <sample_name>_reads_perMatchedCB.txt (file): output file containing unique CB (after UP filter), number of reads associated with each CB
        <sample_name>_reads_perMatchedCB_SB.txt (file): output file containing number of reads associated with each CB+SB, CB, SB
        matcher_summary_<sample_name>.txt (file): output file containing UP_site_matches, cb_snRNAseq, cb_matched, cb_matched_reads, Median_reads_per_nuclei, Mean_reads_per_nuclei
    """
    start_time = time.time()

    # d_umi[CB,SB] = set of UMIs
    d_umi = {}
    
    # Parse fastq_R1
    # d_R1[read_id] = R1.CB, R1.UMI
    # c_cb = number of unique CB in barcodes.tsv.gz
    # c_cbr1 = number of unique CB with match in R1 file
    s_R1,d_R1,c_cb,c_cbr1,c_all = parse_r1(sample_name, fastq_R1, sc_platform, f_barcodes, misc, f_whitelist)
    
    # Count validated UP sequences in R2
    c_up = 0
    # Count validated UP sequences in R2 with matched read_id and CB in R1
    c_upcb = 0
    # Count all R2 records
    c_r2 = 0
    # Count unique valid SB read 2
    d_sb = {}
    # Count unique valid CB read 2
    d_cb = {}
    # Count unique valid CB+SB read 2
    d_cbsb = {}
    # Count invalid UP sequence in R2
    c_invUP = 0
    # Count invalid R2
    c_invR2 = 0
    
    print(f"parse_r2|Parsing R2 file: {fastq_R2}...", flush=True) 
    with gzip.open(fastq_R2,"rt") as FILE_R2:
        recsR2 = SeqIO.parse(FILE_R2,"fastq")
        for rec in recsR2:
            c_r2 += 1
            id2 = rec.id
            # Check if the read has the UP tag and is appropriate length
            # 1. sc_platform not eq  TrekkerU_PIP   (dont care about second condition)
            # 2. sc_platform eq TrekkerU_PIP --> id2 in s_R1   (second condition has to be true to enter
            if (sc_platform not in ('TrekkerU_PIP',)) or (id2 in s_R1):
                if validate_up(rec.seq, b_up):
                    c_up += 1
                    # Validate R2 read_id matches an R1 read_id
                    if id2 in d_R1:
                        c_upcb += 1
                        sb = get_sb(rec.seq)
                        cb = d_R1[id2].cb
                        key = cb + ',' + sb
                        
                        # Count unique valid SB
                        if sb not in d_sb:
                            d_sb[sb] = 0
                        d_sb[sb] += 1

                        # Count unique valid CB
                        if cb not in d_cb:
                            d_cb[cb] = 0
                        d_cb[cb] += 1

                        # Count unique valid CB+SB
                        if key not in d_cbsb:
                            d_cbsb[key] = 0
                        d_cbsb[key] += 1

                        # Count unique valid UMI per CB+SB
                        if key not in d_umi:
                            d_umi[key] = set()

                        umi = ""
                        if umi == "":
                            umi = d_R1[id2].umi
                        d_umi[key].add(umi)
                    else:
                        c_invR2 += 1
                else:
                    c_invUP += 1
   
    print(f"trekker_parser|SUMMARY|number of invalid R2 tags (no R1 match for read_id): {c_invR2}")
    print(f"trekker_parser|SUMMARY|number of invalid UP records: {c_invUP}")
    print(f"trekker_parser|SUMMARY|UP validated {c_up}/{c_r2} records from R2 file: {fastq_R2}")
    print(f"trekker_parser|SUMMARY|successfully parsed {len(d_umi)}/{c_r2} records from R2 file: {fastq_R2}")
    print(f"trekker_parser|SUMMARY|memory usage: {get_memory()}")

    # save memory by deleting the entire R1 dictionary
    del d_R1
    gc.collect()

    # save memory by deleting the entire R1 set created for Fluent
    del s_R1
    gc.collect()

    # Write output file, number of read 2 per valid SB
    cout = 0
    print(f"trekker_parser|writing {len(d_cbsb)} records to output file: reads_per_SB_{sample_name}.txt", flush=True)
    with open(os.path.join(misc, f'reads_per_SB_{sample_name}.txt'), 'w') as FILE_OUT:
        for sb, count in d_sb.items():
            cout += 1
            FILE_OUT.write(f'{sb},{count}\n')
    print(f"trekker_parser|SUMMARY|Successfully wrote {cout} records to file: reads_per_SB_{sample_name}.txt", flush=True)
    # save memory by deleting the d_sb dictionary
    del d_sb
    gc.collect()

    # Write output file, number of read 2 per valid CB
    cout = 0
    l_cb = []
    print(f"trekker_parser|writing {len(d_cb)} records to output file: {sample_name}_reads_perMatchedCB.txt", flush=True)
    with open(os.path.join(misc, f'{sample_name}_reads_perMatchedCB.txt'), 'w') as FILE_OUT:
        for cb, count in d_cb.items():
            cout += 1
            FILE_OUT.write(f'{cb},{count}\n')
            l_cb.append(count)
    print(f"trekker_parser|SUMMARY|Successfully wrote {cout} records to file: {sample_name}_reads_perMatchedCB.txt", flush=True)
    # save memory by deleting the d_cb dictionary
    del d_cb
    gc.collect()

    # Write output file, <sample_name>_reads_perMatchedCB_SB.txt
    # Columns: nReads (CB + SB), CB, SB
    cout = 0
    #l_cbsb = []
    print(f"trekker_parser|writing {len(d_cbsb)} records to output file: {sample_name}_reads_perMatchedCB_SB.txt", flush=True)
    with open(os.path.join(misc, f'{sample_name}_reads_perMatchedCB_SB.txt'), 'w') as FILE_OUT:
        for key, n_cbsb in d_cbsb.items():
            cb,sb = key.split(',')
            cout += 1
            FILE_OUT.write(f'{n_cbsb},{cb},{sb}\n')
    print(f"trekker_parser|SUMMARY|Successfully wrote {cout} records to file: {sample_name}_reads_perMatchedCB_SB.txt", flush=True)
    # save memory by deleting the d_cbsb dictionary
    del d_cbsb
    gc.collect()

    # Write output file, df_whitelist_sample_name.txt
    cout = 0
    #l_umi = []
    print(f"trekker_parser|writing {len(d_umi)} records to output file: df_whitelist_{sample_name}.txt", flush=True)
    with open(os.path.join(misc, f'df_whitelist_{sample_name}.txt'), 'w') as FILE_OUT:
        FILE_OUT.write('CB,SB,nUMI\n')
        for key, umi_set in d_umi.items():
            cout += 1
            FILE_OUT.write(f'{key},{len(umi_set)}\n')
            #l_umi.append(len(umi_set))
    print(f"trekker_parser|SUMMARY|successfully wrote {cout} records to output file: df_whitelist_{sample_name}.txt", flush=True)
    # clean up
    del d_umi
    gc.collect()

    # Write output file, matcher_summary_<sample_name>.txt
    print(f"trekker_parser|writing file: matcher_summary_{sample_name}.txt", flush=True)
    if l_cb:
        median_bc_reads_per_nucleus = statistics.median(l_cb)
        mean_bc_reads_per_nucleus = statistics.mean(l_cb)
    else:
        print(f"trekker_parser|WARNING|no matched cell barcodes for sample: {sample_name}, reporting 0 for median/mean reads per nucleus", flush=True)
        median_bc_reads_per_nucleus = 0
        mean_bc_reads_per_nucleus = 0
    with open(os.path.join(misc, f'matcher_summary_{sample_name}.txt'), 'w') as FILE_OUT:
        FILE_OUT.write("Total_reads,UP_site_matches,cb_snRNAseq,cb_matched,cb_matched_reads,Median_reads_per_nuclei,Mean_reads_per_nuclei\n")
        FILE_OUT.write(f"{c_all},{c_up},{c_cb},{c_cbr1},{c_upcb},{median_bc_reads_per_nucleus},{mean_bc_reads_per_nucleus}\n")
    print(f"trekker_parser|SUMMARY|Successfully wrote file: matcher_summary_{sample_name}.txt", flush=True)
    # clean up
    del l_cb
    gc.collect()

    # Log memory and time elapsed
    print(f"trekker_parser|SUMMARY|memory usage: {get_memory()}", flush=True)
    time_elapsed = time.time() - start_time
    if time_elapsed < 60:
        time_elapsed = f"{time_elapsed:.2f} seconds"
    elif time_elapsed >= 60:
        time_elapsed = f"{time_elapsed/60:.2f} minutes"
    print(f"trekker_parser|SUMMARY|time elapsed: {time_elapsed}", flush=True)

def get_whitelist(dir_whitelist, sc_platform, sc_outdir):
    """generate path to whitelist file based on platform type
    Args:
        dir_whitelist (string): path to whitelist directory
        sc_platform (string): type of Trekker tile ('tenx', 'TrekkerC', 'TrekkerCX', 'TrekkerU_C', 'TrekkerU_CX', 'TrekkerU_M', 'TrekkerFX_FLEX', 'TrekkerR', 'Rhapsody', 'TrekkerU_R', 'TrekkerQ_S', 'TrekkerQ_P', 'TrekkerU_PIP', 'TrekkerU_F_TEST', 'TrekkerU_RATAC', 'TrekkerU_RVDJ')
    Returns:
        f_whitelist (string): path to whitelist file
    """
    f_whitelist = ''
    if sc_platform in ('tenx', 'TrekkerC'):
        f_whitelist = os.path.join(dir_whitelist, "Aug_2023_TrekkerC_3M-february-2018.txt.gz")
    elif sc_platform in ('TrekkerCX',):
        f_whitelist = os.path.join(dir_whitelist, "Apr_2024_TrekkerCX_3M-3pgex-may-2023.txt.gz")
    else:
        print(f"get_whitelist|[ERROR]: invalid platform type: {sc_platform}", flush=True)
        print("Troubleshooting:", flush=True)
        print(" \u2022 The platform string selects which barcode whitelist to load, so an unknown value cannot be handled", flush=True)
        print(" \u2022 It is set by the pipeline to TrekkerFX_FLEX for Flex runs; a different value means the samplesheet was hand-edited", flush=True)
        print(" \u2022 Re-run the stage through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis --force`", flush=True)
        sys.exit(1)
    
    if not os.path.isfile(f_whitelist):
        print(f"get_whitelist|[ERROR]: barcode whitelist file not found: {f_whitelist}", flush=True)
        print("Troubleshooting:", flush=True)
        print(" \u2022 This whitelist ships with the Trekker profiling scripts; a missing file means the checkout is incomplete", flush=True)
        print(f" \u2022 Restore it with `git checkout {os.path.basename(f_whitelist)}` or `git pull`", flush=True)
        print(f" \u2022 Expected under: {dir_whitelist}", flush=True)
        sys.exit(1)
    
    return f_whitelist

def main():
    """Main function of fastq_parser.py, parses samplesheet.csv for sample_name, and paths to fastq_R1, fastq_R2
    CLI Args:
        sample_name:
        fastq_R1:
        fastq_R2:
        sc_outdir:
        dir_whitelist:
        sc_platform:
        misc:
    """
    b_up = False
    if len(sys.argv) == 8:
        sample_name=sys.argv[1]
        fastq_R1=sys.argv[2]
        fastq_R2=sys.argv[3]
        sc_outdir=sys.argv[4]
        dir_whitelist=sys.argv[5]
        sc_platform=sys.argv[6]
        misc=sys.argv[7]
    elif len(sys.argv) == 9:
        sample_name=sys.argv[1]
        fastq_R1=sys.argv[2]
        fastq_R2=sys.argv[3]
        sc_outdir=sys.argv[4]
        dir_whitelist=sys.argv[5]
        sc_platform=sys.argv[6]
        misc=sys.argv[7]
        if sys.argv[8] == "--up":
            b_up = True
    else:
        print(f"fastq_parser|[ERROR]: incorrect number of command line arguments ({len(sys.argv)-1})", flush=True)
        print("Troubleshooting:", flush=True)
        print(" \u2022 This script is normally invoked by the Flex spatial-profiling stage, not by hand", flush=True)
        print(" \u2022 Run it through the pipeline with `./slidr --bcl <BCL_ID> --spatial-analysis`", flush=True)
        print(f" \u2022 Received: {sys.argv[1:]}", flush=True)
        sys.exit(1)
    
    if sc_platform in ('TrekkerQ_P', 'TrekkerQ_S'):
       fastq_R1, fastq_R2 = swap_fastqs(fastq_R1, fastq_R2)

    # Create full path to whitelist file
    if sc_platform in ('tenx', 'TrekkerC', 'TrekkerCX'):
       f_whitelist = get_whitelist(dir_whitelist, sc_platform, sc_outdir)
       print(f_whitelist)
    else:
       f_whitelist = None

    # Create full path to barcodes.tsv.gz file
    if sc_platform in ('TrekkerFX_FLEX',):
        f_barcodes_input = os.path.join(sc_outdir+"/barcodes.tsv.gz")
        f_barcodes = os.path.join(misc+"/"+sample_name+"_barcodes_trekkerFX.tsv.gz")
        cleanup_barcodes_TrekkerFX(f_barcodes_input, f_barcodes)
    elif sc_platform in ('TrekkerU_IL',):
        f_barcodes = os.path.join(sc_outdir+"/trekkerinterim/barcodes.tsv.gz")
        print(f_barcodes)
    else:
        f_barcodes = os.path.join(sc_outdir+"/barcodes.tsv.gz")

    # Execute trekker_parser()
    trekker_parser(sample_name, fastq_R1, fastq_R2, f_barcodes, sc_platform, misc, b_up, f_whitelist)

if __name__ == "__main__":
    main()

