"""
matching bead barcode 
"""
import os
from collections import Counter
import editdistance
from multiprocessing import Pool
import logging
import argparse


# The bead barcode/position files are generated from the puck whitelist by splitspatialbarcodes.py,
# so every malformed-input failure below points at the same few causes and offers the same fixes.
PUCK_HINTS = (
    "Troubleshooting:\n"
    " • BeadBarcodes.txt and BeadLocations.txt are generated from the puck whitelist by splitspatialbarcodes.py\n"
    " • BeadLocations.txt must be exactly two comma-separated lines: all x values, then all y values\n"
    " • BeadBarcodes.txt must hold one barcode per line, one per coordinate pair\n"
    " • Check the `Puck ID` metadata column names the puck this sample was actually run on\n"
    " • Delete the cached puck under the run's flex/pucks/ directory to force a fresh download, then re-run --spatial-analysis --force"
)


def get_barcode_position(bead_bc_file,bead_pos_file):
    def low_cpx_bc(bc,max_homo):
        return Counter(list(bc)).most_common(1)[0][1]>max_homo
    tmp = open(bead_pos_file).readlines()
    if len(tmp) < 2:
        raise ValueError(
            f"{bead_pos_file} must contain at least 2 lines (x-coords, y-coords), found {len(tmp)}\n"
            f"{PUCK_HINTS}"
        )
    # the line-count guard above doesn't catch a present-but-malformed line (blank, truncated, or
    # containing a non-numeric token) — parse explicitly so those raise a clear error naming the
    # file and which line, instead of a bare "could not convert string to float" from deep inside
    # the list comprehension
    def _parse_coord_line(line, which):
        try:
            return [float(it) for it in line.strip().split(",")]
        except ValueError as exc:
            raise ValueError(
                f"{bead_pos_file} has a malformed {which}-coordinate line: {exc}\n"
                f"{PUCK_HINTS}"
            )
    pos_x = _parse_coord_line(tmp[0], "x")
    pos_y = _parse_coord_line(tmp[1], "y")
    bc_list = [it.strip().replace(",","") for it in open(bead_bc_file).readlines()]
    if not bc_list:
        raise ValueError(f"No bead barcodes found in {bead_bc_file}\n{PUCK_HINTS}")
    if not (len(bc_list) == len(pos_x) == len(pos_y)):
        raise ValueError(
            f"Bead barcode/position count mismatch: {len(bc_list)} barcodes, "
            f"{len(pos_x)} x-coords, {len(pos_y)} y-coords in {bead_bc_file}/{bead_pos_file}\n"
            f"{PUCK_HINTS}"
        )
    max_homo = int(0.8*len(list(bc_list[0])))  # assume all have the same length
    bc_pos_dict = {}
    for bc,coordx,coordy in zip(bc_list,pos_x,pos_y):
        if not low_cpx_bc(bc,max_homo):
            bc_pos_dict[bc] = (coordx,coordy)
    return bc_pos_dict, bc_list


def build_6mer_dist(bc_list):
    start_km = {}
    mid_km = {}
    end_km = {}
    for bc in bc_list:
        start_km.setdefault(bc[:6] , []).append(bc)
        mid_km.setdefault(bc[4:10], []).append(bc)
        end_km.setdefault(bc[-6:] , []).append(bc)
    return start_km,mid_km,end_km


def barcode_matching(bc_pos_dict,spatial_bc_list,max_dist=1):
    bc_matching_dict = {}
    def get_sel_bc(bc):
        res = []
        if bc[:6] in start_km:
            res += start_km[bc[:6]]
        if bc[-6:] in end_km:
            res += end_km[bc[-6:]]
        if bc[4:10] in mid_km:
            res += mid_km[bc[4:10]]
        return set(res)
    exact_match = 0
    fuzzy_match =0
    bc_ref_list = list(bc_pos_dict.keys())
    start_km,mid_km,end_km = build_6mer_dist(bc_ref_list)
    for bc in spatial_bc_list:
        if bc in bc_pos_dict:
            exact_match += 1
            bc_matching_dict[bc] = bc
        else:
            sel_bc = get_sel_bc(bc)
            if len(sel_bc)>0:
                fz = [(it, editdistance.eval(it, bc)) for it in sel_bc]
                fz = [it for it in fz if it[1]<=max_dist]
                fz.sort(key=lambda x:x[1])
                if len(fz)==0:
                    continue
                if len(fz)>1 and fz[0][1]==fz[1][1]:
                    last_base_dist0 = editdistance.eval(fz[0][0][-1], bc[-1])
                    last_base_dist1 = editdistance.eval(fz[1][0][-1], bc[-1])
                    fuzzy_match += 1
                    if last_base_dist0 > last_base_dist1:  # higher error rate in the last base of the barcode
                        bc_matching_dict[bc] = fz[1][0]
                    else:
                        # last_base_dist0 < last_base_dist1, or a genuine tie: keep the first
                        # candidate (already the lowest overall edit distance) instead of dropping it
                        bc_matching_dict[bc] = fz[0][0]
                else:
                    fuzzy_match += 1
                    bc_matching_dict[bc] = fz[0][0]
    return bc_matching_dict,exact_match,fuzzy_match

def split_list(li, trunk_size=3000):
    new_li = []
    new_st = 0
    if len(li)<trunk_size:
        return [li]
    while True:
        if (new_st+trunk_size)<len(li):
            new_li.append(li[new_st:(new_st+trunk_size)])
            new_st += trunk_size
        else:
            new_li.append(li[new_st:] )
            break
    return new_li


def write_barcode_match(bc_matching_dict,bc_pos_dict,output_file):
    """_summary_

    Args:
        bc_matching_dict (dict): key:Illumina barcode, value: matched barcode
        bc_pos_dict (dict): key:puck barcode, value: (x, y)
        output_file (string): absolute path of output file
    """
    with open(output_file,"w") as fo:
        fo.write("Illumina_barcode,matched_beadbarcode,xcoord,ycoord\n")
        for bc in bc_matching_dict:
            fo.write("{},{},{},{}\n".format(bc,
                                          bc_matching_dict[bc],
                                          bc_pos_dict[bc_matching_dict[bc]][0],
                                          bc_pos_dict[bc_matching_dict[bc]][1]))


def main(puck_id,barcode_file,bead_barcode_dir,output_file):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    bead_bc_file = os.path.join(bead_barcode_dir,puck_id,"BeadBarcodes.txt")
    bead_pos_file = os.path.join(bead_barcode_dir,puck_id,"BeadLocations.txt")
    bc_pos_dict,_ = get_barcode_position(bead_bc_file,bead_pos_file)
    logging.info("Read {} bead barcodes from puckid folder.".format(len(bc_pos_dict)) )
    #spatial_bc_list = open(barcode_file).readline().strip().split(",")
    spatial_bc_list = [line.strip().split(",")[0] for line in open(barcode_file)]
    print(spatial_bc_list[:10])
    logging.info("Number of barcode from Illumina reads: {}".format(len(spatial_bc_list)))
    logging.info("First 5 barcode: {}".format(spatial_bc_list[:5]))
    logging.info("Start matching...")
    spatial_bc_list_split = split_list(spatial_bc_list)
    with Pool(8) as p:
        res_list = p.starmap(barcode_matching, [(bc_pos_dict, it) for it in spatial_bc_list_split ] )  # match bead barcode from Illumina reads to bead barcode from puck sequencing
        bc_matching_dict = {}
        for d,_,_ in res_list:
            for k in d:
                bc_matching_dict[k] = d[k]
        exact_match = sum(it[1] for it in res_list)
        fuzzy_match = sum(it[2] for it in res_list)
    logging.info("Matching finished...")
    logging.info("Number of exact match: {}".format(exact_match))
    logging.info("Number of fuzzy match: {}".format(fuzzy_match))
    logging.info("Writing output to file: {}".format(output_file))
    write_barcode_match(bc_matching_dict,bc_pos_dict,output_file)
    logging.info("Done!")

def get_args():
    parser = argparse.ArgumentParser(description='Match bead barcode acquired from Illumina reads to the Trekker tile spatial barcode whitelist.')
    parser.add_argument(
        "-b", "--barcode_file",
        help="the file that contain barcodes from Illumina reads, the first line is barcode separated by comma.",
        type=str,
        required=True
        )
    parser.add_argument(
        "-o", "--outcsv",
        help="output csv file .",
        type=str,
        required=True
        )
    parser.add_argument(
        "-i", "--puckid",
        help="the Tile ID such as U0001_001.",
        type=str,
        required=True
        )
    parser.add_argument(
        "-d", "--bead_barcode_dir",
        help="the Tile ID such as U0001_001. will look for the folder in misc, e.g. misc/U0001_001",
        type=str,
        required=True
        )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    main(args.puckid,args.barcode_file,args.bead_barcode_dir,args.outcsv)
    # python bead_matching.py -b rc_dge.csv -i U0001_001 -o matching_result.csv

