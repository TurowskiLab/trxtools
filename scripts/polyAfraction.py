#!/usr/bin/env python
"""
polyAfraction.py
^^^^^^^^^^^^^^^^^^^^^^

This script generates fraction (eg. polyA) for genomic profiles from BigWig files. Bot fwd and rev files will be processed.

Usage:
    polyAfraction.py -f file_polyA_fwd.bw

Options:
    -f FILE, --sam_file FILE
        BigWig file to be processed. Requires file with suffix "_fwd.bw" or "_rev.bw". Either can be used.
        The script will automatically detect the suffix and process the file accordingly.
        The script will also look for the corresponding file with the other suffix (if not provided).
        For example, if you provide "file_polyA_fwd.bw", it will look for "file_polyA_rev.bw" in the same directory.
    -s SUFFIX, --suffix SUFFIX
        Suffix for the input file (default: polyA).
    -h, --help
        Show this help message and exit.

Example:
    polyAfraction.py -f file_polyA_fwd.bw

"""

import os, argparse
import trxtools as tt
import pyBigWig
import pandas as pd
import numpy as np

usage = "his script generates fraction (eg. polyA) for genomic profiles from BigWig files. To  \n" 
usage += "Example usage: polyAfraction.py -f file_polyA_fwd.bw"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="bw_file", help="BigWig file to be processed. Requires file with suffix '_fwd.bw' or '_rev.bw'. Either can be used. The script will automatically detect the suffix and process the file accordingly. The script will also look for the corresponding file with the other suffix (if not provided). For example, if you provide 'file_polyA_fwd.bw', it will look for 'file_polyA_rev.bw' in the same directory.",
                 metavar="FILE", default=None)
files.add_argument("-s", dest="suffix", help="Suffix for the input file (default: polyA)", default="polyA",)
args = parser.parse_args()

### main code ###
path = os.getcwd()+"/"
filename = args.bw_file.replace("_fwd.bw","").replace("_rev.bw","")
name = filename.replace("_"+args.suffix,"")

#### FWD ####
print("Processing forward strand BigWig file: ", args.bw_file)
bw_A = pyBigWig.open(path+name+"_"+args.suffix+"_fwd.bw")
print("Corresponding strand BigWig file: ", name+"_rev.bw")
bw_B = pyBigWig.open(path+name+"_fwd.bw")

chroms = bw_A.chroms()
chroms_touple = [i for i in chroms.items()] 
chroms_touple.sort(key=lambda x: int(x[1]), reverse=True)  # Sort by chromosome len - longest first
# print(chroms_touple)
# chroms_touple = [c for c in chroms_touple if c[0] in ["chrM"]]  # Include listed chromosomes only, e.g. ["chrM"] for mitochondrial DNA

paths = {}
for c in chroms_touple:
    chr, len_chr = c[0], c[1]
    print("Processing chromosome:", chr, "with length:", len_chr)
    A_series = np.array(bw_A.values(chr,0,len_chr))
    if np.isnan(A_series).all():
        print("Warning: No data found for chromosome", chr)
        output_fraction = np.full(len_chr, np.nan)
    else:
        raw = np.array(bw_B.values(chr,0,len_chr))
        output_fraction = A_series/raw
    
    df_to_save = pd.DataFrame(output_fraction, columns=[chr], index=np.arange(len_chr)+1)  # Index starts from 1
    temp_filename = path + "temp_" + name + "_fractionOF" + args.suffix + "_fwd_" + chr + ".pkl.gz"
    df_to_save.to_pickle(temp_filename, compression="gzip")
    paths[chr] = temp_filename

print("Saving BigWig file with fraction of "+args.suffix+"...")
bw_name = path + name + "_fractionOF"+args.suffix+"_fwd.bw"
l = tt.SAMgeneral.saveBigWig(paths=paths,suffix="",bw_name=bw_name,chroms=chroms_touple,pkl=True)
print("BigWig file saved as: ", bw_name)

# Remove all temporary files created and stored in paths
for temp_file in paths.values():
    if os.path.exists(temp_file):
        os.remove(temp_file)

# #### REV ####

print("Processing reverse strand BigWig file: "+name+"_rev.bw")
bw_A = pyBigWig.open(path+name+"_"+args.suffix+"_rev.bw")
bw_B = pyBigWig.open(path+name+"_rev.bw")

chroms = bw_A.chroms()
chroms_touple = [i for i in chroms.items()] 
chroms_touple.sort(key=lambda x: int(x[1]), reverse=True)  # Sort by chromosome len - longest first
# print(chroms_touple)
# chroms_touple = [c for c in chroms_touple if c[0] in ["chrM"]]  # Include listed chromosomes only, e.g. ["chrM"] for mitochondrial DNA

paths = {}
for c in chroms_touple:
    chr, len_chr = c[0], c[1]
    print("Processing chromosome:", chr, "with length:", len_chr)
    A_series = np.array(bw_A.values(chr,0,len_chr))
    if np.isnan(A_series).all():
        print("Warning: No data found for chromosome", chr)
        output_fraction = np.full(len_chr, np.nan)
    else:
        raw = np.array(bw_B.values(chr,0,len_chr))
        output_fraction = A_series/raw
    
    df_to_save = pd.DataFrame(output_fraction, columns=[chr], index=np.arange(len_chr)+1)  # Index starts from 1
    temp_filename = path + "temp_" + name + "_fractionOF" + args.suffix + "_fwd_" + chr + ".pkl.gz"
    df_to_save.to_pickle(temp_filename, compression="gzip")
    paths[chr] = temp_filename

print("Saving BigWig file with fraction of "+args.suffix+"...")
bw_name = path + name + "_fractionOF"+args.suffix+"_rev.bw"
l = tt.SAMgeneral.saveBigWig(paths=paths,suffix="",bw_name=bw_name,chroms=chroms_touple,pkl=True)
print("BigWig file saved as: ", bw_name)

# Remove all temporary files created and stored in paths
for temp_file in paths.values():
    if os.path.exists(temp_file):
        os.remove(temp_file)