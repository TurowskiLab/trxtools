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
bw_A = pyBigWig.open(path+name+"_"+args.suffix+"_fwd.bw")
bw_B = pyBigWig.open(path+name+"_fwd.bw")

chroms = bw_A.chroms()
chroms_touple = [i for i in chroms.items()]

paths = {}
for chr in chroms:
    len_chr = int(chroms[chr])
    A_series = np.array(bw_A.values(chr,0,len_chr))
    raw = np.array(bw_B.values(chr,0,len_chr))
    fraction_polyA = A_series/raw
    paths[chr] = fraction_polyA

bw_name = path + name + "_fractionOF"+args.suffix+"_fwd.bw"
l = tt.SAMgeneral.saveBigWig(paths=paths,suffix="",bw_name=bw_name,chroms=chroms_touple,pkl=False)
print("BigWig file saved as: ", bw_name)

#### REV ####
bw_A = pyBigWig.open(path+name+"_"+args.suffix+"_rev.bw")
bw_B = pyBigWig.open(path+name+"_rev.bw")

chroms = bw_A.chroms()
paths = {}
for chr in chroms:
    len_chr = int(chroms[chr])
    A_series = np.array(bw_A.values(chr,0,len_chr))
    raw = np.array(bw_B.values(chr,0,len_chr))
    fraction_polyA = A_series/raw
    paths[chr] = fraction_polyA

bw_name = path + name + "_fractionOF"+args.suffix+"_rev.bw"
l = tt.SAMgeneral.saveBigWig(paths=paths,suffix="",bw_name=bw_name,chroms=chroms_touple,pkl=False)
print("BigWig file saved as: ", bw_name)