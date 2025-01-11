#!/usr/bin/env python
"""
csv2pickle.py
^^^^^^^^^^^^^
This script reads a CSV file and saves it as a pickle file. It supports optional gzip compression.

Usage:
    python csv2pickle.py -f <filename.csv> --gzip

Arguments:

    -f FILE, --csv_file FILE
        Path to the input CSV file.

    --gzip
        Enable gzip compression for the output pickle file.

Notes:
    The script reads the CSV file from the current working directory and saves the pickle file in the same directory.

"""

import os, argparse
import pandas as pd

usage = "Read CSV file and save as a pickle. Usage: python csv2pickle.py -f <filename.csv> --gzip"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="csv_file", help="SAM file",
                 metavar="FILE", default=None)
files.add_argument("--gzip", dest="compression", help="gzip compression",action="store_true")
args = parser.parse_args()

path = os.getcwd()+"/"
filename=args.csv_file

df = pd.read_csv(path+filename, sep=",",index_col=0,header=0)
if args.compression == True:
    df.to_pickle(path + filename.replace(".csv", ".pcl.gz"),compression='gzip')
else:
    df.to_pickle(path+filename.replace(".csv",".pcl"))