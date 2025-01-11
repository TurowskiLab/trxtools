#!/usr/bin/env python
"""
mergeSalmon.py
^^^^^^^^^^^^^^

This script merges multiple Salmon quantification files into a single tab-delimited file. It takes the name of the Salmon output directory and the column to be used (TPM or NumReads). The script can also filter the files to be loaded, clear some strings from the experiment name, and add some strings to the experiment name. The output is saved in the current directory with the name 'merged.tab'.

Usage:
    python mergeSalmon.py -i 'Salmon_output' -u TPM -f 'sample' -a 'RNAseq'

Options:
    -i
      String of signs to be found in Salmon output directory

    -o
      Name of output file (default: 'merged')

    -u 
      Column to be used {TPM, NumReads} (default: 'NumReads')

    -f 
      String of signs to be found in Salmon output directory. Optional as additional filter (default: None)

    -c  
      String of signs to be found in Salmon output directory, will be cleared (default: None)

    -a 
      String of signs to be added to experiment name (default: None)

"""
import trxtools as tt
import os, argparse
import pandas as pd

usage = "Script merges multiple Salmon quantification files to one tab file \n"
usage += "The script takes the name of the Salmon output directory and the column to be used (TPM or NumReads) \n"
usage += "The script can also filter the files to be loaded, clear some strings from the experiment name, and add some strings to the experiment name \n"
usage += "The script saves the output in the current directory with the name 'merged.tab'"
usage += "Example: python mergeSalmon.py -i 'Salmon_output' -u TPM -f 'sample' -a 'RNAseq'"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-i", dest="input",type=str, help="String of signs to be found in Salmon output directory")

files.add_argument("-o", dest="output",type=str, help="Name of output file",default="merged")
files.add_argument("-u", dest="use",type=str, help="column to be used {TPM,NumReads}",default='NumReads')

files.add_argument("-f", dest="filter",type=str, help="String of signs to be found in Salmon output directory. Optional as additional filter", default=None)
files.add_argument("-c", dest="clear",type=str, help="String of signs to be found in Salmon output directory, will be cleared",default=None)
files.add_argument("-a", dest="add",type=str, help="String of signs to be added to experiment name",default=None)

args = parser.parse_args()

print("Start at "+tt.methods.timestamp())

path = os.getcwd()+"/"
if args.filter==None: args.filter=args.input
if args.clear==None: args.clear=""
if args.add==None: args.add=""

df = tt.methods.readSalmon(nameElem=args.input, path=path, toLoad=args.filter, toClear=[args.clear], toAdd=args.add, column=args.use)
df.to_csv(path+args.output+".tab",sep="\t")

print("Done at "+tt.methods.timestamp())