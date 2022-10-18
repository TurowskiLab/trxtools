#!/usr/bin/env python

import trxtools as tt
import os, argparse
import pandas as pd

usage = "merges Salmon quantification files to one tab file"

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