#!/usr/bin/env python
from TTools.SAMgenome import sam2genome, sam2genome3end
import os, argparse
import pandas as pd
import TTools as tt

usage = "Generates genomic profiles from SAM files."

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="sam_file", help="SAM file",
                 metavar="FILE", default=None)
files.add_argument("-u", dest="use", help="Generate profiles using reads (read) or the 3' ends (end)",choices=["read","end"], default="end" )
# files.add_argument("--del", dest="deletions", help="Generate additional profiles for deletions",action="store_true")
# files.add_argument("-e", dest="expand", help="For deletions position can be expanded by the value on each side (e=5 gives 10 nt long)", type=int, default=5)
files.add_argument("-n", dest="noncoded", help="Process and save non-coded ends. Can be used ONLY with: -u end",
                        action="store_true",default=False)
files.add_argument("-p", dest="pickle", help="Save output as pickle DataFrame",action="store_true",default=False)
files.add_argument("-c", dest="toClear", help="String of signs to be cleared from the name of SAM file", type=str, default='_comp_flexbar_STARAligned.out')
files.add_argument("--chunks", dest="chunks", help="Divide list of genes into a chunks of given size", type=int, default=0)
args = parser.parse_args()

print("Start at "+tt.methods.timestamp())

path = os.getcwd()+"/"

if args.use=="read":
    sam2genome(filename=args.sam_file, path=path, toClear=args.toClear,
                pickle=args.pickle, chunks=args.chunks)
elif args.use=="end":
    sam2genome3end(filename=args.sam_file, path=path, toClear=args.toClear,
                pickle=args.pickle, chunks=args.chunks, noncoded=args.noncoded)

print("Done at "+tt.methods.timestamp())