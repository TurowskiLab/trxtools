#!/usr/bin/env python
from TTools.SAMgenome import sam2genome, sam2genome3end, sam2genome5end
import os, argparse
import TTools as tt

usage = "Generates genomic profiles from SAM files."

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="sam_file", help="SAM file",
                 metavar="FILE", default=None)
files.add_argument("-u", dest="use", help="Generate profiles using reads (read) or the 3' ends (end)",choices=["read","3end",'5end'], default="3end" )
# files.add_argument("--del", dest="deletions", help="Generate additional profiles for deletions",action="store_true")
# files.add_argument("-e", dest="expand", help="For deletions position can be expanded by the value on each side (e=5 gives 10 nt long)", type=int, default=5)
files.add_argument("-n", dest="noncoded_pA", help="Save non-coded polyA ends. Can be used ONLY with: -u 3end",
                        action="store_true",default=False)
files.add_argument("-r", dest="noncoded_raw", help="Save all non-coded ends. Can be used ONLY with: -u 3end",
                        action="store_true",default=False)
files.add_argument("--csv", dest="csv", help="Save output as csv. Default is pickle",action="store_true",default=False)
files.add_argument("-c", dest="toClear", help="String of signs to be cleared from the name of SAM file", type=str, default='_comp_flexbar_STARAligned.out')
files.add_argument("--chunks", dest="chunks", help="Divide list of genes into a chunks of given size", type=int, default=0)
args = parser.parse_args()

print("Start at "+tt.methods.timestamp())

path = os.getcwd()+"/"

if args.csv==False:
    pickle = True
elif args.csv==True:
    pickle = False

if args.use=="read":
    sam2genome(filename=args.sam_file, path=path, toClear=args.toClear,
                pickle=pickle, chunks=args.chunks)
elif args.use=="3end":
    sam2genome3end(filename=args.sam_file, path=path, toClear=args.toClear,
                pickle=pickle, chunks=args.chunks, noncoded_pA=args.noncoded_pA, noncoded_raw=args.noncoded_raw)
elif args.use=="5end":
    sam2genome5end(filename=args.sam_file, path=path, toClear=args.toClear,
                pickle=pickle, chunks=args.chunks)

print("Done at "+tt.methods.timestamp())