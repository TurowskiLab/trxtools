#!/usr/bin/env python
"""
SAM2profilesGenomic.py
^^^^^^^^^^^^^^^^^^^^^^

This script generates genomic profiles from SAM files. It provides various options to customize the generation of profiles based on the input SAM file.

Usage:
    SAM2profilesGenomic.py -f file.sam -u 3end -n -s polyA

Options:
    -f FILE, --sam_file FILE
        SAM file to be processed.

    -u {read,3end,5end,del}, --use {read,3end,5end,del}
        Generate profiles using reads (read) or the 3' ends (3end), 5' ends (5end), or deletions (del).

    -n, --noncoded
        Save non-coded ends. Can be used ONLY with: -u 3end.

    -s {polyA}, --noncoded_suffix {polyA}
        Select non-coded ends. Can be used ONLY with: -u 3end.

    -e EXPAND, --expand EXPAND
        Will expand position of each deletion by +/- the value. Works ONLY with: -u del.

    -c TOCLEAR, --toClear TOCLEAR
        String of signs to be cleared from the name of SAM file.

    --chunks CHUNKS
        Divide list of genes into chunks of given size.

Example:
    SAM2profilesGenomic.py -f file.sam -u 3end -n -s polyA

"""

from trxtools.sam.SAMgenome import sam2genome
import os, argparse
import trxtools as tt

usage = "Generates genomic profiles from SAM files. To  \n" 
usage += "For instrictions and additional options use SAM2profilesGenomic.py --help \n" 
usage += "Example usage: SAM2profilesGenomic.py -f file.sam -u 3end -n -s polyA"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="sam_file", help="SAM file",
                 metavar="FILE", default=None)
files.add_argument("-u", dest="use", help="Generate profiles using reads (read) or the 3' ends (end)",
                        type=str,choices=["read","3end","5end","del"], default="3end")
# files.add_argument("--del", dest="deletions", help="Generate additional profiles for deletions",action="store_true")
# files.add_argument("-e", dest="expand", help="For deletions position can be expanded by the value on each side (e=5 gives 10 nt long)", type=int, default=5)
files.add_argument("-n", dest="noncoded", help="Save non-coded ends. Can be used ONLY with: -u 3end",
                        action="store_true",default=False)
files.add_argument("-s", dest="noncoded_suffix", help="Select non-coded ends. Can be used ONLY with: -u 3end",
                        type=str,choices=["polyA"],default="polyA")
files.add_argument("-e", dest="expand", help="Will expand position of each deletion by +/- the value. Works ONLY with: -u del",
                        type=int,default=0)
files.add_argument("-c", dest="toClear", help="String of signs to be cleared from the name of SAM file", type=str, default='_comp_flexbar_STARAligned.out')
files.add_argument("--chunks", dest="chunks", help="Divide list of genes into a chunks of given size", type=int, default=0)
args = parser.parse_args()

# print("Start at "+tt.methods.timestamp())

path = os.getcwd()+"/"

sam2genome(filename=args.sam_file,
                path=path, 
                toClear=args.toClear,
                chunks=args.chunks, 
                use=args.use,
                noncoded=args.noncoded, # for 3end only
                noncoded_suffix=args.noncoded_suffix, # for 3end and noncoded==True
                expand=args.expand) # for del