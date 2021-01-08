#!/usr/bin/env python

import TTools as tt
import os, argparse
import pandas as pd

usage = "Generates profiles for a given list of genes from SAM files "

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="sam_file", help="SAM file",
                 metavar="FILE", default=None)
files.add_argument("-l", dest="list_file", help="list of genes file",
                 metavar="FILE", default=None)
files.add_argument("-d", dest="details_file", help="file with transcriptome details.",
                   type=str, default='/homes/tollervey/COVID19/databases/transcriptome_details_biomart_tRNA_rRNA_LongNAME_tRNA_corrected.tab')
files.add_argument("--del", dest="deletions", help="Generate additional profiles for deletions",action="store_true")
files.add_argument("-e", dest="expand", help="For deletions position can be expanded by the value on each side (e=5 gives 10 nt long)", type=int, default=5)
files.add_argument("-c", dest="toClear", help="String of signs to be cleared from the name of SAM file", type=str, default='_comp_flexbar_STARAligned.out')
args = parser.parse_args()


print("Start at "+tt.methods.timestamp())
df_details = pd.read_csv(args.details_file, sep="\t",index_col=0,header=0)

path = os.getcwd()+"/"

with open(path+args.list_file) as f:
    gene_list = f.read().splitlines()

tt.SAM2profiles.sam2profiles(filename=args.sam_file, path=path,
             geneList=gene_list, toClear=args.toClear,df_details=df_details,
             deletions=args.deletions ,expand=args.expand)

print("Done at "+tt.methods.timestamp())