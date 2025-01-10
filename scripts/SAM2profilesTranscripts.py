#!/usr/bin/env python
"""
SAM2profilesTranscripts.py
^^^^^^^^^^^^^^^^^^^^^^^^^^

This script generates profiles from SAM files for a list of given transcripts.
It reads a SAM file, a list of genes, and a file with transcriptome details,
and processes the data to create profiles. The output can be saved as a pickle
DataFrame if specified.

Usage:
    SAM2profilesTranscripts.py -f file.sam -l list_of_genes.txt -d transcriptome_details.tab

Options:
    -f FILE, --sam_file FILE
        SAM file to be processed.
    -l FILE, --list_file FILE
        File containing a list of genes.
    -d FILE, --details_file FILE
        File with transcriptome details. Default is '/homes/tollervey/COVID19/databases/transcriptome_details_biomart_tRNA_rRNA_UPDATEJan2021.tab'.
    --del
        Generate additional profiles for deletions.
    -p, --pickle
        Save output as a pickle DataFrame.
    -e INT, --expand INT
        For deletions, position can be expanded by the value on each side (e=5 gives 10 nt long). Default is 5.
    -c STRING, --toClear STRING
        String of signs to be cleared from the name of SAM file. Default is '_comp_flexbar_STARAligned.out'.
    --chunks INT
        Divide list of genes into chunks of given size. Default is 0.

Example:
    SAM2profilesTranscripts.py -f file.sam -l list_of_genes.txt -d transcriptome_details.tab

"""

import trxtools as tt
import os, argparse
import pandas as pd

usage = "Generates profiles from SAM files for a list of given transcripts."
usage += "For instrictions and additional options use SAM2profilesTranscripts.py --help \n"
usage += "Example usage: SAM2profilesTranscripts.py -f file.sam -l list_of_genes.txt -d transcriptome_details.tab"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="sam_file", help="SAM file",
                 metavar="FILE", default=None)
files.add_argument("-l", dest="list_file", help="list of genes file",
                 metavar="FILE", default=None)
files.add_argument("-d", dest="details_file", help="file with transcriptome details.",
                   type=str, default='/homes/tollervey/COVID19/databases/transcriptome_details_biomart_tRNA_rRNA_UPDATEJan2021.tab')
files.add_argument("--del", dest="deletions", help="Generate additional profiles for deletions",action="store_true")
files.add_argument("-p", dest="pickle", help="SAve output as pickle DataFrame",action="store_true")
files.add_argument("-e", dest="expand", help="For deletions position can be expanded by the value on each side (e=5 gives 10 nt long)", type=int, default=5)
files.add_argument("-c", dest="toClear", help="String of signs to be cleared from the name of SAM file", type=str, default='_comp_flexbar_STARAligned.out')
files.add_argument("--chunks", dest="chunks", help="Divide list of genes into a chunks of given size", type=int, default=0)
args = parser.parse_args()

print("Start at "+tt.methods.timestamp())
df_details = pd.read_csv(args.details_file, sep="\t",index_col=0,header=0)

path = os.getcwd()+"/"

with open(path+args.list_file) as f:
    gene_list = f.read().splitlines()

tt.sam.SAMtranscripts.sam2profiles(filename=args.sam_file, path=path,
             geneList=gene_list, toClear=args.toClear,df_details=df_details,
             deletions=args.deletions ,expand=args.expand,pickle=args.pickle,
             chunks=args.chunks)

print("Done at "+tt.methods.timestamp())