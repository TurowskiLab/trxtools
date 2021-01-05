import TTools as tt
import os, argparse
from argparse import RawTextHelpFormatter
import pandas as pd

usage = "Generates profiles for a given list of genes from SAM files "

parser = argparse.ArgumentParser(usage=usage, formatter_class=RawTextHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="sam_file", help="Provide the path to your SAM file",
                 metavar="FILE", default=None)
files.add_argument("-l", dest="list_file", help="Provide the path to your list of genes file",
                 metavar="FILE", default='gene_list.list')
files.add_argument("-d", dest="details_file", help="Provide the path to file with transcriptome details.",
                   metavar="FILE", default='/homes/tollervey/COVID19/databases/transcriptome_details_biomart_tRNA_rRNA_LongNAME.tab')
files.add_argument("-s", dest="sequence_file", help="Provide the path to your sequence file.",
                   metavar="FILE", default='/homes/tollervey/COVID19/databases/transcriptome_h38_September2020/biomart_transcriptome_all_short_header_rRNA_tRNA.tab')
files.add_argument("-m", dest="method", help="Method how to generate profile: {read, middle, deletion}",
                   type=str, default='read')
files.add_argument("-e", dest="expand", help="For method {middle, deletion} position can be expanded by the value on each side (e=5 gives 10 nt long)", type=int, default=0)
files.add_argument("-c", dest="toClear", help="String of signs to be cleared from the name of SAM file", type=str, default='_comp_flexbar_STARAligned.out')
args = parser.parse_args()

df_details = pd.read_csv(args.details_file, sep="\t",index_col=0,header=0)
df_details = pd.read_csv(args.sequence_file, sep="\t",names=['name','sequence'],index_col=0)

path = os.getcwd()

with open(path+"/"+args.list_file) as f:
    gene_list = f.read().splitlines()

tt.SAM2profiles.sam2profiles(filename=args.sam_file, path=path,
             geneList=gene_list, toClear=args.toClear,
            df_details=df_details,df_sequences=df_details,
            how=args.method,expand=args.expand)