#!/usr/bin/env python

import os, argparse, gzip
import pandas as pd
from Bio import SeqIO
import trxtools as tt

### ARGPARSE ###

usage = "Convert FASTA file to sliding windows as sequence and its reverse complement (RC). Output is saved in a folder which name starts with 'window', with the window size. Each chromosome is saved in a separate file. If tab_output is selected, the output is saved as tab-separated file, otherwise as fasta files."

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
params = parser.add_argument_group('Parameters')
params.add_argument("-f","--file", dest="fasta_file", help="FASTA file", metavar="FILE", default=None, required=True)
params.add_argument("-g","--gzip", dest="compression", help="gzip compression",action="store_true")
params.add_argument("-w","--window", dest="window", help="sliding window size", type=int, default=65)
params.add_argument("-t","--temp", dest="temp", help="temperature", type=int, default=37)
params.add_argument("-s","--strand", dest="strand", choices=["plus", "minus", 'both'], help="strand", default="both")
params.add_argument("--tab", dest="tab_output", help="Save output as tab", action="store_true")
args = parser.parse_args()

path = os.getcwd()+"/"

### FUNCTIONS ###

def fasta_to_dataframe(fasta_file, compression=False):
    if compression:
        with gzip.open(fasta_file, "rt") as handle:
            records = SeqIO.parse(handle, "fasta")
    else:
        records = SeqIO.parse(fasta_file, "fasta")
    
    data = [(record.id, str(record.seq)) for record in records]
    df = pd.DataFrame(data, columns=["id", "sequence"])
    return df

def saveOutput(df, out_path, subname, tab_output=False):
    file_out_path = out_path + "/" + subname

    if tab_output:
        df.to_csv(file_out_path+".tab", sep='\t', header=False)
    else:
        with open(file_out_path+".fasta", "w") as output_handle:
            for index, row in df.iterrows():
                output_handle.write(">"+index+"\n"+row['sequence']+"\n")  

### MAIN ###

print("Reading FASTA file")
df1_fasta = fasta_to_dataframe(path + args.fasta_file, compression=args.compression)

window = args.window
print("Sliding window size: ", window)
out_path = path+"window"+str(window)
os.makedirs(out_path, exist_ok=True)

for n, (chromosome, sequence) in df1_fasta.iterrows():
    print("Processing chromosome: ", chromosome)
    if args.strand == 'plus' or args.strand == 'both':
        df2_to_save = tt.nascent.prepareNascent(sequence, name=chromosome, 
                                                  strand="plus", window=window, temp=args.temp)
        saveOutput(df2_to_save, out_path, subname=chromosome+"_plus", tab_output=args.tab_output)

    if args.strand == 'minus' or args.strand == 'both':
        df2_to_save = tt.nascent.prepareNascent(sequence, name=chromosome, 
                                                  strand="minus", window=window, temp=args.temp)
        saveOutput(df2_to_save, out_path, subname=chromosome+"_minus", tab_output=args.tab_output)

    # df2_to_save = sliding_sequence(sequence, name=chromosome, window=window)

print("Done")