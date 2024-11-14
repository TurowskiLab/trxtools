#!/usr/bin/env python

import os, argparse, gzip
import pandas as pd
from Bio import SeqIO

### ARGPARSE ###

usage = "convert FASTA file to sliding windows"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
params = parser.add_argument_group('Parameters')
params.add_argument("-f","--file", dest="fasta_file", help="FASTA file", metavar="FILE", default=None, required=True)
params.add_argument("-g","--gzip", dest="compression", help="gzip compression",action="store_true")
params.add_argument("-w","--window", dest="window", help="sliding window size", type=int, default=65)
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

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = [complement.get(base, base) for base in seq]
    return ''.join(bases[::-1])

def sliding_sequence(sequence=str(), name='str', window=int()):
    df_temp = pd.DataFrame()
    index_temp = list()
    list_temp = list()

    #sliding window for the sequence
    for i in range(window,len(sequence)+1):
        string = sequence[i-window:i]
        list_temp.append(string)
        index_temp.append(name+'_nascent_len'+str(window)+'_pos'+str(i))

    #sliding window for reverse complement of the sequence
    sequence_revcomp = reverse_complement(sequence)
    for i in range(window,len(sequence_revcomp)+1):
        string_rc = sequence_revcomp[i-window+1:i+1] #correction to get exactly the same position for both strands
        list_temp.append(string_rc)
        index_temp.append(name+'_nascentRC_len'+str(window)+'_pos'+str(len(sequence_revcomp)-i))
    
    df_temp['string'] = pd.Series(list_temp, index=index_temp, dtype='str')
    
    return df_temp

### MAIN ###

print("Reading FASTA file")
df1_fasta = fasta_to_dataframe(path + args.fasta_file, compression=args.compression)

window = args.window
print("Sliding window size: ", window)
out_path = path+"window"+str(window)
os.makedirs(out_path, exist_ok=True)

for n, (chromosome, sequence) in df1_fasta.iterrows():
    print("Processing chromosome: ", chromosome)
    df2_to_save = sliding_sequence(sequence, name=chromosome, window=window)
    df2_to_save.to_csv(out_path+"/"+chromosome+".tab", sep="\t", header=False)

print("Done")