#!/usr/bin/env python
"""
genomeNascentFolding.py
^^^^^^^^^^^^^^^^^^^^^^^

This script processes a given FASTA file to generate sliding windows of sequences and their reverse complements. 
It then folds these sequences using RNAfold and parses the folding output to calculate the free energy (dG) of the folded sequences.
The results are saved in WIG and BigWig formats for visualization.

Usage:
    python genomeNascentFolding.py -f genome.fasta -w 65 -t 37 -s both

Arguments:

    -f, --file
      Path to the input FASTA file (required).

    -g, --gzip
      Enable gzip compression for the input file.

    -w, --window
      Size of the sliding window (default: 65).

    -t, --temp
      Temperature for RNA folding (default: 37).
    
    -s, --strand
      Strand to process ('plus', 'minus', or 'both'; default: 'both').

Output:
    The script generates sliding window sequences and their reverse complements, folds them using RNAfold, 
    and saves the folding free energy (dG) values in WIG and BigWig formats. The output files are saved in the same 
    directory as the input file, with appropriate naming conventions.

"""
import os, argparse, gzip
import pandas as pd
from Bio import SeqIO
import trxtools as tt
import re

### ARGPARSE ###

usage = "Convert FASTA file to sliding windows as sequence and its reverse complement (RC). \n"
usage += "Output is saved in a folder which name starts with 'window', with the window size. Each chromosome is saved in a separate file. If tab_output is selected, the output is saved as tab-separated file, otherwise as fasta files."
usage += "The output is saved in the same folder as the input file."
usage += "Example usage: python genomeNascentFolding.py -f genome.fasta -w 65 -t 37 -s both"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
params = parser.add_argument_group('Parameters')
params.add_argument("-f","--file", dest="fasta_file", help="FASTA file", metavar="FILE", default=None, required=True)
params.add_argument("-g","--gzip", dest="compression", help="gzip compression",action="store_true")
params.add_argument("-w","--window", dest="window", help="sliding window size", type=int, default=65)
params.add_argument("-t","--temp", dest="temp", help="temperature", type=int, default=37)
params.add_argument("-s","--strand", dest="strand", choices=["plus", "minus", 'both'], help="strand", default="both")
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

def saveSlidedFastaOutput(df, out_path, subname):
    file_out_path = out_path + "/" + subname

    with open(file_out_path+".fasta", "w") as output_handle:
        for index, row in df.iterrows():
            output_handle.write(">"+index+"\n"+row['sequence']+"\n")  

def parseFoldingOutput(input_file):
    pattern = re.compile(r'\(\s*(-?\d+\.\d+)\)')
    
    f = open(input_file, "r")
    
    parsed_dG = {}
    curr_seq_id = None
    count = 0

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            curr_seq_id = line[1:]
            count = 1

        elif count == 1: count = 2

        elif count == 2:
            match = pattern.search(line)
            if match:
                parsed_dG[curr_seq_id] = float(match.group(1))
                count = 0
            else:
                print(line)
    return parsed_dG

def processFolding(name="", working_path="",
                    window=int(), temp=int()):
    chr, strand = name.split("_")
    sequence = df1_fasta[df1_fasta['id'] == chr]['sequence'].values[0]
    length = len(sequence)


    ## Creating slided FASTA files
    print("Processing: ", chr, strand, "(sliding sequences)")
    df2_to_save = tt.nascent.prepareNascent(sequence, name=chr, 
                                                  strand=strand, window=window, temp=args.temp)
    saveSlidedFastaOutput(df2_to_save, working_path, subname=chr+"_"+strand)


    ## Folding sequences
    print("Processing: ", chr, strand, "(folding)") 
    outFolding = working_path + "/" + chr + "_" + strand + "_RNAfold.output"

    command =  "RNAfold --noPS -T " + str(args.temp) + " -i " + working_path + "/" + chr + "_" + strand + ".fasta"
    command += " > " + outFolding
    tt.methods.bashCommand(command)
    
    if os.path.getsize(outFolding) == 0:
        KeyError("ERROR: RNAfold failed for ", chr, strand, " Please check is RNAfold is installed and working properly.")

    ## Parsing folding output
    parseddG = parseFoldingOutput(outFolding)
    s1 = pd.Series(parseddG, name='dG')
    df_dG = pd.DataFrame(s1).reset_index()
    df_t1 = df_dG['index'].str.split("_", expand=True)
    df_t1.columns = ['chromosome','strand','pos','temp','len']

    #assign name elements to the df
    df_dG['chrom'] = df_t1['chromosome']
    df_dG['strand'] = df_t1['strand'].str.replace('nascentRC','-').str.replace('nascent','+')
    df_dG['win'] = df_t1['len'].str.replace('win','').astype(int)
    df_dG['pos'] = df_t1['pos'].str.replace('pos','').astype(int)

    #save the output
    df_dG.to_pickle(working_path + "/" + chr + "_" + strand + "_dG.pkl")
    print("Done: ", chr, strand, "(folding)")

    return df_dG, chr, length

def saveWIG(data=pd.DataFrame(),name='test',path=str(),chr_len=dict()):
    full_path = path+name+"_dG.wig"

    for chrom in chr_len.keys():
        print("Saving: ", chrom)
        #preparing body of chromosome
        df_saving = pd.DataFrame(pd.Series([0]*chr_len[chrom], name='empty'))
        df_saving['chrom'] = chrom
        df_saving['dG'] = data[data['chrom']==chrom]['dG']
        df_saving = df_saving.drop('empty',axis=1).fillna(0)

        #saving header
        # if chrom=="Mito": chrom = "chrM" #to visualize with the IGV
        with open(full_path,'w') as output_file:
            output_file.write("fixedStep chrom="+chrom+" start=1 step=1\n")
        #saving body of chromosome
        df_final = df_saving.drop('chrom',axis=1)
        with open(full_path,'a') as output_file:
            df_final.to_csv(output_file, index=None, header=None, sep="\t")

def WIGtoBigWig(wig_file,chrom_sizes_file):
    command = "wigToBigWig " + wig_file + " " + chrom_sizes_file + " " + wig_file.replace(".wig",".bw")
    tt.methods.bashCommand(command)

### MAIN ###

print("Reading FASTA file")
df1_fasta = fasta_to_dataframe(path + args.fasta_file, compression=args.compression)

## working folder
window = args.window
print("Sliding window size: ", window)
working_path = path+"window"+str(window)+"temp"+str(args.temp)
os.makedirs(working_path, exist_ok=True)

## list of working files
workingList = []
chr_len = {}
df_dG_all = pd.DataFrame()

if args.strand == 'plus' or args.strand == 'both':
     workingList = workingList + [i+"_plus" for i in df1_fasta['id'].tolist()]
if args.strand == 'minus' or args.strand == 'both':
     workingList = workingList + [i+"_minus" for i in df1_fasta['id'].tolist()]

for i in workingList:
    df_dG, chr, length = processFolding(name=i, working_path=working_path,
                    window=args.window, temp=args.temp)
    chr_len[chr] = length
    df_dG_all = pd.concat([df_dG_all, df_dG],ignore_index=True)

print("Saving WIG files")
df_dG_all = df_dG_all.drop(['index','win'],axis=1)
#selecting and indexing strands
df_t2_strand = df_dG_all[df_dG_all['strand'] == "plus"]
df_t3_strandRC = df_dG_all[df_dG_all['strand'] == "minus"]
df_t2_strand = df_t2_strand.drop("strand", axis=1).set_index('pos').sort_index()
df_t3_strandRC = df_t3_strandRC.drop("strand", axis=1).set_index('pos').sort_index()

#saving as WIG files
n = args.fasta_file.replace('.fasta','_')
saveWIG(df_t2_strand, name=n+'plus_window'+str(window),path=os.getcwd()+"/",chr_len=chr_len)
saveWIG(df_t3_strandRC, name=n+'minus_window'+str(window),path=os.getcwd()+"/",chr_len=chr_len)

# converting to BigWig
print("Converting WIG to BigWig")

with open(os.getcwd() + "/" + n + "chr_len.tab", "w") as f:
    for key, value in chr_len.items():
        f.write(f"{key}\t{value}\n")

WIGtoBigWig(wig_file=os.getcwd()+"/"+n+"plus_window"+str(window)+"_dG.wig",
            chrom_sizes_file=os.getcwd() + "/" + n + "chr_len.tab")
WIGtoBigWig(wig_file=os.getcwd()+"/"+n+"minus_window"+str(window)+"_dG.wig",
            chrom_sizes_file=os.getcwd() + "/" + n + "chr_len.tab")


print("Done.")