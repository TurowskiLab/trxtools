#!/usr/bin/env python

import os, argparse
from platform import win32_edition
import pandas as pd

usage = "Generate GTF for given genome with fixed-size windows"

parser = argparse.ArgumentParser(usage=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
files = parser.add_argument_group('Options for input files')
files.add_argument("-f", dest="len_file", help="File with genome length",
                 metavar="FILE", default=None)
files.add_argument("-w", dest="window", help="Window size",default=100)
args = parser.parse_args()

path = os.getcwd()+"/"
filename=args.len_file

df = pd.read_csv(path+filename, sep="\t", names=["chrName",'length'], index_col=0)

def gtfFromLength(s1=pd.Series,w=100,strand="+"):
    
    df_out = pd.DataFrame()
    for i,length in s1.iteritems():
        starts = pd.Series([1]+[i for i in range(w+1,length-w,w)])
        stops  = pd.Series([i for i in range(w,length-w,w)]+[length])         
        df_temp = pd.DataFrame({
            'seqname' : i,
            'source' : 'window'+str(w),
            'feature' : 'exon',
            'start' : starts,
            'end' : stops,
            'score' : ".",
            'strand' : strand,
            'frame' : ".",
            'attribute' : 'gene_id '+i+'|window|'+starts.astype(str)+'|'+stops.astype(str)+'|'+strand+'; gene_name '+i+'|window|'+starts.astype(str)+'|'+stops.astype(str)+'|'+strand
        })

        df_out = pd.concat([df_out, df_temp],ignore_index=True)

    return df_out

w = int(args.window) #window
df_plus = gtfFromLength(df['length'], w=w, strand="+")
df_minus = gtfFromLength(df['length'], w=w, strand="-")

df_out = pd.concat([df_plus,df_minus])

df_out.to_csv(path+filename.replace('.len','')+"_window"+str(w)+".gtf",sep='\t',index=False,header=False)