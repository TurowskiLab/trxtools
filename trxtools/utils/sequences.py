import subprocess, time, random
import os, sys, re, itertools
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from trxtools.plotting import clusterClusterMap

################################################
#############        DNA/RNA sequences

def calGC(dataset=pd.DataFrame(), calFor=['G','C']):
    '''Returns GC content in a given string - uses ['nucleotide'] column

    :param dataset: DataFrame() with "nucleotide" column
    :type dataset: pandas.DataFrame
    :param calFor: list of nucleotides to calculate GC content
    :type calFor: list

    :return: fraction of GC content, float
    :rtype: float
    '''
    return float(len(dataset[dataset['nucleotide'].isin(calFor)]))/float(len(dataset))

alt_map = {'ins':'0'}
complement_DNA = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
complement_RNA = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}

def reverse_complement_DNA(seq):
    '''Reverse complement of DNA sequence

    :param seq: DNA sequence to reverse complement
    :type seq: str

    :return: Reverse complement of the input sequence
    :rtype: str

    :example:

    >>> reverse_complement_DNA("ATCG")
    'CGAT'
    '''

    if "U" in seq: return str()
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement_DNA.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def reverse_complement_RNA(seq):
    '''Reverse complement of RNA sequence

    :param seq: RNA sequence to reverse complement
    :type seq: str

    :return: Reverse complement of the input sequence
    :rtype: str

    :example:

    >>> reverse_complement_RNA("AUCG")
    'CGAU'
    '''

    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement_RNA.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def reverse_complement(seq):
    '''Reverse complement of DNA or RNA sequence. Identifies RNA by presence of 'U' in the sequence.

    :param seq: DNA or RNA sequence to reverse complement
    :type seq: str

    :return: Reverse complement of the input sequence
    :rtype: str

    :example:

    >>> reverse_complement("AUCG")
    'CGAU'
    '''
    if "U" in seq:
        return reverse_complement_RNA(seq)
    else:
        return reverse_complement_DNA(seq)

def randomDNAsingle(length=int(), letters="CGTA"):
    '''Random generator of nucleotide sequence

    :param length: Length of the sequence
    :type length: int
    :param letters: str with letters that will be used
    :type letters: str

    :return: Random sequence of the given length
    :rtype: str

    :example:

    >>> randomDNAsingle(10,"CGTA")
    'CGTACGTAAC'
    '''

    return''.join(random.choices(letters, k=length))

def randomDNAall(length=int(), letters="CGTA"):
    '''Generates all possible random sequences of a given length

    :param length: Length of the sequence
    :type length: int
    :param letters: str with letters that will be used
    :type letters: str

    :return: List of all possible random sequences of the given length
    :rtype: list

    :example:

    >>> randomDNAall(2,"CGTA")
    ['CC', 'CG', 'CT', 'CA', 'GC', 'GG', 'GT', 'GA', 'TC', 'TG', 'TT', 'TA', 'AC', 'AG', 'AT', 'AA']
    '''

    output_list = []
    for i in itertools.product(list(letters), repeat=length):
        output_list.append("".join(i))
    return output_list

def rollingGC(s=pd.Series, window=10): #rolling window, smoothing of data
    '''Calculates GC from sequence, uses 'boxcar' window

    :param s: Series containing sequence
    :type s: pandas.Series
    :param window: window size for GC calculation
    :type window: int

    :return: Series with GC calculated, center=False
    :rtype: pandas.Series
    '''
    #be aware of win_type and center parameters
    return s.replace(['G','C'],1).replace(['T','A'],0).rolling(window=window, win_type='boxcar',center=False).mean()

def letterContent(seq=str(), letter="A"):
    '''Calculates the content of a given letter in a sequence

    :param seq: Sequence to calculate the content of the given letter in
    :type seq: str
    :param letter: letter to calculate the content of in the sequence, defaults to "A"
    :type letter: str

    :return: Content of the given letter in the sequence, round to two decimal places
    :rtype: float

    :example:

    >>> letterContent("ATCG", "A")
    0.25
    '''

    return round(len([i for i in seq if i==letter])/len(seq),2)

def DNA_string_positions(df, string, name_col='name', seq_col='sequence'):
    '''Finds all occurences of a given substring in a string (i.e. DNA sequence) and reports their positions

    :param df: Input dataframe
    :type df: pandas.DataFrame
    :param string: Query string
    :type string: str
    :param name_col: name of column containing sequence names, defaults to 'name'
    :type name_col: str, optional
    :param seq_col: name of column containing sequences, defaults to 'sequence'
    :type seq_col: str, optional

    :return: Dataframe in long format containing positionsof found string
    :rtype: pandas.DataFrame
    '''

    out_df = pd.DataFrame({name_col:[], 'length':[]})
    for index, row in df.iterrows():
        positions = [m.start() for m in re.finditer(string, row[seq_col])]
        results = pd.DataFrame(
            {
                'name': [row[name_col]] * len(positions),
                'start': positions,
                'length': [len(string)] * len(positions)
            }
        )
        out_df = pd.concat([out_df, results], ignore_index=True)
    return out_df

def is_inside(inner_start, inner_end, outer_start, outer_end):
    '''Check if one region is inside another string.

    :param inner_start: Start position of the inner region
    :type inner_start: int
    :param inner_end: End position of the inner region
    :type inner_end: int
    :param outer_start: Start position of the outer region
    :type outer_start: int
    :param outer_end: End position of the outer region
    :type outer_end: int

    :return: True if the inner region is inside the outer region, False otherwise
    :rtype: bool

    :example:

    >>> is_inside(1, 3, 0, 5)
    True
    >>> is_inside(1, 3, 0, 2)
    False
    '''

    left_border = ((inner_start == outer_start) & (inner_end < outer_end))
    inside = ((inner_start > outer_start) & (inner_end < outer_end))
    right_border = ((inner_start > outer_start) & (inner_end == outer_end))
    return (left_border | inside | right_border)

def nested_region_cleanup(df, start_col='start', end_col='end'):
    '''Helper function to remove shorter regions nested in longer ones

    :param df: input dataframe
    :type df: pandas.DataFrame
    :param start_col: name of column containing region start positions, defaults to 'start'
    :type start_col: str, optional
    :param end_col: name of column containing region end positions, defaults to 'end'
    :type end_col: str, optional
    :param length_col: name of column containing region lengths, defaults to 'length'
    :type length_col: str, optional

    :return: dataframe with nested regions removed
    :rtype: pandas.DataFrame
    '''

    drop_series = []
    out_df = df.copy()
    # out_df['length'] = out_df[end_col] - out_df[start_col]
    out_df = out_df.sort_values('length')
    for index, row in out_df.iterrows():   
        #in each iteration compare all other rows to current row
        for index2, row2 in out_df.iterrows():
            if is_inside(row2[start_col], row2[end_col], row[start_col], row[end_col]):
                drop_series.append(index2)
                # print(index, index2)
    drop_series = pd.Series(drop_series).unique()
    # print(drop_series)
    out_df = out_df.drop(drop_series)
    return out_df
        

def DNA_stretch_positions(df, char, min_len, max_len, name_col='name', seq_col='sequence'):
    '''Finds all character stretches of given length in a string (i.e. DNA sequence) and reports their positions. Wrapper for DNA_string_positions()

    :param df: Input dataframe
    :type df: pandas.DataFrame
    :param char: Query character (base)
    :type char: str
    :param min_len: minimum length of stretch to look for
    :type min_len: int
    :param min_len: maximum length of stretch to look for
    :type min_len: int
    :param name_col: name of column containing sequence names, defaults to 'name'
    :type name_col: str, optional
    :param seq_col: name of column containing sequences, defaults to 'sequence'
    :type seq_col: str, optional

    :return: Dataframe in long format containing positions of found string
    :rtype: pandas.DataFrame
    '''
 
    out_df = pd.DataFrame({name_col: [], 'length': []})
    for i in range(min_len, max_len+1):
        string = i*char
        out_df = pd.concat([out_df, DNA_string_positions(
            df, string, name_col=name_col, seq_col=seq_col)], ignore_index=True)
    # out_df =  pd.DataFrame(out_df.groupby(['name', 'start'])['length'].max()).reset_index()
    out_df['end'] = out_df['start'] + out_df['length']
    return out_df.groupby(name_col, group_keys=False).apply(nested_region_cleanup)

def find_pol3_terminators(df, min_T, max_T, name_col='name', seq_col='sequence'):
    '''Finds all Pol3 terminators (T stretches) of given length in a string (i.e. DNA sequence) and reports their positions. Wrapper for DNA_string_positions()

    :param df: Input dataframe
    :type df: pandas.DataFrame
    :param min_T: minimum length of T stretch to look for
    :type min_T: int
    :param min_T: maximum length of T stretch to look for
    :type min_T: int
    :param region_col: name of column containing sequence names, defaults to 'region'
    :type region_col: str, optional
    :param seq_col: name of column containing sequences, defaults to 'sequence'
    :type seq_col: str, optional

    :return: Dataframe in long format containing positionsof found string
    :rtype: pandas.DataFrame
    '''

    out_df = DNA_stretch_positions(df, char='T', min_len=min_T, max_len=max_T, name_col=name_col, seq_col=seq_col)
    return out_df