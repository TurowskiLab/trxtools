import pandas as pd
import numpy as np
import os, collections, shutil, re
import trxtools as tt
import pyBigWig

####################################################
#                 support functions                #
####################################################

def groupCIGAR(cigar_string=""):
    '''Split CIGAR string to list of tuples

    :param cigar_string:
    :type cigar_string: str
    :return: list of tuples ``[( ),()]``
    :rtype: list

    >>> groupCIGAR("3S44M1S1H")
    [('3', 'S'), ('44', 'M'), ('1', 'S'), ('1', 'H')]
    '''

    return re.findall(r'(\d+)([A-Z]{1})', cigar_string)

def stripCIGAR(match=[], to_strip="H"):
    '''Removes H from output of groupCIGAR

    :param match: output of groupCIGAR, defaults to []
    :type match: list
    :param to_strip: CIGAR mark to be stripped, defaults to "H"
    :type to_strip: str, optional
    :return: modified list of tuples
    :rtype: list
    '''

    return [(n,l) for (n,l) in match if l != to_strip]

def tostripCIGARfive(match=[]):
    '''Calculates length of soft-clipped nucleotides at the 5' end of the read

    :param match: output of groupCIGAR, defaults to []
    :type match: list
    :return: number of substituted nucleotides
    :rtype: int

    >>> tostripCIGARfive([('3', 'S'), ('44', 'M'), ('1', 'S'), ('1', 'H')])
    3
    '''

    if "S" in match[0]:
        return int(match[0][0])
    else:
        return None


def tostripCIGARthree(match=[]):
    '''Nucleotides without alignment at the 3' end of the read

    :param match: output of groupCIGAR, defaults to []
    :type match: list
    :return: number of substituted nucleotides
    :rtype: int

    >>> tostripCIGARthree([('3', 'S'), ('44', 'M'), ('1', 'S')])
    1
    '''

    if "S" in match[-1]:
        return int(match[-1][0])
    else:
        return None


def stripSubstitutions(match):
    '''Strip substutiotns on both ends

    :param match: output of groupCIGAR, defaults to []
    :type match: list
    :return: list of tuples
    :rtype: list
    '''

    if "S" in match[0]:
        match = match[1:]
    if "S" in match[-1]:
        match = match[:-1]
    return match


def countRead(i=tuple()):
    '''Takes tuple (position,CIGARstring) and returns list of mapped positions

    :param i: tuple (first position of the read ,CIGAR string)
    :type i: tuple
    :return: list of mapped positions
    :rtype: np.array

    >>> countRead((400,"3S15M1D9M2S"))
    array([400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412,
       413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424])
    '''

    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    output = np.arange(length, dtype=int)
    return output + position


def countMiddle(i=tuple(), expand=0):
    '''Takes tuple (position,CIGARstring) and returns list with the middle of mapped read

    :param i: tuple (first position of the read ,CIGAR string)
    :type i: tuple
    :param expand: number of nucleotides to expand for each side, defaults to 0
    :type expand: int, optional
    :return: list of mapped positions
    :rtype: np.array

    >>> countMiddle((400,"3S15M1D9M2S"))
    array([412])
    >>> countMiddle((400,"3S15M1D9M2S"),expand=3)
    array([409, 410, 411, 412, 413, 414])
    '''

    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    mid = round(length / 2)
    if expand > 0:
        output = np.arange(mid - expand, mid + expand)
        return output + position
    else:
        return np.arange(mid, mid + 1) + position

def countDeletion(i=tuple(), expand=0):
    '''Takes tuple (position,CIGARstring) and returns list of mapped deletions. Each deletion is counted once, but expand parameter can merge some longer deletions (common expansion).

    :param i: tuple (first position of the read ,CIGAR string)
    :type i: tuple
    :param expand: number of nucleotides to expand for each side, defaults to 0
    :type expand: int, optional
    :return: list of mapped positions
    :rtype: np.array

    >>> countDeletion((400,"3S15M1D9M2S"))
    array([415])
    >>> countDeletion((400,"3S15M1D9M2S"),expand=3)
    array([412, 413, 414, 415, 416, 417, 418])
    '''

    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    output = np.arange(length, dtype=int)
    outputMask = np.zeros(length, dtype=int)
    cumsum = 0
    for i, x in match:
        i = int(i)
        if x == 'D' and expand == 0:
            outputMask[cumsum:cumsum + i] = 1
        elif x == 'D' and expand > 0:
            outputMask[cumsum - expand:cumsum + i + expand] = 1
        cumsum += i
    output = (output + position) * outputMask
    return output[output != 0]

# def parseNoncoded(d=dict(), minLen=3):
#     '''Parse dict with non-coded ends and returns structured DataFrame

#     :param d: dictionary with list of tuples for each chromosme ``{"chrI" : [], "chrII" : []}``, defaults to dict()
#     :type d: dict
#     :param minLen: minimal length for non-coded end to keep, defaults to 3
#     :type minLen: int, optional
#     :return: DataFrame with parsed non-coded ends
#     :rtype: DataFrame

#     >>> parseNoncoded({"chrI":[(40, 'AAA'), (35, 'AACAA')]})
#        index  AAA  AACAA   chr
#     0     40  1.0    NaN  chrI
#     1     35  NaN    1.0  chrI
#     '''

#     df_output = pd.DataFrame()
#     for name in d.keys():
#         df_temp = pd.DataFrame(d[name], columns=['position','end'])
#         df_temp = df_temp[df_temp['end'].str.len() >= minLen].sort_values('end') #keep only minLen ends
        
#         df_short = pd.DataFrame()
#         for n,df in df_temp.groupby('end'):
#             # df_short = df_short.append(df.groupby('position')['end'].count().sort_index().rename(n))
#             s = df.groupby('position')['end'].count().sort_index().rename(n)
#             df_short = pd.concat([df_short,s],axis=1) #append
#             # df_short = pd.concat([df_short,s]) #append
#         df_short['chr'] = name
#         df_output = pd.concat([df_output, df_short.reset_index()])

#     return df_output.reset_index(drop=True)

# def parseNoncodedList(l=[], minLen=3):
#     '''Parse list with non-coded ends and returns structured DataFrame

#     :param l: list of tuples ``[(int,str)]``, defaults to lits()
#     :type l: list
#     :param minLen: minimal length for non-coded end to keep, defaults to 3
#     :type minLen: int, optional
#     :return: DataFrame with parsed non-coded ends
#     :rtype: DataFrame

#     >>> parseNoncodedList([(40, 'AAA'), (35, 'AACAA')])
#         AAA	AACAA
#     35	NaN	1.0
#     40	1.0	NaN
#     '''
#     df_temp = pd.DataFrame(l, columns=['position','end'])
#     df_temp = df_temp[df_temp['end'].str.len() >= minLen].sort_values('end') #keep only minLen ends
    
#     df_short = pd.DataFrame()
#     for n,df in df_temp.groupby('end'):
#         s = df.groupby('position')['end'].count().sort_index().rename(n)
#         df_short = pd.concat([df_short,s],axis=1) #append

#     return df_short

# def selectPolyA(df=pd.DataFrame()):
#     '''Select only polyA non-coded ends containinig ``"AAA"`` and "A"-content above 75%

#     :param df: output of parseNoncoded or parseNoncodedList
#     :type df: DataFrame
#     :return: modified DataFrame
#     :rtype: DataFrame
#     '''
#     #select columns with at least "AAA" and "A" content above 0.75
#     cols = df.columns[(df.columns.str.contains("AAA"))&((df.columns.to_series().apply(tt.methods.letterContent)>0.75))]
#     if 'chr' in df.columns.values:
#         return df[['index','chr']+cols.tolist()]
#     else:
#         return df[cols.tolist()]
    
def selectNoncodedAndProfile(l=[],minLen=3,tail="AAA",letter="A",content=0.80):
    '''Select sequences with non-coded ends and profile their occurrence.

    :param l: list of tuples where each tuple contains a position of 3' end of read match and a sequence string, defaults to []
    :type l: list
    :param minLen: minimum length of the sequence to be considered, defaults to 3
    :type minLen: int, optional
    :param tail: the non-coded end sequence to look for, defaults to "AAA"
    :type tail: str, optional
    :param letter: the letter to check the content of in the sequence, defaults to "A"
    :type letter: str, optional
    :param content: the minimum content threshold of the letter in the sequence, defaults to >=0.80
    :type content: float, optional

    :return: a pandas Series profiling the occurrence of selected sequences
    :rtype: pd.Series

    :example:

    >>> l = [(1, "AAA"), (2, "AACAA")]
    >>> selectNoncodedAndProfile(l, minLen=3,tail="AAA",letter="A",content=0.80)
    1    1
    '''
    l_output = []
    for i in l:
        if tail in i[1]: #check if contains non-coded end
            if tt.methods.letterContent(i[1],letter)>=content: #check if content of letter is above threshold
                l_output.append(i[0])
    
    profile = pd.Series(collections.Counter(l_output), dtype="int").sort_index().astype(float)  # faster that using zip and numpy

    return profile
    

def selectEnds(df=pd.DataFrame(),ends="polyA"):
    '''Wrapper for functions selecting non-coded ends

    :param df: output of parseNoncoded or parseNoncodedList
    :type df: DataFrame
    :param ends: type of ends, currently only "polyA" is availible, defaults to "polyA"
    :type ends: str
    :return: runs selectPolyA
    :rtype: DataFrame
    '''
    if ends=="polyA":
        return selectPolyA(df=df)

def noncoded2profile(df_input=pd.DataFrame(), df_details=pd.DataFrame()):
    '''Turns non-coded ends into profile

    :param df_input: output of parseNoncoded function
    :type df_input: DataFrame
    :param df_details: chromosome lengths
    :type df_details: DataFrame
    :return: DataFrame with profiles
    :rtype: DataFrame
    '''
    
    df_output = pd.DataFrame()
    for i,df in df_input.groupby('chr'):
        df = df.drop('chr',"columns").set_index('index')
        profile = df.sum(1).astype(float)
        length = df_details.loc[i]['length']
        # profile = profile.reindex(pd.RangeIndex(length + 1)).fillna(0.0)  # fills spaces with 0 counts
        
        # df_output = df_output.append(profile.rename(i))
        df_output = pd.concat([df_output,profile.rename(i)],axis=0, join='outer')
        
    return df_output

def noncoded2profile1(df=pd.DataFrame(), length=int()):
    '''Turns non-coded ends into profile

    :param df_input: output of parseNoncodedList function
    :type df_input: DataFrame
    :param df_details: chromosome lengths
    :type df_details: DataFrame
    :return: Series with profiles
    :rtype: Series
    '''
    
    profile = df.sum(1).sort_index().astype(float)
    # profile = profile.reindex(pd.RangeIndex(length + 1)).fillna(0.0)  # fills spaces with 0 counts
    return profile

def saveBigWig(paths=dict(),suffix=str(),bw_name=str(),chroms=list(),pkl=True):
    '''Save gzip pickle data to BigWig
    This function saves data from gzip pickle files to a BigWig file.

    :param paths: Dictionary with paths to gzip pickle files, where keys are file paths and values are file contents.
    :type paths: dict
    :param suffix: Suffix to be removed from chromosome names in the paths, defaults to an empty string.
    :type suffix: str, optional
    :param bw_name: Name of the output BigWig file.
    :type bw_name: str
    :param chroms: List of tuples with chromosome names and their lengths.
    :type chroms: list of tuples

    :return: Success message indicating the suffix and that the file was saved successfully.
    :rtype: str

    :example:

    >>> paths = {'/path/to/file_chr1.pkl.gz': 'data1', '/path/to/file_chr2.pkl.gz': 'data2'}
    >>> suffix = '3end'
    >>> bw_name = 'output.bw'
    >>> chroms = [('chr1', 1000000), ('chr2', 2000000)]
    >>> saveBigWig(paths, suffix, bw_name, chroms)
    '3end saved successfully'
    '''

    bw = pyBigWig.open(bw_name, "w")
    bw.addHeader(chroms)
    if pkl == True:
        for p in paths.keys():
            c = p.replace(suffix,"")
            df = pd.read_pickle(paths[p],compression='gzip')
            try:
                stops = df.index.to_numpy()
                starts = stops-1
                vals=df[c].to_numpy()
                bw.addEntries([c] * len(starts), starts, ends=stops, values=vals)
            except: #dealing with potential problems - unclear the origin of the problems
                print(c+" using except")
                df = df.sort_index()
                stops = df.index.to_numpy()
                starts = list(stops-1)
                if len(starts)==0: continue
                vals=df[c].to_list()
                bw.addEntries([c] * len(starts), starts, ends=stops.tolist(), values=vals)
    else:
        for p in paths.keys():
            c = p.replace(suffix,"")
            vals = paths[p]
            vals = pd.Series(vals).dropna()
            starts = vals.index.to_numpy()
            stops = starts+1
            # stops = np.arange(1,chr_len+1)
            # starts = stops-1
            bw.addEntries([c] * len(starts), starts, ends=stops, values=vals.to_numpy())
            
    bw.close()
    return (suffix+" saved succesfully")

def selectSortPaths(paths={},chroms=[],suffix=""):
    '''Select and sort paths based on chromosome names and suffix
    This function filters and sorts the paths dictionary based on the provided chromosome names and suffix.

    :param paths: Dictionary with paths to gzip pickle files, where keys are file paths and values are file contents.
    :type paths: dict
    :param chroms: List of tuples with chromosome names and their lengths.
    :type chroms: list of tuples
    :param suffix: Suffix to be matched in the file paths, defaults to an empty string.
    :type suffix: str, optional

    :return: Dictionary with filtered and sorted paths.
    :rtype: dict

    :example:

    >>> paths = {'/path/to/file_chr1_3end.pkl.gz': 'data1', '/path/to/file_chr2_3end.pkl.gz': 'data2'}
    >>> chroms = [('chr1', 1000000), ('chr2', 2000000)]
    >>> suffix = '3end'
    >>> selectSortPaths(paths, chroms, suffix)
    {'/path/to/file_chr1_3end.pkl.gz': 'data1', '/path/to/file_chr2_3end.pkl.gz': 'data2'}
    '''

    out_dict = {}
    for i in dict(chroms).keys():
        for p in paths.keys():
            if (p.endswith(suffix)) and (i+"_" in p):
                out_dict[p] = paths[p]
    return out_dict