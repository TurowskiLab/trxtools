import pandas as pd
import numpy as np
import os, collections, shutil, re
import TTools as tt

####################################################
#                 support functions                #
####################################################

def groupCIGAR(cigar_string=""):
    '''Split CIGAR string to list of tuples. 
    Takes ``"3S44M1S1H"`` and returns 
    ``[('3', 'S'), ('44', 'M'), ('1', 'S'), ('1', 'H')]``

    :param cigar_string: CIGAR string ""
    :type cigar_string: str
    :return: list of tuples
    :rtype: list
    '''

    return re.findall(r'(\d+)([A-Z]{1})', cigar_string)

def stripCIGAR(match=[], to_strip="H"):
    '''Removes H from output of groupCIGAR'''
    return [(n,l) for (n,l) in match if l != to_strip]

def tostripCIGARfive(match=[]):
    '''Nucleotides without alignment at the 5' end of the read

    :param match: output from groupCIGAR (list of tuples)
    :type match: list
    :return: number of substituted nucleotides
    :rtype: int
    '''
    if "S" in match[0]:
        return int(match[0][0])
    else:
        return None


def tostripCIGARthree(match=[]):
    '''Nucleotides without alignment at the 3' end of the read

    :param match: output from groupCIGAR (list of tuples)
    :type match: list
    :return: number of substituted nucleotides
    :rtype: int
    '''
    if "S" in match[-1]:
        return int(match[-1][0])
    else:
        return None


def stripSubstitutions(match):
    '''Strip substutiotns on both ends

    :param match: output from groupCIGAR (list of tuples)
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

def parseNoncoded(d=dict(), minLen=3):
    df_output = pd.DataFrame()
    for name in d.keys():
        df_temp = pd.DataFrame(d[name], columns=['position','end'])
        df_temp = df_temp[df_temp['end'].str.len() >= minLen].sort_values('end') #keep only minLen ends
        df_short = pd.DataFrame()
        for n,df in df_temp.groupby('end'):
            df_short = df_short.append(df.groupby('position')['end'].count().sort_index().rename(n))
        
        df_short = df_short.T
        df_short['chr'] = name
        df_output = pd.concat([df_output, df_short.reset_index()])

    return df_output.reset_index(drop=True)

def selectPolyA(df=pd.DataFrame()):
    #select columns with at least "AAA" and "A" content above 0.75
    cols = df.columns[(df.columns.str.contains("AAA"))&((df.columns.to_series().apply(tt.methods.letterContent)>0.75))]
    return df[['index','chr']+cols.tolist()]