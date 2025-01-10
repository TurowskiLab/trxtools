import subprocess, time, random
import os, sys, re, itertools
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import trxtools as tt

### METHODS ###
# various python methods

################################################
#############        bash

def listPaths(folder=".", suffix=str(), nameElem=None):
    '''List all file paths in a folder with a specific suffix and optional name element.

    :param folder: Folder to list files from, defaults to "."
    :type folder: str, optional
    :param suffix: Suffix of the files to list, defaults to an empty string
    :type suffix: str, optional
    :param nameElem: Optional element that must be in the file name, defaults to None
    :type nameElem: str, optional

    :return: List of file paths that match the criteria
    :rtype: list

    :example:

    >>> listPaths(folder=".", suffix=".txt", nameElem="file")
    ['file1.txt', 'file2.txt']
    '''
    if nameElem:
        return [folder + file for file in os.listdir(folder) if file.endswith(suffix) and nameElem in file]
    else:
        return [folder + file for file in os.listdir(folder) if file.endswith(suffix)]

def bashCommand(bashCommand=str()):
    '''Run command in bash using subprocess.call()
    
    :param bashCommand: str with command

    :return: None
    '''
    pipes = subprocess.Popen(bashCommand,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
            err_msg = "%s. Code: %s" % (std_err.strip(), pipes.returncode)
            raise Exception(err_msg)

################################################
#############        DNA/RNA

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
################################################
#############        importing data and handling files

def read_list(filepath=str()):
    '''Read list from file. Each row becomes item in the list.
    
    :param filepath: Path to the file
    :type filepath: str

    :return: List of items from the file
    :rtype: list

    :example:

    >>> read_list("file.txt")
    ['item1', 'item2']
    '''

    txt_file = open(filepath, "r")
    content_list = txt_file.read().splitlines()
    txt_file.close()
    return content_list

def read_tabFile(nameElem="", path="", toLoad="", toClear=[], toAdd="", df=None, overwrite=False):
    '''Read tab files with common first column

    :param nameElem: part of a filename which is present in all files
    :type nameElem: str
    :param path: path to directory with files
    :type path: str
    :param toLoad: to be present in file name
    :type toLoad: str
    :param toClear: part of filename to be removed from file name
    :type toClear: str
    :param toAdd: string to be added to file name
    :type toAdd: str
    :param df: DataFrame, to be appended; default=None
    :type df: pandas.DataFrame
    :param overwrite: boolean, allows for overwriting during appending, default = False
    :type overwrite: bool

    :return: dataframe with all files, where columns are values from 'name' column
    :rtype: pandas.DataFrame
    '''

    # list files with STAT mapping
    l1_mapping = [f for f in os.listdir(path) if nameElem in f]
    if toLoad:
        l1_mapping = [f for f in l1_mapping if toLoad in f]

    # check input dataframe
    if isinstance(df, pd.DataFrame):
        namesInUse = df.columns.tolist()
        df_output = df.copy()
    else:
        if df == None:
            df_output = pd.DataFrame()
            namesInUse = []
        else:
            exit("df is not a DataFrame")

    for f in l1_mapping:
        tempDF = pd.read_csv(path + f, sep='\t', names=['name', 'value'])
        tempDF = tempDF.set_index('name').dropna()

        # clear names
        name = f.replace(nameElem, '')
        for c in toClear:
            name = name.replace(c, '')
        name = name + toAdd

        # overwrite warninig
        if name in namesInUse:
            if overwrite == False:
                return print(name + " exits in input df. Use overwrite=True to ignore.")

        # adding to dataframe
        df_output = pd.concat([df_output, tempDF['value'].rename(name)],axis=1)

    return df_output.reindex(sorted(df_output.columns), axis=1)


def read_STARstats(path="", toClear=[], toAdd="", df=None, overwrite=False):
    '''Reads multiple STAR Log final output files to one DataFrame

    :param path: Path to directory with files
    :type path: str
    :param toClear: List of strings to be removed from file names
    :type toClear: list
    :param toAdd: String to be added to file names
    :type toAdd: str
    :param df: DataFrame to be appended; default=None
    :type df: pandas.DataFrame, optional
    :param overwrite: Boolean, allows for overwriting during appending, default=False
    :type overwrite: bool

    :return: DataFrame with all files, where columns are values from 'name' column
    :rtype: pandas.DataFrame
    '''

    return read_tabFile(nameElem='_STARLog.final.out', path=path,
                        toClear=toClear, toAdd=toAdd, df=df, overwrite=overwrite)


def read_HTSeq_output(path="", toLoad="classes", toClear=[], toAdd="", df=None, overwrite=False):
    '''Reads multiple HTSeq tab files to one DataFrame

    :param path: Path to directory with files
    :type path: str
    :param toClear: List of strings to be removed from file names
    :type toClear: list
    :param toAdd: String to be added to file names
    :type toAdd: str
    :param df: DataFrame to be appended; default=None
    :type df: pandas.DataFrame, optional
    :param overwrite: Boolean, allows for overwriting during appending, default=False
    :type overwrite: bool

    :return: DataFrame with all files, where columns are values from 'name' column
    :rtype: pandas.DataFrame
    '''

    return read_tabFile(nameElem='_hittable.tab', path=path, toLoad=toLoad,
                        toClear=toClear, toAdd=toAdd, df=df, overwrite=overwrite)

def readSalmon(nameElem="", path="", toLoad="", toClear=[], toAdd="", column='NumReads', df=None, overwrite=False):
    '''Reads multiple Salmon quant.sf files to one DataFrame

    :param nameElem: Element to load, present in all files
    :type nameElem: str
    :param path: Path to directory with files
    :type path: str
    :param toLoad: Additional parameter for filtering, defaults to nameElem
    :type toLoad: str, optional
    :param toClear: List of strings to be removed from file names
    :type toClear: list
    :param toAdd: String to be added to file names
    :type toAdd: str
    :param column: Column to extract from quant.sf file, defaults to 'NumReads'
    :type column: str, optional
    :param df: DataFrame to be appended; default=None
    :type df: pandas.DataFrame, optional
    :param overwrite: Boolean, allows for overwriting during appending, default=False
    :type overwrite: bool

    :return: DataFrame with all files, where columns are values from 'Name' column
    :rtype: pandas.DataFrame
    '''

    # list files with STAT mapping
    l1_mapping = [f for f in os.listdir(path) if nameElem in f]
    if toLoad:
        l1_mapping = [f for f in l1_mapping if toLoad in f]

    # check input dataframe
    if isinstance(df, pd.DataFrame):
        namesInUse = df.columns.tolist()
        df_output = df.copy()
    else:
        if df == None:
            df_output = pd.DataFrame()
            namesInUse = []
        else:
            exit("df is not a DataFrame")

    for f in l1_mapping:
        tempDF = pd.read_csv(path + f+"/quant.sf", sep='\t')
        tempDF = tempDF.set_index('Name').dropna()

        # clear names
        name = f.replace(nameElem, '')
        for c in toClear:
            name = name.replace(c, '')
        name = name + toAdd

        # overwrite warninig
        if name in namesInUse:
            if overwrite == False:
                return print(name + " exits in input df. Use overwrite=True to ignore.")

        # adding to dataframe
        df_output = pd.concat([df_output, tempDF[column].rename(name)],axis=1)

    return df_output.reindex(sorted(df_output.columns), axis=1)

def loadGTF(gtf_path=""):
    '''Load GTF file into a DataFrame

    :param gtf_path: Path to the GTF file, defaults to ""
    :type gtf_path: str, optional

    :return: DataFrame containing GTF data
    :rtype: pandas.DataFrame
    '''

    names = ['chr','source','type','start','end','score','strand','phase','attributes']
    df_GTF_ncRNA = pd.read_csv(gtf_path,sep="\t",names=names,index_col=False)
    return df_GTF_ncRNA

def read_featureCount(nameElem="", path="", toLoad="", toClear=[], toAdd="", df=None, overwrite=False):
    '''Read featureCount files with common first column

    :param nameElem: Element to load, present in all files
    :type nameElem: str
    :param path: Path to directory with files
    :type path: str
    :param toLoad: Additional parameter for filtering, defaults to nameElem
    :type toLoad: str, optional
    :param toClear: List of strings to be removed from file names
    :type toClear: list
    :param toAdd: String to be added to file names
    :type toAdd: str
    :param df: DataFrame to be appended; default=None
    :type df: pandas.DataFrame, optional
    :param overwrite: Boolean, allows for overwriting during appending, default=False
    :type overwrite: bool

    :return: DataFrame with all files, where columns are values from 'name' column
    :rtype: pandas.DataFrame
    '''

    # list files
    l1_mapping = [f for f in os.listdir(path) if nameElem in f and '.summary' not in f]
    if toLoad:
        l1_mapping = [f for f in l1_mapping if toLoad in f]

    # check input dataframe
    if isinstance(df, pd.DataFrame):
        namesInUse = df.columns.tolist()
        df_output = df.copy()
    else:
        if df == None:
            df_output = pd.DataFrame()
            namesInUse = []
        else:
            exit("df is not a DataFrame")   
            
    for f in l1_mapping:
        tempDF = pd.read_csv(path + f, sep='\t', index_col=0, header=0, comment="#")

        # clear names
        name = f.replace(nameElem, '')
        for c in toClear:
            name = name.replace(c, '')
        name = name + toAdd

        # overwrite warninig
        if name in namesInUse:
            if overwrite == False:
                return print(name + " exits in input df. Use overwrite=True to ignore.")

        # adding to dataframe
        last_name = tempDF.columns.values[-1:][0]
        df_output = pd.concat([df_output, tempDF[last_name].rename(name)],axis=1)

    return df_output.reindex(sorted(df_output.columns), axis=1)

def read_DEseq(p):
    '''WARNING: Not tested properly. May not work as expected.
    Read DESeq2 output file and add gene names.

    :param p: Path to the DESeq2 output file
    :type p: str

    :return: DataFrame with DESeq2 results and gene names
    :rtype: pandas.DataFrame
    '''

    df = pd.read_csv(p, header=0, index_col=0)
    #adding gene name
    df.index = df.index.str.split(".").str[0]
    df["gene_name"] = df00_genename['Gene name']
    df["gene_name"] = df["gene_name"].astype(str)
    df["gene_name"][df["gene_name"]=='nan'] = df[df["gene_name"] =='nan'].index.tolist()
    
    return df[~df['padj'].isnull()]

def enriched(df, padj=0.05, fc=2):
    '''Filter DataFrame for enriched genes based on adjusted p-value and log2 fold change.

    :param df: DataFrame containing gene expression data
    :type df: pandas.DataFrame
    :param padj: Adjusted p-value threshold, defaults to 0.05
    :type padj: float, optional
    :param fc: Log2 fold change threshold, defaults to 2
    :type fc: float, optional

    :return: Filtered DataFrame with enriched genes
    :rtype: pandas.DataFrame
    '''

    return df[(df['padj'] < padj) & (df['log2FoldChange'] > fc)]

def bed2len(bed=pd.DataFrame()):
    '''Convert a bed DataFrame to a length DataFrame.

    :param bed: pandas DataFrame containing bed data.
    :type bed: pandas.DataFrame

    :return: pandas DataFrame with length data.
    :rtype: pandas.DataFrame

    :example:

    >>> bed = pd.DataFrame({'chr': ['chr1', 'chr1'], 'start': [100, 200], 'stop': [150, 250], 'region': ['region1', 'region2']})
    >>> bed2len(bed)
    region
    region1    50
    region2    50
    dtype: int64
    '''

    bed = bed.T[:4].T
    bed.columns = ['chr','start','stop','region']
    bed = bed.set_index('region')
    return bed['stop']-bed['start']

################################################
#############        handling multiple experiments

def define_experiments(paths_in, whole_name=False, strip='_hittable_reads.txt'):
    '''Parse file names and extract experiment name from them

    :param paths_in: List of file paths
    :type paths_in: list
    :param whole_name: Whether to use the whole file name as the experiment name, defaults to False
    :type whole_name: bool, optional
    :param strip: String to strip from the file name, defaults to '_hittable_reads.txt'
    :type strip: str, optional

    :return: List of experiment names and list of paths
    :rtype: tuple
    '''

    experiments = list()
    paths = list()
    for path in paths_in:
        paths.append(path.strip())
        file_path = path.split('/')
        file_name = file_path[len(file_path)-1]
        if whole_name == False: name = "_".join(file_name.split('_')[0:3]) #take fist three elements of file name as experiment name
        elif whole_name == True: name = file_name.rstrip(strip)
        experiments.append(name)
        if len(experiments) != len(paths):
            exit("No. of experiments is not equal to no. of paths")
    return experiments, paths

def expNameParser(name, additional_tags=list(), order='b_d_e_p'):
    '''Function handles experiment name; recognizes AB123456 as experiment date; BY4741 or HTP or given string as bait protein

    :param name: Experiment name
    :type name: str
    :param additional_tags: List of additional tags
    :type additional_tags: list
    :param order: Order of elements in the output, defaults to 'b_d_e_p'
    :type order: str, optional

    :return: Reordered experiment name
    :rtype: str
    '''

    tag_list = ['HTP', 'HTG', 'HTF', 'BY4741'] + additional_tags
    output_dict = {'b': str(), 'd': str(), 'e': str(), 'p': list()}  # bait; details; experiment; prefix
    name_elements = name.split('_')
    for e in name_elements:
        tag_in_e = [tag for tag in tag_list if tag in e]
        if tag_in_e and len(tag_in_e) >= 1:
            output_dict['b'] = e  # bait
            try:
                output_dict['d'] = name.split(tag_in_e[0], 1)[1].strip('_')  # details
            except:
                output_dict['d'] = 'wt'
                print('WARNING: wt added for '+name)
        elif re.search(r"[a-zA-Z][a-zA-Z]\d{6}", e) or re.search(r"[a-zA-Z][a-zA-Z][a-zA-Z]\d{6}", e):
            output_dict['e'] = e  # experiment name
            try:
                output_dict['p'] = name.split(e, 1)[0].strip('_')  # prefix
            except:
                output_dict['p'] = ''

    if len(output_dict['b']) < 3 or len(output_dict['e']) < 8:
        print(output_dict)
        sys.exit("ERROR: Can not define experiment or bait for "+name)

    return_list = list()
    for out in order.split('_'):
        return_list.append(output_dict[out])

    return '_'.join(return_list).strip('_')


def cleanNames(data, strings=[]):
    '''Cleans the names in the given data by removing specified strings.

    :param data: The data to be cleaned. It can be either a dictionary or a pandas DataFrame.
    :type data: dict or pandas.DataFrame
    :param strings: A list of strings to be removed from the names. Defaults to an empty list.
    :type strings: list, optional

    :return: The cleaned data with names modified according to the specified strings.
    :rtype: dict or pandas.DataFrame
    '''
    
    if isinstance(data, dict):
        for string in strings:
            data = {key.replace(string, ''): value for key, value in data.items()}
        return data
    elif isinstance(data, pd.DataFrame):
        for string in strings:
            data = data.rename(columns=lambda x: x.replace(string, ''))
        return data
    else:
        return data

     
def indexOrder(df=pd.DataFrame(), additional_tags=list(), output='root', order='b_d_e_p'):
    '''Apply expNameParser to whole DataFrame

    :param df: DataFrame where names of columns are names of experiments
    :type df: pandas.DataFrame
    :param additional_tags: List of additional tags to consider in expNameParser
    :type additional_tags: list
    :param output: Not used in the function, kept for compatibility
    :type output: str
    :param order: Order of elements in the output, defaults to 'b_d_e_p'
    :type order: str, optional

    :return: DataFrame with new names
    :rtype: pandas.DataFrame
    '''

    df = cleanNames(df, additional_tags=additional_tags)
    df.columns = [expNameParser(f, additional_tags=additional_tags, order=order) for f in list(df.columns.values)]
    return df.sort_index(axis=1)

def filterExp(datasets, let_in=[''], let_out=['wont_find_this_string'],verbose=False):
    '''Returns object with filtered columns/keys.

    :param datasets: DataFrame or dict with exp name as a key
    :param let_in: List of elements of name to filter in
    :type let_in: list
    :param let_out: List of elements of name to filter out
    :type let_out: list
    :param verbose: If True, prints the filtered columns/keys, defaults to False
    :type verbose: bool, optional

    :return: DataFrame or dict with filtered columns/keys
    :rtype: pandas.DataFrame or dict
    '''

    #for DataFrame()
    if isinstance(datasets, pd.DataFrame):
        output_df = pd.DataFrame()
        for f in [d for d in list(datasets.columns.values) if all(i in d for i in let_in) and all(o not in d for o in let_out)]:
            output_df[f]=datasets[f]
        if verbose: print(output_df.columns.tolist())
        return output_df

    #for dict()
    elif isinstance(datasets, dict):
        output_dict = dict()
        for f in [d for d in list(datasets.keys()) if all(i in d for i in let_in) and all(o not in d for o in let_out)]:
            output_dict[f]=datasets[f]
        if verbose: print(list(output_dict.keys()))
        return output_dict


def parseCRACname(s1=pd.Series):
    '''Parse CRAC name into ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3'] using this order.
     "_" is used to split the name

    :param s1: Series containing CRAC names
    :type s1: pandas.Series

    :return: DataFrame with parsed CRAC names
    :rtype: pandas.DataFrame
    '''

    df = s1.str.split("_", expand=True)
    df.columns = ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3']
    df['expFull'] = df['expID'] + "_" + df['expDate']
    df['sampleRep'] = None
    df['sample'] = None
    for i, row in df.iterrows():
        if row['condition3'] == None:
            n1 = row['protein'] + "_" + row['condition1'] + "_" + row['condition2']
        else:
            n1 = row['protein'] + "_" + row['condition1'] + "_" + row['condition2'] + "_" + row['condition3']

        if row['condition3'] == None:
            n2 = row['protein'] + "_" + row['condition1']
        else:
            n2 = row['protein'] + "_" + row['condition1'] + "_" + row['condition2']

        df.loc[i]['sampleRep'] = n1
        df.loc[i]['sample'] = n2
    return df


def groupCRACsamples(df=pd.DataFrame, use='protein', toDrop=[]):
    '''Parse CRAC names and annotates them using one of the following features:
    ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3', 'sample', 'sampleRep']

    :param df: DataFrame with CRAC names as index
    :type df: pandas.DataFrame
    :param use: Feature to annotate, choose from ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3', 'sample', 'sampleRep'], defaults to 'protein'
    :type use: str, optional
    :param toDrop: List of words in CRAC name that will qualify the sample for rejection, defaults to []
    :type toDrop: list, optional

    :return: DataFrame with added column ['group']
    :rtype: pandas.DataFrame
    '''

    df2 = parseCRACname(df.index.to_series())

    df2['group'] = None
    for group_name, df_temp in df2.groupby(use):
        df2['group'].loc[df_temp.index.tolist()] = group_name

    df['group'] = df2['group']

    if toDrop:
        dropping = [i for i in df.index.tolist() if any(d in i for d in toDrop)]
        return df.drop(dropping)
    else:
        return df

def expStats(input_df=pd.DataFrame(), smooth=True, window=10, win_type='blackman'):
    '''Returns DataFrame with 'mean', 'median', 'min', 'max' and quartiles if more than 2 experiments

    :param input_df: DataFrame containing the input data
    :type input_df: pandas.DataFrame
    :param smooth: Whether to apply smoothing window, defaults to True
    :type smooth: bool, optional
    :param window: Smoothing window size, defaults to 10
    :type window: int, optional
    :param win_type: Type of smoothing window, defaults to 'blackman'
    :type win_type: str, optional

    :return: DataFrame with calculated statistics
    :rtype: pandas.DataFrame
    '''

    working_df, result_df = pd.DataFrame(), pd.DataFrame()

    #smoothing
    if smooth == True:
        for f in input_df.columns.values:
            print(f)
            working_df[f]=input_df[f].rolling(window, win_type=win_type, center=True).mean()
    else:
        working_df = input_df.copy()

    #calculating stats
    for function in ['mean', 'median', 'min', 'max']: result_df[function]=getattr(working_df, function)(axis=1) #calculates using pandas function listed in []
    if len(working_df.columns) > 2: #calculating quartiles only in more than two experiments
        result_df['q1'], result_df['q3'] = working_df.quantile(q=0.25, axis=1), working_df.quantile(q=0.75, axis=1)

    return result_df

################################################
#############       statistics and analysis

def normalize(df=pd.DataFrame, log2=False, pseudocounts=0.1, CPM=True):
    '''Normalize the DataFrame by adding pseudocounts, dividing by the sum, and optionally applying log2 transformation and CPM scaling.

    :param df: DataFrame to be normalized
    :type df: pandas.DataFrame
    :param log2: Whether to apply log2 transformation, defaults to False
    :type log2: bool, optional
    :param pseudocounts: Value to add to each element to avoid log of zero, defaults to 0.1
    :type pseudocounts: float, optional
    :param CPM: Whether to scale the data to counts per million, defaults to True
    :type CPM: bool, optional

    :return: Normalized DataFrame
    :rtype: pandas.DataFrame
    '''

    df = df.add(pseudocounts)
    df = df / df.sum()
    if CPM==True:
        df = df.multiply(1000000)
    
    if log2==True:
        return df.apply(np.log2)
    else:
        return df

def quantileCategory(s1=pd.Series(dtype=float), q=4):
    '''Quantile-based discretization function based on pandas.qcut function.

    :param s1: Series to be discretized
    :type s1: pandas.Series
    :param q: Number of quantiles, defaults to 4
    :type q: int, optional

    :return: Series with quantile categories
    :rtype: pandas.Series

    :example:

    >>> quantileCategory(pd.Series([1, 2, 3, 4, 5]), q=2)
    0       0
    1       0
    2       1
    3       1
    4       1
    dtype: int64
    '''

    temp_df = pd.DataFrame()
    temp_df['data'] = s1
    temp_df['quantiles'] = 200
    quantiles = pd.qcut(s1, q, retbins=True)
    for i in range(0, len(quantiles[1]) - 1):
        temp_df['quantiles'][temp_df.data >= quantiles[1][i]] = i
    return temp_df['quantiles']

def runPCA(data=pd.DataFrame(), n_components=2):
    '''Run PCA analysis and re-assigns column names and index names

    :param data: DataFrame containing the input data
    :type data: pandas.DataFrame
    :param n_components: Number of principal components to compute, defaults to 2
    :type n_components: int, optional

    :return: Tuple consisting of DataFrame with PCA results and a list of explained variance ratios
    :rtype: tuple

    :example:

    >>> runPCA(pd.DataFrame([[1, 2], [3, 4], [5, 6]]), n_components=2)
    (          PC1       PC2
    0 -4.242641  0.000000
    1  0.000000  0.000000
    2  4.242641  0.000000, [100.0, 0.0])
    '''

    # x = StandardScaler().fit_transform(df1_codone_composition)

    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(data)
    principalDf = pd.DataFrame(data=principalComponents,
                               columns=['PC' + str(i + 1) for i in np.arange(0, n_components)])

    values = pca.explained_variance_ratio_

    finalDf = pd.concat([principalDf, pd.Series(data.index)], axis=1)
    finalDf.rename(columns={0: 'name'}, inplace=True)
    finalDf = finalDf.set_index('name')

    return finalDf, [round(i*100,2) for i in values.tolist()]

def addCluster(df=pd.DataFrame(), n=10):
    '''Assigns n clusters to the data using KMeans algorithm

    :param df: DataFrame containing the data to be clustered
    :type df: pandas.DataFrame
    :param n: Number of clusters to form, defaults to 10
    :type n: int, optional

    :return: DataFrame with an additional 'cluster' column indicating cluster assignment
    :rtype: pandas.DataFrame

    :example:

    >>> df = pd.DataFrame({'x': [1, 2, 3, 4, 5], 'y': [5, 4, 3, 2, 1]})
    >>> addCluster(df, n=2)
       x  y  cluster
    0  1  5        1
    1  2  4        1
    2  3  3        0
    3  4  2        0
    4  5  1        0
    '''

    if 'cluster' in df.columns:
        df = df.drop('cluster', 1)
    kmeans = KMeans(n_clusters=n, random_state=0).fit(df)
    df['cluster'] = kmeans.labels_
    df = df.sort_values('cluster', ascending=True)

    # summary
    tt.plotting.clusterClusterMap(df)

    return df

################################################
#############        other

def timestamp():
    '''Returns current timestamp as a string

    :return: timestamp in a format 'YYYYMMDD_HHMMSS'
    :rtype: str
    '''

    return str(time.strftime("%Y%m%d_%H%M%S"))

def timestampRandomInt():
    '''Returns current timestamp with a random integer as a string

    :return: timestamp in a format 'YYYYMMDD_HHMMSS_RANDOMINT'
    :rtype: str
    '''
    
    return str(time.strftime("%Y%m%d_%H%M%S"))+"_"+str(random.randint(0,1000))

###         OLD         ###
# def getRefFile(file_from_options, file_type):
#     """
#     Sorting out source of gtf, fasta or tab path, in order (1) from options parser, (2) from ~/bin/default.aml
#      or (3) from environmental variable $xxx_PATH
#     :param file_from_options: file path from options parser, can be an empty string
#                 file_type: 'GTF', 'FASTA', 'TAB'
#     :return: path to the file
#     """
#     file_types = {'GTF' : 'GTF_PATH', 'FASTA' : 'FASTA_PATH', 'TAB' : 'TAB_PATH'}
#     if file_from_options:
#         return file_from_options
#     elif 'default.aml' in os.listdir(os.getenv("HOME")+'/bin/'):
#         default = yaml.load(open(os.getenv("HOME")+'/bin/default.aml'))
#         # print "# Using "+file_type+" file from ~/bin/default.aml"
#         return default[file_types[file_type]]
#     else:
#         if file_types[file_type] in os.environ:
#             # print "# Using "+file_type+" file from $GTF_PATH variable"
#             return os.environ[file_types[file_type]]
#         else:
#             exit('Provide '+file_type+' file path using -g or setup default.aml in your bin folder')
#
# def getGTF(gtf_from_options):
#     return getRefFile(gtf_from_options, 'GTF')
#
# def getFASTA(fasta_from_options):
#     return getRefFile(fasta_from_options, 'FASTA')
#
# def getTAB(tab_from_options):
#     return getRefFile(tab_from_options, 'TAB')
