from trxtools.utils import bash
from trxtools.utils import files
from trxtools.utils import names
from trxtools.utils import stats
from trxtools.utils import sequences
import warnings
import pandas as pd

### METHODS ###
# various python methods

################################################
#############        bash

def listPaths(folder=".", suffix=str(), nameElem=None):
    '''List all file paths in a folder with a specific suffix and optional name element.
    
    Warning: This function is deprecated. Please use bash.listPaths() instead.

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

    
    warnings.warn(
        "listPaths() is deprecated. Please use utils.bash.listPaths() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    return bash.listPaths(folder=folder, suffix=suffix, nameElem=nameElem)


def bashCommand(bashCommand=str()):
    '''Run command in bash using subprocess.call()
    
    :param bashCommand: str with command

    :return: None
    '''
    
    warnings.warn(
        "bashCommand() is deprecated. Please use utils.bash.bashCommand() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return bash.bashCommand(bashCommand=bashCommand)

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
    warnings.warn(
        "calGC() is deprecated. Please use utils.sequences.calGC() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    return sequences.calGC(dataset=dataset, calFor=calFor)

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
    warnings.warn(
        "reverse_complement_DNA() is deprecated. Please use utils.sequences.reverse_complement_DNA() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    return sequences.reverse_complement_DNA(seq=seq)


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
    warnings.warn(
        "reverse_complement_RNA() is deprecated. Please use utils.sequences.reverse_complement_RNA() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.reverse_complement_RNA(seq=seq)

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
    warnings.warn(
        "reverse_complement() is deprecated. Please use utils.sequences.reverse_complement() instead.",
        DeprecationWarning,
        stacklevel=2
    )   
    return sequences.reverse_complement(seq=seq)

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
    warnings.warn(
        "randomDNAsingle() is deprecated. Please use utils.sequences.randomDNAsingle() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.randomDNAsingle(length=length, letters=letters)

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
    warnings.warn(
        "randomDNAall() is deprecated. Please use utils.sequences.randomDNAall() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.randomDNAall(length=length, letters=letters)


def rollingGC(s=pd.Series, window=10): #rolling window, smoothing of data
    '''Calculates GC from sequence, uses 'boxcar' window

    :param s: Series containing sequence
    :type s: pandas.Series
    :param window: window size for GC calculation
    :type window: int

    :return: Series with GC calculated, center=False
    :rtype: pandas.Series
    '''
    warnings.warn(
        "rollingGC() is deprecated. Please use utils.sequences.rollingGC() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.rollingGC(s=s, window=window)

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
    warnings.warn(
        "letterContent() is deprecated. Please use utils.sequences.letterContent() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.letterContent(seq=seq, letter=letter)


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

    warnings.warn(
        "DNA_string_positions() is deprecated. Please use utils.sequences.DNA_string_positions() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.DNA_string_positions(df=df, string=string, name_col=name_col, seq_col=seq_col)

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
    warnings.warn(
        "is_inside() is deprecated. Please use utils.sequences.is_inside() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.is_inside(inner_start, inner_end, outer_start, outer_end)


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

    warnings.warn(
        "nested_region_cleanup() is deprecated. Please use utils.sequences.nested_region_cleanup() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.nested_region_cleanup(df=df, start_col=start_col, end_col=end_col)
        

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
    warnings.warn(
        "DNA_stretch_positions() is deprecated. Please use utils.sequences.DNA_stretch_positions() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return sequences.DNA_stretch_positions(df=df, char=char, min_len=min_len, max_len=max_len, name_col=name_col, seq_col=seq_col)


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

    warnings.warn(
        "find_pol3_terminators() is deprecated. Please use utils.sequences.find_pol3_terminators() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    return sequences.find_pol3_terminators(df=df, min_T=min_T, max_T=max_T, name_col=name_col, seq_col=seq_col)

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
    warnings.warn(
        "read_list() is deprecated. Please use utils.files.read_list() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.read_list(filepath=filepath)

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

    warnings.warn(
        "read_tabFile() is deprecated. Please use utils.files.read_tabFile() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.read_tabFile(nameElem=nameElem, path=path, toLoad=toLoad,
                              toClear=toClear, toAdd=toAdd, df=df, overwrite=overwrite)


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
    warnings.warn(
        "read_STARstats() is deprecated. Please use utils.files.read_STARstats() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.read_STARstats(path=path, toClear=toClear,
                             toAdd=toAdd, df=df, overwrite=overwrite)

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
    warnings.warn(
        "read_HTSeq_output() is deprecated. Please use utils.files.read_HTSeq_output() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.read_HTSeq_output(path=path, toLoad=toLoad,
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

    warnings.warn(
        "readSalmon() is deprecated. Please use utils.files.readSalmon() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.readSalmon(nameElem=nameElem, path=path, toLoad=toLoad,
                            toClear=toClear, toAdd=toAdd, column=column, df=df, overwrite=overwrite)

def loadGTF(gtf_path=""):
    '''Load GTF file into a DataFrame

    :param gtf_path: Path to the GTF file, defaults to ""
    :type gtf_path: str, optional

    :return: DataFrame containing GTF data
    :rtype: pandas.DataFrame
    '''

    warnings.warn(
        "loadGTF() is deprecated. Please use utils.files.loadGTF() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.loadGTF(gtf_path=gtf_path)


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
    warnings.warn(
        "read_featureCount() is deprecated. Please use utils.files.read_featureCount() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.read_featureCount(nameElem=nameElem, path=path, toLoad=toLoad,
                              toClear=toClear, toAdd=toAdd, df=df, overwrite=overwrite)

def read_DEseq(p):
    '''WARNING: Not tested properly. May not work as expected.
    Read DESeq2 output file and add gene names.

    :param p: Path to the DESeq2 output file
    :type p: str

    :return: DataFrame with DESeq2 results and gene names
    :rtype: pandas.DataFrame
    '''
    warnings.warn(
        "read_DEseq() is deprecated. Please use utils.files.read_DEseq() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.read_DEseq(p=p)


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

    warnings.warn(
        "bed2len() is deprecated. Please use utils.files.bed2len() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.bed2len(bed=bed)

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

    warnings.warn(
        "define_experiments() is deprecated. Please use utils.names.define_experiments() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.define_experiments(paths_in=paths_in, whole_name=whole_name, strip=strip)


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

    warnings.warn(
        "expNameParser() is deprecated. Please use utils.names.expNameParser() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.expNameParser(name=name, additional_tags=additional_tags, order=order)

def cleanNames(data, strings=[]):
    '''Cleans the names in the given data by removing specified strings.

    :param data: The data to be cleaned. It can be either a dictionary or a pandas DataFrame.
    :type data: dict or pandas.DataFrame
    :param strings: A list of strings to be removed from the names. Defaults to an empty list.
    :type strings: list, optional

    :return: The cleaned data with names modified according to the specified strings.
    :rtype: dict or pandas.DataFrame
    '''
    warnings.warn(
        "cleanNames() is deprecated. Please use utils.names.cleanNames() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.cleanNames(data=data, strings=strings)

     
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

    warnings.warn(
        "indexOrder() is deprecated. Please use utils.names.indexOrder() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.indexOrder(df=df, additional_tags=additional_tags, output=output, order=order)


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

    warnings.warn(
        "filterExp() is deprecated. Please use utils.names.filterExp() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.filterExp(datasets=datasets, let_in=let_in, let_out=let_out, verbose=verbose)



def parseCRACname(s1=pd.Series):
    '''Parse CRAC name into ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3'] using this order.
     "_" is used to split the name

    :param s1: Series containing CRAC names
    :type s1: pandas.Series

    :return: DataFrame with parsed CRAC names
    :rtype: pandas.DataFrame
    '''

    warnings.warn(
        "parseCRACname() is deprecated. Please use utils.names.parseCRACname() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.parseCRACname(s1=s1)


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

    warnings.warn(
        "groupCRACsamples() is deprecated. Please use utils.names.groupCRACsamples() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return names.groupCRACsamples(df=df, use=use, toDrop=toDrop)


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

    warnings.warn(
        "expStats() is deprecated. Please use utils.stats.expStats() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return stats.expStats(input_df=input_df, smooth=smooth, window=window, win_type=win_type)


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

    warnings.warn(
        "normalize() is deprecated. Please use utils.stats.normalize() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return stats.normalize(df=df, log2=log2, pseudocounts=pseudocounts, CPM=CPM)

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

    warnings.warn(
        "quantileCategory() is deprecated. Please use utils.stats.quantileCategory() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return stats.quantileCategory(s1=s1, q=q)

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

    warnings.warn(
        "runPCA() is deprecated. Please use utils.stats.runPCA() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return stats.runPCA(data=data, n_components=n_components)


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

    warnings.warn(
        "addCluster() is deprecated. Please use utils.stats.addCluster() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return stats.addCluster(df=df, n=n)

################################################
#############        other

def timestamp():
    '''Returns current timestamp as a string

    :return: timestamp in a format 'YYYYMMDD_HHMMSS'
    :rtype: str
    '''
    warnings.warn(
        "timestamp() is deprecated. Please use utils.files.timestamp() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.timestamp()

def timestampRandomInt():
    '''Returns current timestamp with a random integer as a string

    :return: timestamp in a format 'YYYYMMDD_HHMMSS_RANDOMINT'
    :rtype: str
    '''
    warnings.warn(
        "timestampRandomInt() is deprecated. Please use utils.files.timestampRandomInt() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return files.timestampRandomInt()