import pandas as pd
import pyBigWig
import warnings

### ALIAS FUNCTIONS
### these will be removed in a later release!

def read_bed(bed_path):
    warnings.warn(
        'read_bed() will be renamed to readBED() in a future release. Update your code to silence this warning.',
        FutureWarning
        )
    return readBED(bed_path)

def matrix_from_bigwig(bw_path, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    warnings.warn(
        'matrix_from_bigwig() will be renamed to matrixFromBigWig() in a future release. Update your code to silence this warning.',
        FutureWarning
    )
    return matrixFromBigWig(
        bw_path=bw_path,
        bed_df=bed_df,
        flank_5=flank_5,
        flank_3=flank_3,
        fill_na=fill_na,
        pseudocounts=pseudocounts,
        normalize_libsize=normalize_libsize,
        align_3end=align_3end
        )
def get_multiple_matrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    warnings.warn(
        'get_multiple_matrices() will be renamed to getMultipleMatrices() in a future release. Update your code to silence this warning.',
        FutureWarning
    )
    return getMultipleMatrices(
        bw_paths_plus=bw_paths_plus,
        bw_paths_minus=bw_paths_minus,
        bed_df=bed_df,
        flank_5=flank_5,
        flank_3=flank_3,
        fill_na=fill_na,
        pseudocounts=pseudocounts,
        normalize_libsize=normalize_libsize,
        align_3end=align_3end
    )


### level -3
def readBED(bed_path):
    """
    Simple BED file parser

    :param bed_path: Path to BED file
    :type bed_path: str
    :return: Contents of BED file in DataFrame form
    :rtype: pandas.DataFrame
    """
    bed_df = pd.read_csv(bed_path, sep='\t', header=None)
    if len(bed_df.columns) > 6:
        bed_df.drop([6,7,8,9,10,11], axis=1, inplace=True)
    # Set score column to 0
    bed_df[4] = 0
    return bed_df

### level -2
def get_bw_data(bed_row, bw, flank_5=0, flank_3=0, align_3end=False):
    # Retrieve BigWig scores for positions in a given region, optionally including flanks of given length.
    if bed_row[5] == '+':
        outseries = pd.Series(bw.values(bed_row[0], bed_row[1]-flank_5, bed_row[2]+flank_3))
    if bed_row[5] == '-':
        # flip orientation on minus strand
        # keep the index so pandas doesn't flip it again
        outseries = pd.Series(bw.values(bed_row[0], bed_row[1]-flank_3, bed_row[2]+flank_5))
        rev_index = outseries.index
        outseries = outseries.iloc[::-1]
        outseries.index = rev_index
    if align_3end:
        outseries.index = outseries.index - flank_5 - (bed_row[2] - bed_row[1])
    else:
        outseries.index = outseries.index - flank_5
    return outseries

### level -1    
def bed_split_strands(bed_df):
    # Split dataframe (from read_bw) by strand
    bed_plus = bed_df[bed_df[5] == "+"]
    bed_minus = bed_df[bed_df[5] == "-"]
    return bed_plus, bed_minus

def matrixFromBigWig(bw_path, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    """
    Get matrix with BigWig scores for all regions in a bed df from a single BigWig file.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.
    Flanks are strand-aware.

    :param bw_path: Path to target BigWig file
    :type bw_path: str
    :param bed_df: DataFrame with BED data (use read_bed())
    :type bed_df: pandas.DataFrame
    :param flank_5: length of 5'flank to extend BED regions by, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true replace NaN values with 0 (pybigwig returns positions with 0 coverage as NaN), defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to retrieved scores, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: Whether to normalize output to library size, this is usually the right thing to do, defaults to True
    :type normalize_libsize: bool, optional
    :param align_3end: If true, position 0 in the resulting matrix will be set at the target region's 3'end instead of 5'end, defaults to False
    :type align_3end: bool, optional
    :return: DataFrame with the result score matrix
    :rtype: pandas.DataFrame
    """    
    bw = pyBigWig.open(bw_path)
    out_df = bed_df.apply(get_bw_data, bw=bw, flank_5=flank_5, flank_3=flank_3, align_3end=align_3end, axis=1)
    # if out_df.columns[0] != 0:
    #     out_df.columns = out_df.columns[::-1]
    out_df['region'] = bed_df[3]
    out_df = out_df.set_index('region')
    if fill_na:
        out_df = out_df.fillna(0)
    if isinstance(pseudocounts, float):
        out_df = out_df+pseudocounts
    if normalize_libsize:
        libsize = bw.header()['sumData']
        out_df = out_df/libsize
    bw.close()
    out_df = out_df.reset_index()
    return out_df

def join_strand_matrices(plus_dict, minus_dict):
    # Combine score matrices for regions on + and - strands.

    
    # /home/jm/repos/trxtools/trxtools/metaprofiles.py:143: FutureWarning:
    # The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. 
    # In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. 
    # To retain the old behavior, exclude the relevant entries before the concat operation.
    
    out_dict={}
    for key in plus_dict.keys():
        key_min = key.replace("plus", "minus")
        key_min = key.replace("fwd", "rev")
        if key_min not in minus_dict.keys():
            raise Exception("Keys in dictionaries do not match")
            break
        else:
            # This raises this FutureWarning:
            # /home/jm/repos/trxtools/trxtools/metaprofiles.py:143: FutureWarning:
            # The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. 
            # In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. 
            # To retain the old behavior, exclude the relevant entries before the concat operation.
            out_dict[key] = pd.concat([plus_dict[key], minus_dict[key_min]], ignore_index=True)
    return out_dict

### level 0
def getMultipleMatrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    """
    Get score matrices for positions in given regions (with optional flanks) from multiple BigWig files.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.

    :param bw_paths_plus: list of paths to BigWig files (+ strand)
    :type bw_paths_plus: list
    :param bw_paths_minus: list of paths to BigWig files (- strand)
    :type bw_paths_minus: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
    :type flank_5: int, optional
    :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
    :type flank_3: int, optional
    :param fill_na: _description_, defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to datapoints, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
    :type normalize_libsize: bool, optional
    :param align_3end: if True, align profiles at the 3'end of features. Otherwise align at 5'end, defaults to False
    :type align_3end: bool, optional
    :return:  A dictionary containing score matrices for individual BigWig files. Dictionary keys are BigWig file names.
    :rtype: dict
    """    
    bed_plus, bed_minus = bed_split_strands(bed_df)
    plus_dict = {}
    for bw in bw_paths_plus:
        plus_dict[bw] = matrixFromBigWig(bw_path=bw, bed_df=bed_plus, flank_5=flank_5, flank_3=flank_3, fill_na=fill_na, pseudocounts=pseudocounts, normalize_libsize=normalize_libsize, align_3end=align_3end)
    minus_dict = {}
    for bw in bw_paths_minus:
        minus_dict[bw] = matrixFromBigWig(bw_path=bw, bed_df=bed_minus, flank_5=flank_5, flank_3=flank_3, fill_na=fill_na, pseudocounts=pseudocounts, normalize_libsize=normalize_libsize, align_3end=align_3end)
    return join_strand_matrices(plus_dict, minus_dict)

def metaprofile(matrix_dict, agg_type='mean', normalize_internal=False):
    '''
    Calculate metaprofiles from score matrices by aggregating each position in all regions. These can then be plotted with your favorite lib.

    :param matrix_dict: Dict containing score matrices returned by get_multiple_matrices()
    :type matrix_dict: dict
    :param agg_type: Type of aggregation to use. Available: 'mean', 'median', 'sum', defaults to 'mean'
    :type agg_type: str, optional
    :param normalize: if true, normalize each profile internally (i.e. gene-wise) (x/sum(x)) before aggregating, defaults to False
    :type normalize: bool, optional
    :raises Exception: _description_
    :return: dataframe containing metaprofile values for each position in each of the input matrices (i.e. bigwig files)
    :rtype: pandas.DataFrame
    '''    
    if agg_type not in ['mean', 'median', 'sum']:
        raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")
    if normalize_internal:
        dropped = {key: value.drop('region', axis=1).div(value.sum(axis=1,numeric_only=True),axis=0) for key, value in matrix_dict.items()}
        return pd.DataFrame({key: value.agg(agg_type,numeric_only=True) for key, value in dropped.items()})
    else:
        return pd.DataFrame({key: value.agg(agg_type, numeric_only=True) for key, value in matrix_dict.items()})
    
def regionScore(bw_paths_plus, bw_paths_minus, bed_df, agg_type='sum', flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True):
    """
    Calculate coverage or statistic for multiple regions in multiple BigWig files.

    :param bw_paths_plus: list of paths to BigWig files (+ strand)
    :type bw_paths_plus: list
    :param bw_paths_minus: list of paths to BigWig files (- strand)
    :type bw_paths_minus: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param agg_type: operation to perform on region scores. Available options: 'sum', 'mean', 'median', defaults to 'sum'
    :type agg_type: str, optional
    :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
    :type flank_5: int, optional
    :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true, NaNs will be replaced with zeroes (recommended, as pyBigWig reports positions with 0 score as NaN), defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to datapoints, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
    :type normalize_libsize: bool, optional
    :return:  DataFrame with calculated scores. Rows are regions/genes, columns are BigWig files.
    :rtype: pandas.DataFrame
    """    
    if agg_type not in ['mean', 'median', 'sum']:
        raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")
    # Get score matrices
    outdict = getMultipleMatrices(
        bw_paths_plus=bw_paths_plus,
        bw_paths_minus=bw_paths_minus,
        bed_df=bed_df,
        flank_5=flank_5,
        flank_3=flank_3,
        fill_na=fill_na,
        pseudocounts=pseudocounts,
        normalize_libsize=normalize_libsize,
        align_3end=False
        )
    # Aggregate all positions per gene using chosen function
    outdict = {key:value.set_index('region').agg(func=agg_type, axis=1) for key, value in outdict.items()}
    out_df = pd.DataFrame(outdict)
    return out_df