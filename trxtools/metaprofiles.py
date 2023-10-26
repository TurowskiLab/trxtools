import pandas as pd
import pyBigWig

### level -3
def read_bed(bed_path):
    bed_df = pd.read_csv(bed_path, sep='\t', header=None)
    if len(bed_df.columns) > 6:
        bed_df.drop([6,7,8,9,10,11], axis=1, inplace=True)
    # Set score column to 0
    bed_df[4] = 0
    return bed_df

### level -2
def get_bw_data(bed_row, bw, flank_5=0, flank_3=0, align_3end=False):
    '''
    Retrieve BigWig scores for positions in a given region, optionally including flanks of given length.
    '''
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
    '''
    Split a BED dataframe by strand.
    
    :returns: bed_plus, bed_minus
    :rtype: pandas.DataFrame
    '''
    bed_plus = bed_df[bed_df[5] == "+"]
    bed_minus = bed_df[bed_df[5] == "-"]
    return bed_plus, bed_minus

def matrix_from_bigwig(bw_path, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    '''
    Get matrix with bigwig scores for all regions in a bed df.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.
    '''
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
    '''
    Combine score matrices for regions on + and - strands.
    '''
    out_dict={}
    for key in plus_dict.keys():
        key_min = key.replace("plus", "minus")
        key_min = key.replace("fwd", "rev")
        if key_min not in minus_dict.keys():
            raise Exception("Keys in dictionaries do not match")
            break
        else:
            out_dict[key] = pd.concat([plus_dict[key], minus_dict[key_min]], ignore_index=True)
    return out_dict

### level 0
def get_multiple_matrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    """Get score matrices for positions in given regions (with optional flanks) from multiple BigWig files.
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
        plus_dict[bw] = matrix_from_bigwig(bw_path=bw, bed_df=bed_plus, flank_5=flank_5, flank_3=flank_3, fill_na=fill_na, pseudocounts=pseudocounts, normalize_libsize=normalize_libsize, align_3end=align_3end)
    minus_dict = {}
    for bw in bw_paths_minus:
        minus_dict[bw] = matrix_from_bigwig(bw_path=bw, bed_df=bed_minus, flank_5=flank_5, flank_3=flank_3, fill_na=fill_na, pseudocounts=pseudocounts, normalize_libsize=normalize_libsize, align_3end=align_3end)
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