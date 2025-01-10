import pandas as pd
import numpy as np
import pyBigWig
import warnings
import trxtools as tt
from pybedtools import BedTool

### ALIAS FUNCTIONS
### these will be removed in a later release!

def read_bed(bed_path):
    """Simple BED file parser

    :param bed_path: Path to BED file
    :type bed_path: str
    :return: Contents of BED file in DataFrame form
    :rtype: pandas.DataFrame

    .. warning:: This function is deprecated and will be removed in future versions.
    """
    warnings.warn(
        'read_bed() will be renamed to readBED() in a future release. Update your code to silence this warning.',
        FutureWarning
        )
    return readBED(bed_path)

def matrix_from_bigwig(bw_path, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, align_3end=False):
    '''WARNINIG: This function is replaced by matrixFromBigWig() \n 
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
    :param align_3end: If true, position 0 in the resulting matrix will be set at the target region's 3'end instead of 5'end, defaults to False
    :type align_3end: bool, optional

    :return: DataFrame with the result score matrix
    :rtype: pandas.DataFrame

    .. warning:: This function is deprecated and will be removed in future versions.
    '''
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
        align_3end=align_3end
        )

def get_multiple_matrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    '''WARNINIG: This function is replaced by getMultipleMatrices() \n
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

    .. warning:: This function is deprecated and will be removed in future versions.
    '''
    
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
    '''Simple BED file parser

    :param bed_path: Path to BED file
    :type bed_path: str
    :return: Contents of BED file in DataFrame form
    :rtype: pandas.DataFrame
    '''
    bed_df = pd.read_csv(bed_path, sep='\t', header=None)
    if len(bed_df.columns) > 6:
        bed_df.drop([6,7,8,9,10,11], axis=1, inplace=True)
    # Set score column to 0
    bed_df[4] = 0
    return bed_df

### level -2
def get_bw_data(bed_row, bw, flank_5=0, flank_3=0, align_3end=False):
    '''Retrieve BigWig scores for positions in a given region, optionally including flanks of given length.

    :param bed_row: A row from a BED file, containing chromosome, start, end, and strand information.
    :type bed_row: list or pandas.Series
    :param bw: A pyBigWig object to retrieve values from.
    :type bw: pyBigWig
    :param flank_5: Length of the 5' flank to include, defaults to 0.
    :type flank_5: int, optional
    :param flank_3: Length of the 3' flank to include, defaults to 0.
    :type flank_3: int, optional
    :param align_3end: Whether to align the series to the 3' end, defaults to False.
    :type align_3end: bool, optional
    :return: A pandas Series containing the BigWig scores for the specified region.
    :rtype: pandas.Series
    '''

    # Retrieve BigWig scores for positions in a given region, optionally including flanks of given length.
    if bed_row[5] == '+':
        outseries = pd.Series(bw.values(bed_row[0], int(bed_row[1]-flank_5), int(bed_row[2]+flank_3)))
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
    '''Splits a BED dataframe into two separate dataframes based on strand information.

    :param bed_df: A dataframe containing BED format data. The strand information is expected to be in the 6th column (index 5).
    :type bed_df: pandas.DataFrame
    :return: A tuple containing two dataframes:
        - bed_plus (pandas.DataFrame): Dataframe containing entries with the '+' strand.
        - bed_minus (pandas.DataFrame): Dataframe containing entries with the '-' strand.
    :rtype: tuple
    '''
    # Split dataframe (from read_bw) by strand
    bed_plus = bed_df[bed_df[5] == "+"]
    bed_minus = bed_df[bed_df[5] == "-"]
    return bed_plus, bed_minus

def matrixFromBigWig(bw_path, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, align_3end=False):
    '''Get matrix with BigWig scores for all regions in a bed df from a single BigWig file.
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
    :param align_3end: If true, position 0 in the resulting matrix will be set at the target region's 3'end instead of 5'end, defaults to False
    :type align_3end: bool, optional
    :return: DataFrame with the result score matrix
    :rtype: pandas.DataFrame
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
    bw.close()
    out_df = out_df.reset_index()
    return out_df

def join_strand_matrices(plus_dict, minus_dict):
    '''Combine score matrices for regions on + and - strands.

    :param plus_dict: Dictionary containing score matrices for regions on the + strand.
    :type plus_dict: dict
    :param minus_dict: Dictionary containing score matrices for regions on the - strand.
    :type minus_dict: dict
    :raises Exception: If keys in dictionaries do not match.
    :return: Dictionary containing combined score matrices for regions on both strands.
    :rtype: dict
    '''
    # Combine score matrices for regions on + and - strands.

    
    # /home/jm/repos/trxtools/trxtools/metaprofiles.py:143: FutureWarning:
    # The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. 
    # In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. 
    # To retain the old behavior, exclude the relevant entries before the concat operation.
    
    out_dict={}
    for key in list(plus_dict.keys()):
        key_min = key.replace("plus", "minus")
        key_min = key_min.replace("fwd", "rev")
        if key_min not in minus_dict.keys():
            raise Exception("Keys in dictionaries do not match: " + str(key_min))
        else:
            # This raises this FutureWarning:
            # /home/jm/repos/trxtools/trxtools/metaprofiles.py:143: FutureWarning:
            # The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. 
            # In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. 
            # To retain the old behavior, exclude the relevant entries before the concat operation.
            out_dict[key] = pd.concat([plus_dict[key], minus_dict[key_min]], ignore_index=True)
    return out_dict

def peak2matrice(bed_df=pd.DataFrame, peak_file_path='', 
                 g='references/GRCh38.primary_assembly.genome.cleaned.len', 
                 flank_5=0, flank_3=0, fillna=True):
    '''Get score matrices for positions in given regions (with optional flanks) from a single nroadPeak or narrowPeak file. 
    Matrix rows correspond to regions in BED. Columns correpond to nucleotide positions in regions + flanks. Flanks are strand-aware.
    
    :param bed_df: DataFrame with BED data (can use read_bed())
    :type bed_df: pandas.DataFrame
    :param peak_file_path: Path to peak file (broadPeak or narrowPeak format)
    :type peak_file_path: str
    :param g: Path to genome file (chromosome lengths), defaults to 'references/GRCh38.primary_assembly.genome.cleaned.len'
    :type g: str
    :param flank_5: length of 5'flank to extend BED regions by, strand-aware, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, strand-aware, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true replace NaN values with 0 (pybigwig returns positions with 0 coverage as NaN), defaults to True
    :type fill_na: bool, optional
    :return: DataFrame with the result score matrix
    :rtype: pandas.DataFrame
    '''
    
    peak_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pvalue', 'qValue','peak'] #broadPeak do not have 'peak' column
    peak_df = pd.read_csv(peak_file_path, sep='\t', header=None, names=peak_columns)
    a = BedTool.from_dataframe(peak_df)
    
    bed_df.index=bed_df.iloc[:, 3].rename('region')
    
    output_df = pd.DataFrame()
    for i,row in bed_df.iterrows():
        b = BedTool.from_dataframe(pd.DataFrame(row[:7]).T)
        b_extended = b.slop(s=True, l=flank_5, r=flank_3, g=g)
        
        peaks = a.intersect(b_extended).to_dataframe(names=peak_columns)
        
        start,end,strand = b_extended.to_dataframe()[['start','end','strand']].loc[0].tolist()
#         print(end-start)
        profile = pd.Series(index=np.arange(-flank_5,end-start-flank_5),data=0,name=i, dtype=float)

        ###code here how to convert selected beds into Series with the data
       
        for p,row in peaks.iterrows():
            p_start, p_end, p_strand = row['start'], row['end'], row['strand']
            # for "+" strand
            if strand == "+" and (p_strand == "." or p_strand == "+"):
                profile[p_start-start:p_end-start] = row['signalValue']
            
            # for "-" strand
            elif strand == "-" and (p_strand == "." or p_strand == "-"):
                profile[end-p_start:end-p_end] = row['signalValue'] ### !!! ### to be tested
            
        output_df = pd.concat([output_df, profile], axis=1)
        if fillna==True:
            output_df = output_df.fillna(0.0)
    output_df = output_df.T
    output_df.index.name = "region"
    return output_df.reset_index()

### level 0
def getMultipleMatrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, 
                        fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    
    '''Get score matrices for positions in given regions (with optional flanks) from multiple BigWig files.
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
    '''

    bed_plus, bed_minus = bed_split_strands(bed_df)
    plus_dict = {}
    minus_dict = {}
    for bw_plus, bw_minus in zip(bw_paths_plus, bw_paths_minus):
        try:
            plus_dict[bw_plus] = matrixFromBigWig(bw_path=bw_plus, bed_df=bed_plus, flank_5=flank_5, flank_3=flank_3, fill_na=fill_na, pseudocounts=pseudocounts, align_3end=align_3end)
            minus_dict[bw_minus] = matrixFromBigWig(bw_path=bw_minus, bed_df=bed_minus, flank_5=flank_5, flank_3=flank_3, fill_na=fill_na, pseudocounts=pseudocounts, align_3end=align_3end)
            if normalize_libsize:
                with pyBigWig.open(bw_plus) as bwp:
                    libsize_plus = int(bwp.header()['sumData'])
                with pyBigWig.open(bw_minus) as bwm:
                    libsize_minus = int(bwm.header()['sumData'])
                libsize = libsize_plus+libsize_minus
                plus_dict[bw_plus] = plus_dict[bw_plus].set_index('region').div(libsize).reset_index()
                minus_dict[bw_minus] = minus_dict[bw_minus].set_index('region').div(libsize).reset_index()
        except:
            e = str(bw_plus)+ " " + str(bw_minus)
            print(f"Error: {e}")
            # continue
    return join_strand_matrices(plus_dict, minus_dict)

def getMultipleMatricesFromPeak(peak_paths=[], bed_df=pd.DataFrame,
                 g='references/GRCh38.primary_assembly.genome.cleaned.len', 
                 flank_5=0, flank_3=0):
    '''Get score matrices for positions in given regions (with optional flanks) from multiple peak files.
    Matrix rows correspond to regions in BED. Columns correpond to nucleotide positions in regions + flanks.
    :param peak_paths: list of paths to peak files (broadPeak or narrowPeak format)
    :type peak_paths: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param g: Path to genome file (chromosome lengths), defaults to 'references/GRCh38.primary_assembly.genome.cleaned.len'
    :type g: str
    :param flank_5: length of 5'flank to extend BED regions by, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to 0
    :type flank_3: int, optional
    :return: dictionary of score matrices for each peak file
    :rtype: dict
    '''

    out_dict = {}
    for path in peak_paths:
        name = path.split("/")[-1].replace('_peaks.narrowPeak','').replace('_peaks.broadPeak','')
        out_dict[name] = peak2matrice(bed_df=bed_df, peak_file_path=path, 
                 g=g, flank_5=flank_5, flank_3=flank_3)
    
    return out_dict

def metaprofile(matrix_dict, agg_type='mean', normalize_internal=False, subset=None):
    '''Calculate metaprofiles from score matrices by aggregating each position in all regions.

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

    if isinstance(subset, pd.DataFrame):
        matrix_dict = {key: value[value['region'].isin(subset.index)] for key, value in matrix_dict.items()}
    elif isinstance(subset, list):
        matrix_dict = {key: value[value['region'].isin(subset)] for key, value in matrix_dict.items()}


    if normalize_internal:
        dropped = {key: value.drop('region', axis=1).div(value.sum(axis=1,numeric_only=True),axis=0) for key, value in matrix_dict.items()}
        return pd.DataFrame({key: value.agg(agg_type,numeric_only=True) for key, value in dropped.items()})

    else:
        return pd.DataFrame({key: value.agg(agg_type, numeric_only=True) for key, value in matrix_dict.items()})

def regionScore(bw_paths_plus, bw_paths_minus, bed_df, agg_type='sum', flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True):
    '''Calculate coverage or statistic for multiple regions in multiple BigWig files.

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
    '''

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

### level 1
def binMultipleMatrices(mm={}, bins=[50, 10, 50], bed_df=pd.DataFrame(), flank_5=None, flank_3=None):
    '''Bin multiple matrices of tRNA profiles into a single dataframe
      
    :param mm: dictionary of matrices of gene profiles 
    :type mm: dict
    :param bins: list of three integers, defaults to [50, 10, 50]
    :type bins: list, optional
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions, defaults to pd.DataFrame()
    :type bed_df: pandas.DataFrame, optional
    :param flank_5: length of 5'flank to extend BED regions by, defaults to None
    :type flank_5: int, optional
    :type bed_df: pandas.DataFrame, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to None
    :type flank_3: int, optional
    :type bed_df: pandas.DataFrame, optional
    :raises ValueError: if bins is not a list of three integers
    :raises ValueError: if input_value is not an integer or a DataFrame
    :return: dictionary of binned matrices of tRNA profiles
    :rtype: dict
    '''

    if len(bins) != 3:
        raise ValueError(f"{bins} takes three numbers")
        
    def process_input(input_value):
        if isinstance(input_value, pd.DataFrame):
            return tt.methods.bed2len(input_value)
        elif isinstance(input_value, int):
            return pd.Series(index=bed_df.iloc[:, 3].rename('region'), data=input_value) #uses region as index
        else:
            raise ValueError(f"{input_value} should be either an integer or a DataFrame")

    # Concatenate all the results into length_df
    length_df = pd.concat([process_input(flank_5),
                           tt.methods.bed2len(bed_df),
                           process_input(flank_3)], axis=1)

    results_mm = {} 
    # Iterate over the datasets (keys) of the mm dictionary
    for bw_key, value in mm.items():
        results_df = pd.DataFrame(columns=bed_df.iloc[:, 3].rename('region'))

        # bin individual profiles of tRNAs
        for region,row in value.set_index('region').iterrows():
            results_df[region] = tt.profiles.binCollect3(s1=row,
                                         lengths=length_df.loc[region].tolist(),
                                         bins=bins)
        
        results_df.index = np.arange(-bins[0],bins[1]+bins[2])
        results_mm[bw_key] = results_df.T.reset_index()
    
    return results_mm

## helper functions
def selectSubsetMM(matrix_dict, subset=None):
    '''Select a subset of regions from the matrix dictionary.

    :param matrix_dict: Dictionary containing score matrices.
    :type matrix_dict: dict
    :param subset: DataFrame or list containing the subset of regions to select, defaults to None
    :type subset: Union[pd.DataFrame, list], optional
    :return: Dictionary containing the selected subset of score matrices.
    :rtype: dict
    '''
    if isinstance(subset, pd.DataFrame):
        return {key: value[value['region'].isin(subset.index)] for key, value in matrix_dict.items()}
    elif isinstance(subset, list):
        return {key: value[value['region'].isin(subset)] for key, value in matrix_dict.items()}

# def regionScore2(bw_paths_plus, bw_paths_minus, bed_df, agg_type='sum', flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True):
#     """
#     Calculate coverage or statistic for multiple regions in multiple BigWig files.

#     :param bw_paths_plus: list of paths to BigWig files (+ strand)
#     :type bw_paths_plus: list
#     :param bw_paths_minus: list of paths to BigWig files (- strand)
#     :type bw_paths_minus: list
#     :param bed_df: dataframe in BED format containing genomic coordinates of target regions
#     :type bed_df: pandas.DataFrame
#     :param agg_type: operation to perform on region scores. Available options: 'sum', 'mean', 'median', defaults to 'sum'
#     :type agg_type: str, optional
#     :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
#     :type flank_5: int, optional
#     :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
#     :type flank_3: int, optional
#     :param fill_na: If true, NaNs will be replaced with zeroes (recommended, as pyBigWig reports positions with 0 score as NaN), defaults to True
#     :type fill_na: bool, optional
#     :param pseudocounts: pseudocounts to add to datapoints, defaults to None
#     :type pseudocounts: float, optional
#     :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
#     :type normalize_libsize: bool, optional
#     :return:  DataFrame with calculated scores. Rows are regions/genes, columns are BigWig files.
#     :rtype: pandas.DataFrame
#     """    
#     if agg_type not in ['mean', 'median', 'sum']:
#         raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")
    
#     bed_plus, bed_minus = bed_split_strands(bed_df)
    
#     out_df_list = []
#     # Get scores from BigWig
#     for bigwig in zip(bw_paths_plus, bw_paths_minus):
#         out_plus = matrixFromBigWig(
#             bigwig[0],
#             bed_plus,
#             flank_5=flank_5,
#             flank_3=flank_3,
#             fill_na=fill_na,
#             pseudocounts=pseudocounts,
#             normalize_libsize=normalize_libsize
#             ).set_index('region')
#         out_minus = matrixFromBigWig(
#             bigwig[0],
#             bed_minus,
#             flank_5=flank_5,
#             flank_3=flank_3,
#             fill_na=fill_na,
#             pseudocounts=pseudocounts,
#             normalize_libsize=normalize_libsize
#             ).set_index('region')
#         # Merge strands and aggregate
#         out_merged = pd.concat([out_plus, out_minus]).agg(func=agg_type, axis=1)
#         out_merged.columns = [bigwig[0]]
#         print(out_merged)
#         out_df_list.append(out_merged)
#     # Join df's as columns
#     return pd.concat(out_df_list, axis=1)