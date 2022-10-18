from scipy.signal import argrelextrema
import numpy as np
from collections import OrderedDict
import pandas as pd
from warnings import warn
import pyBigWig

### PROFILES (ProfileAnalyser) ###
# functions to deal with profiles of individual genes

#importing data
##concat file
def parseConcatFile(path, gtf, use='reads', RPM=False, ranges=1000):
    '''Parse concat file

    :param path: str with path of the concat file
    :param gtf: pyCRAC.GTF2 object with GTF and TAB files loaded
    :param use: str with name of column tu use ['reads', 'substitutions', 'deletions'], default "reads"
    :param RPM: boolean, default False
    :param ranges: int flanks to be added for the gene, default 0
    :return: dict of DataFrames; using gene name as a key
    '''
    # load concat
    columns = ['gene', 'position', 'nucleotide', 'reads', 'substitutions', 'deletions', 'experiment', 'reads_pM',
               'substitutions_pM', 'deletions_pM']
    df_rawConcat = pd.read_csv(path, sep="\t", names=columns, comment="#")

    # use collumn
    if RPM == True: use = use + "_pM"

    # output dict
    output = {}

    # genes
    gene_list = df_rawConcat['gene'].unique()

    # loop trough genes
    for gene_name in gene_list:
        geneLen = gtf.geneLength(gene_name)
        # preparing dataframe
        df_output = pd.DataFrame()
        df_output['position'] = pd.Series(list(range(-ranges, 0)) + list(range(1, geneLen + ranges + 1)))
        df_output['nucleotide'] = pd.Series(list(gtf.genomicSequence(gene_name, ranges=ranges)))

        # loop trough experiments
        for e, df in df_rawConcat.groupby('experiment'):
            df_output[e] = df.reset_index()[use]

            # check 'nucleotide' collumn
            if df_output['nucleotide'].tolist() != df['nucleotide'].tolist():
                warn("For " + gene_name + " gene list of nucleotides does not match/", UserWarning)

        output[gene_name] = df_output
    return output

##BigWig
def stripBigWignames(files=[]):
    '''Strip "_rev.bw" and "_fwd.bw" form file names

    :param files: list of filenames
    :return: list of unique names
    '''
    # returns uniq names for *fwd.bw and *rev.bw files
    return list(set([f.replace("_rev.bw", "").replace("_fwd.bw", "") for f in files]))


def dictBigWig(files=[], path="", strands=True):
    '''Preloads BigWig files to memory using pyBigWig tools

    :param files: list of files
    :param path: str
    :param strands: boolean, default True
    :return: dict or dict, dict
    '''
    if strands == True:
        # output dict
        output_dict_fwd = OrderedDict()
        output_dict_rev = OrderedDict()

        strippedNames = stripBigWignames(files=files)  # list of common roots for strand specific data
        for f in strippedNames:
            output_dict_fwd[f] = pyBigWig.open(path + f + "_fwd.bw")
            output_dict_rev[f] = pyBigWig.open(path + f + "_rev.bw")
        return output_dict_fwd, output_dict_rev

    if strands == False:
        # output dict
        output_dict = {}
        for f in files:
            output_dict[f.replace('.bw', '')] = pyBigWig.open(path + f)
        return output_dict

def geneFromBigWig(gene_name, gtf, bwFWD={}, bwREV={}, toStrip="", ranges=0):
    '''Pulls genome coverage from BigWig data for a given gene. One BigWig file -> one column.

    :param gene_name: str
    :param gtf: pyCRAC.GTF2 object with GTF and TAB files loaded
    :param bwFWD: dict of pyBigWig objects
    :param bwREV: dict of pyBigWig objects
    :param toStrip: str of name to be stripped
    :param ranges: int flanks to be added for the gene, default 0
    :return: DataFrame
    '''
    # preparing output
    df_t1 = pd.DataFrame()
    df_t1["nucleotide"] = "_".join(gtf.genomicSequence(gene_name, ranges=ranges)).split("_")

    # coordinates from GTF file
    strand, chromosome, coordinates = gtf.strand(gene_name), gtf.chromosome(gene_name), gtf.chromosomeCoordinates(
        gene_name)

    # extracting data
    if strand == "+":
        for name in bwFWD:
            bw = bwFWD[name]
            s1 = pd.Series(bw.values(chromosome, min(coordinates) - ranges, max(coordinates) + ranges))
            df_t1[name.replace(toStrip, '')] = s1

    if strand == "-":
        for name in bwREV:
            bw = bwREV[name]
            s1 = pd.Series(bw.values(chromosome, min(coordinates) - ranges, max(coordinates) + ranges)[::-1])
            df_t1[name.replace(toStrip, '')] = s1

    return df_t1

def FoldingFromBigWig(gene_name, gtf, bwFWD={}, bwREV={}, ranges=0,offset=15,fold="dG65nt@30C"):
    '''Pulls folding information from BigWig folding data for a given gene.

    :param gene_name: str
    :param gtf: pyCRAC.GTF2 object with GTF and TAB files loaded
    :param bwFWD: dict of pyBigWig objects
    :param bwREV: dict of pyBigWig objects
    :param ranges: int flanks to be added for the gene, default 0
    :param offset: int to offset folding data, default 15
    :param fold: name of output column, default="dG65nt@30C"
    :return: DataFrame
    '''
    # preparing output
    df_t1 = pd.DataFrame()

    # coordinates from GTF file
    strand, chromosome, coordinates = gtf.strand(gene_name), gtf.chromosome(gene_name), gtf.chromosomeCoordinates(
        gene_name)

    # extracting data
    if strand == "+":
        for name in bwFWD:
            bw = bwFWD[name]
            s1 = pd.Series(bw.values(chromosome, min(coordinates) - ranges, max(coordinates) + ranges))
            # shifting to the extrustion point
            df_t1[fold+"_add"+str(offset)+"nt"] = pd.Series((offset * [np.nan]) + s1.tolist())

    if strand == "-":
        for name in bwREV:
            bw = bwREV[name]
            s1 = pd.Series(bw.values(chromosome, min(coordinates) - ranges, max(coordinates) + ranges)[::-1])
            # shifting to the extrustion point
            df_t1[fold+"_add"+str(offset)+"nt"] = pd.Series((offset * [np.nan]) + s1.tolist())

    return df_t1

##normalization
def pseudocounts(df=pd.DataFrame, value=0.01, drop=True):
    '''Add pseudocounts to data

    :param df: DataFrame
    :param value: float, default 0.01
    :param drop: boolean, if True drop 'position' and 'nucleotide' columns, default True
    :return: DataFrame
    '''
    # drop additional columns
    cols = df.columns.tolist()
    if 'position' in cols:
        s1_pos = df['position']
        df = df.drop('position', axis=1)
    if 'nucleotide' in cols:
        s2_nt = df['nucleotide']
        df = df.drop('nucleotide', axis=1)

    # pseudocounts
    df = df.fillna(0).replace(0, value)

    # output
    if drop == True:
        return df
    elif drop == False:
        df['position'] = s1_pos
        df['nucleotide'] = s2_nt
        return df


def ntotal(df=pd.DataFrame, drop=True):
    '''Normalize data in DataFrame to fraction of total column

    :param df: DataFrame
    :param drop: boolean, if True drop 'position' and 'nucleotide' columns, default True
    :return: DataFrame
    '''
    # drop additional columns
    cols = df.columns.tolist()
    if 'position' in cols:
        s1_pos = df['position']
        df = df.drop('position', axis=1)
    if 'nucleotide' in cols:
        s2_nt = df['nucleotide']
        df = df.drop('nucleotide', axis=1)

    # normalization
    df = df / df.sum()

    # output
    if drop == True:
        return df
    elif drop == False:
        df['position'] = s1_pos
        df['nucleotide'] = s2_nt
        return df

################################################
#############        major trxtools

def preprocess(input_df=pd.DataFrame(), let_in=[''], let_out=['wont_find_this_string'],
              stats=False, smooth=True , window=10, win_type='blackman'):
    '''Combines methods.filterExp and expStats. Returns DataFrame with choosen experiments, optionally apply smoothing and stats

    :param input_df: DataFrame
    :param let_in: list of words that characterize experiment, default ['']
    :param let_out: list of words that disqualify experiments, default ['wont_find_this_string']
    :param stats: boolean, if True return stats for all experiments, default False
    :param smooth: boolean, if True apply smoothing window, default True
    :param window: int smoothing window, default 10
    :param win_type: str type of smoothing window, default "blackman"
    :return: DataFrame with 'mean', 'median', 'min', 'max' and quartiles if more than 2 experiments
    '''
    # filtering
    selected = [d for d in list(input_df.columns.values) if
                all(i in d for i in let_in) and all(o not in d for o in let_out)]
    print("Experiments: ")
    print(selected)

    # selecting data
    working_df, result_df = pd.DataFrame(), pd.DataFrame()
    for f in selected:
        if smooth == True:
            working_df[f] = input_df[f].rolling(window, win_type=win_type, center=True).mean()
        else:
            working_df[f] = input_df[f]

    # preparing output
    if stats == False:
        return working_df
    else:
        for function in ['mean', 'median', 'min', 'max']: result_df[function] = getattr(working_df, function)(
            axis=1)  # calculates using pandas function listed in []
        if len(selected) > 2:  # calculating quartiles only in more than two experiments
            result_df['q1'], result_df['q3'] = working_df.quantile(q=0.25, axis=1), working_df.quantile(q=0.75, axis=1)
        return result_df

def calculateFDR(data=pd.Series(dtype=float), iterations=100, target_FDR=0.05):
    '''Calculates False Discovery Rate (FDR) for a given dataset.

    :param data: Series
    :param iterations: int, default 100
    :param target_FDR: float, detault 0.05
    :return: Series
    '''

    normalized_data = data / data.sum()  # normalize data
    random_df = pd.DataFrame(np.random.rand(len(data), iterations))  # generating random datasets
    normalized_random_df = random_df / random_df.sum()  # normalize random df

    # getting Discovery Rate (DR)
    DR_df = normalized_random_df.T >= normalized_data  # random discovery rate
    DR = DR_df.sum() / iterations

    # comparing DR to False Discovery Rate (FDR)
    FDR = DR <= target_FDR
    return data * FDR

def findPeaks(s1=pd.Series(dtype=float), window=1, win_type='blackman', order=20):
    '''Find local extrema using SciPy argrelextrema function

    :param s1: Series data to localize peaks
    :param window: int, To smooth data before peak-calling. default 1 (no smoothing)
    :param win_type: str type of smoothing window, default "blackman"
    :param order: int minimal spacing between peaks, argrelextrema order parameter, default  20
    :return: list of peaks
    '''

    # smoothing
    if window > 1:
        s1 = s1.rolling(window, win_type=win_type, center=True).mean()

    # find peaks
    output = argrelextrema(data=s1.to_numpy(), comparator=np.greater, order=order)[0]
    return list(output)

def findTroughs(s1=pd.Series(dtype=float), window=1, win_type='blackman', order=20):
    '''Find local minima using SciPy argrelextrema function

    :param s1: Series data to localize peaks
    :param window: int, To smooth data before trough-calling. default 1 (no smoothing)
    :param win_type: str type of smoothing window, default "blackman"
    :param order: int minimal spacing between min, argrelextrema order parameter, default 20
    :return: list of troughs
    '''

    # smoothing
    if window > 1:
        s1 = s1.rolling(window, win_type=win_type, center=True).mean()

    # find troughs
    output = argrelextrema(data=s1.to_numpy(), comparator=np.less, order=order)[0]
    return list(output)


def compare1toRef(ref, dataset=pd.Series(dtype=float), ranges='mm', heatmap=False, relative=False):
    '''Takes Series and compare this with reference DataFrame()

    :param ref: str with path to csv file or DataFrame
    :param dataset: Series
    :param ranges: mm : min-max or qq : q1-q3
    :param heatmap: boolean, heatmap=False: Dataframe with(reference_above_experiment minimum etc.): rae_min, rae_max, ear_min, ear_max;
            heatmap=True: Series of differences to plot heatmap
    :param relative: boolean, only for heatmap, recalculates differences according to the peak size. Warning: negative values are in range -1 to 0 but positive are from 0 to values higher than 1
    :return: Dataframe (heatmap=False) or Series (heatmap=True)
    '''

    ranges_dict = {'mm': ['min', 'max'], 'qq': ['q1', 'q3']}

    # preparing dataframe and reference
    differences_df, return_df = pd.DataFrame(), pd.DataFrame()

    #handling reference plot
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    differences_df['exp'] = dataset
    differences_df['ref_min'] = reference[ranges_dict[ranges][0]]  # choosing q1 or min
    differences_df['ref_max'] = reference[ranges_dict[ranges][1]]  # choosing q3 or max

    ## finding differences (only positive value indicate difference)
    differences_df['ref_above_exp'] = differences_df['ref_min'] - differences_df['exp']
    differences_df['exp_above_ref'] = differences_df['exp'] - differences_df['ref_max']

    ## filling differences
    for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max', 'diff_ear', 'diff_rae', 'diff']: differences_df[f] = 0
    differences_df['rae_min'][differences_df['ref_above_exp'] > 0] = differences_df['exp']
    differences_df['rae_max'][differences_df['ref_above_exp'] > 0] = differences_df['ref_min']
    differences_df['ear_min'][differences_df['exp_above_ref'] > 0] = differences_df['ref_max']
    differences_df['ear_max'][differences_df['exp_above_ref'] > 0] = differences_df['exp']

    # combining differences for heatmap
    differences_df['diff_ear'] = differences_df['ear_max'] - differences_df['ear_min']
    differences_df['diff_rae'] = differences_df['rae_min'] - differences_df['rae_max']
    differences_df['diff'][differences_df['diff_ear'] > 0] = differences_df['diff_ear']
    differences_df['diff'][differences_df['diff_rae'] < 0] = differences_df['diff_rae']

    if heatmap == True and relative == True:
        differences_df['ref_median'] = reference['median']
        differences_df['rel_diff'] = differences_df['diff'] / differences_df['ref_median']
        return differences_df['rel_diff']  # return Series
    elif heatmap == True:
        return differences_df['diff']  # return Series
    if heatmap == False:
        return_df['exp'] = dataset
        for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: return_df[f] = differences_df[f]
        return return_df  # return Dataframe


def compareMoretoRef(ref, dataset=pd.DataFrame(), ranges='mm'):
    '''Takes Dataframe created by filter_df and compare this with reference DataFrame

    :param ref: str with path to csv file or DataFrame
    :param dataset: Series
    :param ranges: mm : min-max or qq : q1-q3
    :return: Dataframe
    '''
    ranges_dict = {'mm': ['min', 'max'], 'qq': ['q1', 'q3']}
    # preparing dataframe and reference
    differences_df, return_df = pd.DataFrame(), pd.DataFrame()

    #handling reference plot
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    if len(dataset.columns) == 4:  # if only two experiments
        differences_df['exp_min'] = dataset['min']
        differences_df['exp_max'] = dataset['max']
    else:  # if more than two exp
        differences_df['exp_min'] = dataset[ranges_dict[ranges][0]]  # choosing q1 or min
        differences_df['exp_max'] = dataset[ranges_dict[ranges][1]]  # choosing q3 or max
    differences_df['ref_min'] = reference[ranges_dict[ranges][0]]  # choosing q1 or min
    differences_df['ref_max'] = reference[ranges_dict[ranges][1]]  # choosing q3 or max

    ## finding differences (only positive value indicate difference)
    differences_df['ref_above_exp'] = differences_df['ref_min'] - differences_df['exp_max']
    differences_df['exp_above_ref'] = differences_df['exp_min'] - differences_df['ref_max']

    ## filling differences
    for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: differences_df[f] = 0
    differences_df['rae_min'][differences_df['ref_above_exp'] > 0] = differences_df['exp_max']
    differences_df['rae_max'][differences_df['ref_above_exp'] > 0] = differences_df['ref_min']
    differences_df['ear_min'][differences_df['exp_above_ref'] > 0] = differences_df['ref_max']
    differences_df['ear_max'][differences_df['exp_above_ref'] > 0] = differences_df['exp_min']
    #     return_df = dataset
    for f in ['rae_min', 'rae_max', 'ear_min', 'ear_max']: return_df[f] = differences_df[f]
    return return_df  # returns Dataframe


################################################
#############        plotting

# def plot_ChIP(df_sense=pd.DataFrame(), df_anti=pd.DataFrame(), title=None, start=None, stop=None, figsize=(15, 6),
#               ylim=(-0.001, 0.001), s_color='red', as_color='blue', h_lines=list(), lc='black', dpi=150,
#               csv_path='/home/tturowski/notebooks/RDN37_reference_collapsed.csv', color='green'):
#     '''Function creates plot similar to box plot: median, 2 and 3 quartiles and min-max range
#     csv_path: str()
#         Path to CRAC or other reference file
#     title: str()
#     start: int()
#     stop: int()
#     figsize: tuple()
#     ylim: tuple()
#         OY axes lim - def (None,0.01)
#     color: str()
#         plot color
#     h_lines: list()
#         Optional. list() of horizontal lines
#     lc: str()
#         color of horizontal lines
#     Returns
#     -------
#     None
#     '''
#     reference = pd.read_csv(csv_path, index_col=0).drop('nucleotide', 1)
#     s2 = reference[start:stop]
#     # plotting reference dataset
#     fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
#     plt.title(title)
#     ax1.set_xlabel('position')
#     ax1.set_ylabel('fraction of reads')
#     ax1.set_ylim(ylim)
#     # reference plot
#     ax1.plot(s2.index, s2['median'], color=color)
#     if set(['q1', 'q3']).issubset(list(s2.columns.values)):
#         ax1.fill_between(s2.index, s2['q1'], s2['q3'], label='range (2nd-3rd quartile)', color=color, alpha=0.2)
#     if set(['min', 'max']).issubset(list(s2.columns.values)):
#         ax1.fill_between(s2.index, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)
#
#     for i in [i for i in h_lines if i in range(start, stop)]: ax1.axvline(i, color=lc)
#     ax1.axhline(0, color='black', alpha=0.7, ls='dashed', lw=1)
#
#     # ChIP sense
#     c1 = df_sense[start:stop]
#     ax1.plot(s2.index, c1['median'], color=s_color)
#     if set(['q1', 'q3']).issubset(list(c1.columns.values)):
#         ax1.fill_between(s2.index, c1['q1'], c1['q3'], label='range (2nd-3rd quartile)', color=s_color, alpha=0.2)
#     if set(['min', 'max']).issubset(list(c1.columns.values)):
#         ax1.fill_between(s2.index, c1['min'], c1['max'], label='range (min-max)', color=s_color, alpha=0.07)
#     # ChIP antisense
#     c2 = df_anti[start:stop] * -1
#     ax1.plot(s2.index, c2['median'], color=as_color)
#     if set(['q1', 'q3']).issubset(list(c2.columns.values)):
#         ax1.fill_between(s2.index, c2['q1'], c2['q3'], label='range (2nd-3rd quartile)', color=as_color, alpha=0.2)
#     if set(['min', 'max']).issubset(list(c2.columns.values)):
#         ax1.fill_between(s2.index, c2['min'], c2['max'], label='range (min-max)', color=as_color, alpha=0.07)
#
#     ax1.legend()
#     return None


################################################
#############        other

def save_csv(data_ref=pd.DataFrame(), datasets=pd.DataFrame(), path=None):
    '''Saves Dataframe to csv

    :param data_ref: DataFrame with ``['position']`` and ``['nucleotide']`` columns
    :param datasets: DataFrame containinig experimental data only
    :param path: str, Optional: path to save csv. Default None
    :return: DataFrame
    '''

    reference = pd.DataFrame()
    # if 'position'

    reference['position'], reference['nucleotide'] = data_ref['position'], data_ref['nucleotide']
    reference['mean'], reference['std'] = datasets.mean(axis=1), datasets.std(axis=1)
    reference['median'] = datasets.median(axis=1)
    reference['q1'], reference['q3'] = datasets.quantile(q=0.25,axis=1), datasets.quantile(q=0.75, axis=1)
    reference['min'], reference['max'] = datasets.min(axis=1), datasets.max(axis=1)

    if path:
        reference.to_csv(path)  ## reference plot
    return reference