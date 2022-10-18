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

def list_paths_in_current_dir(suffix=str(), stdin=False):
    '''
    :param suffix: str() lists paths in current directory ending with an indicated suffix only
    :param stdin: boolean() if True read from standard input instead current directory
    :return: list() of paths
    '''
    #choosing between curr dir and std input
    if stdin == False:
        where = os.listdir('.')
    elif stdin == True:
        where = sys.stdin

    paths = [f for f in where if os.path.isfile(f) and f.endswith(suffix)]
    return paths

def bashCommand(bashCommand=str()):
    '''Run command in bash using subprocess.call()'''
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
    :return: fraction of GC content, float
    '''
    return float(len(dataset[dataset['nucleotide'].isin(calFor)]))/float(len(dataset))

alt_map = {'ins':'0'}
complement_DNA = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
complement_RNA = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}

def reverse_complement_DNA(seq):
    '''Reverse complement

    :param seq: str
    :return: str
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
    '''Reverse complement

    :param seq: str
    :return: str
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
    '''Reverse complement

    :param seq: str
    :return: str
    '''
    if "U" in seq:
        return reverse_complement_RNA(seq)
    else:
        return reverse_complement_DNA(seq)

def randomDNAsingle(length=int(), letters="CGTA"):
    '''Random generator of nucleotide sequence

    :param length: int()
    :param letters: str() with letters that will be used
    :return: str()
    '''

    return''.join(random.choices(letters, k=length))

def randomDNAall(length=int(), letters="CGTA"):
    '''Generates all possible random sequences of a given length

    :param length: int()
    :param letters: str() with letters that will be used
    :return: list() of str()
    '''

    output_list = []
    for i in itertools.product(list(letters), repeat=length):
        output_list.append("".join(i))
    return output_list

def rollingGC(s=pd.Series, window=10): #rolling window, smoothing of data
    '''Calculates GC from sequence, uses 'boxcar' window

    :param s: Series containing sequence
    :param window: window size for GC calculation
    :return: Series with GC calculated, center=False
    '''
    #be aware of win_type and center parameters
    return s.replace(['G','C'],1).replace(['T','A'],0).rolling(window=window, win_type='boxcar',center=False).mean()

def letterContent(s="", letter="A"):
    return round(len([i for i in s if i==letter])/len(s),2)

    
################################################
#############        importing data

def read_list(filepath=''):
    '''Read list from file. Each row becomes item in the list.
    
    :param filepath: str
    :return: list
    '''

    txt_file = open(filepath, "r")
    content_list = txt_file.read().splitlines()
    txt_file.close()
    return content_list

def read_tabFile(nameElem="", path="", toLoad="", toClear=[], toAdd="", df=None, overwrite=False):
    '''
    Read tab files with common first column
    :param nameElem: str, present in all files
    :param path: str, path to directory with files
    :param toLoad: str, to be present in file name
    :param toClear: str, will be removed from file name
    :param toAdd: str, to be added to file name
    :param df: DataFrame, to be appended; default=None
    :param overwrite: boolean, allows for overwriting during appending, default = False
    :return: DataFrame
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
    '''Reads multiple HTSeq tab files to one DataFrame

    :param path: str, path to directory with files
    :param toClear: str, will be removed from file name
    :param toAdd: str, to be added to file name
    :param df: DataFrame, to be appended; default=None
    :param overwrite: boolean, allows for overwriting during appending, default = False
    :return: DataFrame
    '''
    return read_tabFile(nameElem='_STARLog.final.out', path=path,
                        toClear=toClear, toAdd=toAdd, df=df, overwrite=overwrite)


def read_HTSeq_output(path="", toLoad="classes", toClear=[], toAdd="", df=None, overwrite=False):
    '''Reads multiple HTSeq tab files to one DataFrame

    :param path: str, path to directory with files
    :param toClear: str, will be removed from file name
    :param toAdd: str, to be added to file name
    :param df: DataFrame, to be appended; default=None
    :param overwrite: boolean, allows for overwriting during appending, default = False
    :return: DataFrame
    '''
    return read_tabFile(nameElem='_hittable.tab', path=path, toLoad=toLoad,
                        toClear=toClear, toAdd=toAdd, df=df, overwrite=overwrite)

def readSalmon(nameElem="", path="", toLoad="", toClear=[], toAdd="", column='NumReads', df=None, overwrite=False):
    '''

    :param nameElem: str, elem to load
    :param path: str
    :param toLoad: str, additional param for filtering, by default equal to nameElem
    :param toClear: str
    :param toAdd: str
    :param df: pd.DataFrame
    :param overwrite: boolean, default=False
    :return:
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

def read_featureCount(nameElem="", path="", toLoad="", toClear=[], toAdd="", df=None, overwrite=False):
    '''
    Read tab files with common first column
    :param nameElem: str, present in all files
    :param path: str, path to directory with files
    :param toLoad: str, to be present in file name (optional)
    :param toClear: str, will be removed from file name
    :param toAdd: str, to be added to file name
    :param df: DataFrame, to be appended; default=None
    :param overwrite: boolean, allows for overwriting during appending, default = False
    :return: DataFrame
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

################################################
#############        handling multiple experiments

def define_experiments(paths_in, whole_name=False, strip='_hittable_reads.txt'):
    '''Parse file names and extract experiment name from them

    :param paths_in: str()
    :param whole_name: boolean() default False. As defaults script takes first 'a_b_c'
    :param strip: str() to strip from filename.
    :return: list() of experiment names, list() of paths.
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

    :param name:
    :param additional_tags: list of tags
    :param output: default 'root' ; print other elements when 'all'
    :param order: defoult 'b_d_e_p' b-bait; d-details, e-experiment, p-prefix
    :return: list of reordered name
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

def cleanNames(df=pd.DataFrame(), additional_tags=[]):
    '''Cleans some problems with names if exist

    :param df: DataFrame() where names of columns are name of experiments
    :param additional_tags: list()
    :return: DataFrame() with new names
    '''

    for tag in additional_tags: df.columns = [f.replace(tag, tag+'HTP') for f in list(df.columns.values)]
    df.columns = [f.replace('HTPHTP', 'HTP').replace('HTPHTP', 'HTP') for f in list(df.columns.values)]
    return df

def indexOrder(df=pd.DataFrame(), additional_tags=list(), output='root', order='b_d_e_p'):
    '''Apply expNameParser to whole DataFrame

    :param df: DataFrame() where names of columns are name of experiments
    :param additional_tags: list()
    :param output:
    :param order: str() default 'b_d_e_p' b-bait; d-details, e-experiment, p-prefix
    :return: DataFrame() with new names
    '''

    df = cleanNames(df, additional_tags=additional_tags)
    df.columns = [expNameParser(f, additional_tags=additional_tags, order=order) for f in list(df.columns.values)]
    return df.sort_index(axis=1)

def filterExp(datasets, let_in=[''], let_out=['wont_find_this_string'],verbose=False):
    '''Returns object with filtered columns/keys.

    :param datasets: DataFrame() or dict() with exp name as a key
    :param let_in: list() with elements of name to filter in
    :param let_out: list() with elements of name to filter out
    :return: DataFrame() or dict()
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

    :param s1: Series,
    :return: DataFrame
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
    '''Parse CRAC names and annotates them using on of following features
    ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3', 'sample','sampleRep']

    :param df: DataFrame
    :param use: str, choose from ['expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3', 'sample','sampleRep'], default = 'protein'
    :param toDrop: list of word in CRAC name that will qualify the sample to rejection, default = []
    :return: DataFrame with added column ['group']
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

    :param input_df: DataFrame
    :param smooth: boolean, if True apply smoothing window, default=True
    :param window: int, smoothing window, default 10
    :param win_type: str type of smoothing window, default "blackman"
    :return: DataFrame
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

def normalize(df=pd.DataFrame, log2=False, pseudocounts=0.1):
    '''

    :param df: DataFrame
    :param log2: boolean, default=False
    :param pseudocounts: float, default=0.1
    :return:
    '''
    df = df.add(pseudocounts)
    df = df / df.sum()
    df = df.multiply(1000000)
    if log2==True:
        return df.apply(np.log2)
    else:
        return df

def quantileCategory(s1=pd.Series(dtype=float), q=4):
    '''Quantile-based discretization function based on pandas.qcut function.

    :param s1: Series()
    :param q: int() number of quantiles: 10 for deciles, 5 for quantiles, 4 for quartiles, etc., default q=4
    :return: Series
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

    :param data: DataFrame
    :param n_components: int, default 2
    :return: tuple consisting of DataFrame with PCA results and a list of PC values
    :rtype: tuple
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
    '''
    Assigns n clusters to the data using KMeans algorithm
    :param df: DataFrame
    :param n: no. of clusters, int
    :return:
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
    '''
    :return: timestamp as a str()
    '''
    return str(time.strftime("%Y%m%d_%H%M%S"))

def timestampRandomInt():
    '''
    :return: timestamp and random number as a str()
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
