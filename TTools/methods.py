import subprocess, time, random
import os, sys, re, itertools
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

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
    subprocess.call(bashCommand, shell=True)

################################################
#############        DNA/RNA

def calGC(dataset=pd.DataFrame(), calFor=['G','C']):
    '''
    Returns GC content in a given string - uses ['nucleotide'] column

    :param dataset: DataFrame() with "nucleotide" column
    :return: fraction of GC content
    '''
    return float(len(dataset[dataset['nucleotide'].isin(calFor)]))/float(len(dataset))

alt_map = {'ins':'0'}
complement_DNA = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
complement_RNA = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}

def reverse_complement_DNA(seq):
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
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement_RNA.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def reverse_complement(seq):
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


################################################
#############        handling multiple experiments

def define_experiments(paths_in, whole_name=False, strip='_hittable_reads.txt'):
    '''
    Parse file names and extract experiment name from them

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
    '''
    Function handles experiment name; recognizes AB123456 as experiment date; BY4741 or HTP or given string as bait protein
    :param name:
    :param additional_tags: list of tags
    :param output: default 'root' ; print other elements when 'all'
    :param order: defoult 'b_d_e_p' b-bait; d-details, e-experiment, p-prefix
    :return: reordered name
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

def cleanNames(df=pd.DataFrame(), additional_tags=list()):
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

def filterExp(datasets, let_in=[''], let_out=['wont_find_this_string']):
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
        return output_df

    #for dict()
    elif isinstance(datasets, dict):
        output_dict = dict()
        for f in [d for d in list(datasets.keys()) if all(i in d for i in let_in) and all(o not in d for o in let_out)]:
            output_dict[f]=datasets[f]
        return output_dict

################################################
#############       statistics and analysis

def quantileCategory(s1=pd.Series(), q=4):
    '''Quantile-based discretization function based on pandas.qcut function.

    :param s1: Series()
    :param q: int() number of quantiles, default q=4
    :return: 10 for deciles, 5 for quantiles, 4 for quartiles, etc.
    '''

    temp_df = pd.DataFrame()
    temp_df['data'] = s1
    temp_df['quantiles'] = 200
    quantiles = pd.qcut(s1, q, retbins=True)
    for i in range(0, len(quantiles[1]) - 1):
        temp_df['quantiles'][temp_df.data >= quantiles[1][i]] = i
    return temp_df['quantiles']

def runPCA(data=pd.DataFrame(), n_components=2):
    # x = StandardScaler().fit_transform(df1_codone_composition)

    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(data)
    principalDf = pd.DataFrame(data=principalComponents,
                               columns=['PC' + str(i + 1) for i in np.arange(0, n_components)])

    print(pca.explained_variance_ratio_)

    finalDf = pd.concat([principalDf, pd.Series(data.index)], axis=1)
    finalDf.rename(columns={0: 'name'}, inplace=True)
    finalDf = finalDf.set_index('name')

    return finalDf


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