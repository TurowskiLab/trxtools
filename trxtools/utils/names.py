import sys, re
import pandas as pd



################################################
#############        handling names for multiple experiments

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