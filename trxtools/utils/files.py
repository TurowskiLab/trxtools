import time, random
import os
import pandas as pd



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

def gtf2bed(gtf_file):
    '''Convert GTF file to BED format
    :param gtf_file: Path to the GTF file
    :type gtf_file: str
    
    :return: DataFrame in BED format
    :rtype: pandas.DataFrame'''

    gtf = pd.read_csv(gtf_file, sep="\t", header=None)
    gtf = gtf[gtf[2].isin(['gene'])]

    gtf = pd.concat([
        gtf.rename(columns={0: 'chr', 1: 'source', 2:'feature_type'}),
        gtf[8].str.split("; ", expand=True)[[0,1,2]]
        ], axis=1).rename(columns={0: 'gene_id', 1: 'gene_type', 2:'gene_name'})
    gtf = gtf.drop(columns=[8])
    gtf['gene_id'] = gtf['gene_id'].str.replace('gene_id "', '', regex=False)
    gtf['gene_id'] = gtf['gene_id'].str.replace('"', '', regex=False)
    gtf['gene_type'] = gtf['gene_type'].str.replace('gene_type "', '', regex=False)
    gtf['gene_type'] = gtf['gene_type'].str.replace('"', '', regex=False)
    gtf['gene_name'] = gtf['gene_name'].str.replace('gene_name "', '', regex=False)
    gtf['gene_name'] = gtf['gene_name'].str.replace('"', '', regex=False)

    bed = gtf[['chr', 3, 4, 'gene_id', 5, 6, 'gene_type', 'gene_name']]
    bed.columns = [0, 1, 2, 3, 4, 5, 'gene_type', 'gene_name']
    
    return bed

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
