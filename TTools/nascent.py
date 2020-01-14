import pandas as pd
import TTools as tt
import os, shutil, random

### NASCENT (nascentFolding)

################################################
#############    preparing sequence

def extendingWindow(sequence="", name="name", strand="plus", temp=37, m=7):
    '''
    Returns DataFrame of sequences of all possible lengths between minimum (m) and length of input sequence -1
    :param sequence: str
    :param name: str, default=name
    :param strand: str {plus,minus,both}, default=plus (not tested for others)
    :param temp: int, default=37
    :param m: inf, default=7 - RNAfold does not return any values for length shorter than 7 even at 4 deg C
    :return: DataFrame of sequences with names according to nascent.slidingWindow convention
    '''
    seqList = []
    nameList = []

    for i in range(m, len(sequence)):
        seqList.append(sequence[0:i])
        nameList.append(name + "_" + strand + "_pos" + str(i) + "_temp" + str(temp) + "_win" + str(i))

    return pd.DataFrame(pd.Series(data=seqList, index=nameList, name='sequence'))

def slidingWindow(sequence="", name="name", strand="plus", temp=37, window=100):
    '''
    Slices sequence using sliding window
    :param sequence: str
    :param name: str default="name"
    :param strand: str {"plus","both","minus"} default="plus"
    :param temp: int default=37
    :param window: int default=80
    :return: DataFrame with sliding windows
    '''
    # sliding window for the sequence
    index_temp = list()
    list_temp = list()

    # plus strand
    if strand == "plus" or strand == 'both':
        for i in range(window, len(sequence) + 1):
            string = sequence[i - window:i]
            list_temp.append(string)
            index_temp.append(name + '_plus' + '_pos' + str(i) + '_temp' + str(temp) + "_win" + str(window))

    # minus strand
    if strand == "both" or strand == 'minus':
        sequence_revcomp = tt.methods.reverse_complement_RNA(sequence.replace("T", "U"))
        for i in range(window, len(sequence_revcomp) + 1):
            string_rc = sequence_revcomp[i - window:i]  # correction to get exactly the same position for both strands
            list_temp.append(string_rc)
            index_temp.append(
                name + '_minus' + '_pos' + str(str(len(sequence_revcomp) - i + 1)) + '_temp' + str(temp) + "_win" + str(
                    window))

    # output
    df_temp = pd.DataFrame()
    df_temp['sequence'] = pd.Series(list_temp, index=index_temp)
    return df_temp


def prepareNascent(sequence="", name="name", strand="plus", temp=37, window=100):
    '''
    Combines output of extendingWindow and slidingWindow
    :param sequence: str
    :param name: str, default=name
    :param strand: str {plus,minus,both}, default=plus (not tested for others)
    :param temp: int, default=37
    :param window:
    :return: Dataframe with sequences
    '''
    extended = extendingWindow(sequence=sequence[:window], name=name,
                                         strand=strand, temp=temp)
    slided = slidingWindow(sequence=sequence, name=name,
                                         strand=strand, temp=temp, window=window)
    return extended.append(slided)

def handleInput(data, keepNames=True):
    '''
    input data with columns: "seq" or "sequence" and "name" (optional)
    :param data: str() or list() or Series() or DataFrame()
    :param keepNames: default True, use given names
    :return: DataFrame() where index become name of sequence to fold
    '''

    # string to list
    if isinstance(data, str):
        data = [data]

    # list to Series
    if isinstance(data, list):
        data = pd.Series(data)

    # DataFrame as input
    if isinstance(data, pd.DataFrame):
        # check
        columns = data.columns.values.tolist()
        try:
            col_seq = [name for name in columns if 'seq' in name or 'sequence' in name][0]
        except:
            print('Could not find column with sequence "sequence" or "seq".')
            return None

        if keepNames == False:
            data = data[col_seq]

        elif keepNames == True:
            try:
                col_name = [name for name in columns if 'name' in name][0]
                data = data[[col_seq, col_name]]
                data = data.set_index(col_name)
            except:
                data.index = data.index.rename('name')

    # Series to DataFrame, returns DataFrame
    if isinstance(data, pd.Series):
        if keepNames == True:
            data = pd.DataFrame(data.rename('sequence'))
            data.index = data.index.rename('name')
        elif keepNames == False:
            data = data.rename('sequence').reset_index()
            data['name'] = "seq_" + data['index'].astype(str)
            data = data.drop('index', axis=1)
            data = data.set_index('name')

    return data

################################################
#############        folding

class Fold:
    def __init__(self, tempDir=None):
        if not tempDir:
            self.tempDir = os.getenv("HOME") + '/nascentFolding/temp/'
        else:
            self.tempDir = tempDir

        self.temp = 30
        self.Na = 0.06  # for UNAfold
        self.Mg = 0.008  # for UNAfold

    def bashFolding(self, method="RNA"):
        '''
        Runs RNA folding using bash
        :param method: "RNA" for ViennaRNA or "UNA" for UNAfold"
        :return:
        '''
        ''''''
        tt.methods.bashCommand("tab2fasta.awk temp.tab > temp.fasta")
        if method == "UNA" or method == 'both':
            tt.methods.bashCommand("hybrid-ss-min --tmin=" + str(self.temp) + " --tmax=" + str(self.temp) + " --magnesium=" + str(
                self.Mg) + " --sodium=" + str(self.Na) + " temp.fasta")
        elif method == "RNA" or method == 'both':
            tt.methods.bashCommand("RNAfold -T " + str(
                self.temp) + " -i temp.fasta | sed 's/( -/(-/g' | awk 'BEGIN{RS=\">\"}{print $1\"\t\"$2\"\t\"$3\"\t\"$4}' > RNAfold.output")

    def UNAfold(self, data, saveData=False, temp=None):
        '''
        Calculates dG using UNAfold
        :param data: input data (list, Series or DataFrame with "name" column)
        :param saveData: boolean() default=False
        :param temp: int() default=None
        :return: DataFrame()
        '''

        # check input
        data = handleInput(data)

        # prepare folder
        directory = self.tempDir + tt.methods.timestampRandomInt() + "_" + str(self.temp) + "_UNAfold"
        os.makedirs(directory)
        path = directory + "/"

        # folding
        if temp: self.temp = temp
        data_temp = data
        data_temp.to_csv(path + "temp.tab", sep="\t", header=False)
        os.chdir(path)
        self.bashFolding(method="UNA")
        df_temp = pd.read_csv(path + "temp.fasta.dG", sep="\t")

        # clean folder
        if saveData == False:
            os.chdir(self.tempDir)
            shutil.rmtree(directory)

        # return dG
        data['dG'] = df_temp['-RT ln Z']
        return data

    def RNAfold(self, data, saveData=False, temp=None):
        '''
        Calculates dG using RNAfold (ViennaRNA)
        :param data: input data (list, Series or DataFrame with "name" column)
        :param saveData: boolean() default=False
        :param temp: int() default=None
        :return: DataFrame()
        '''

        # check input
        data = handleInput(data)

        # prepare folder
        directory = self.tempDir + tt.methods.timestampRandomInt() + "_" + str(self.temp) + "_RNAfold"
        os.makedirs(directory)
        path = directory + "/"

        # folding
        if temp: self.temp = temp
        data_temp = data
        data_temp.to_csv(path + "temp.tab", sep="\t", header=False)
        os.chdir(path)
        self.bashFolding(method="RNA")
        df_temp = pd.read_csv(path + "RNAfold.output", sep="\t", names=['name', 'sequence', 'vienna', 'dG'], header=0)

        # clean folder
        if saveData == False:
            os.chdir(self.tempDir)
            shutil.rmtree(directory)

        # return dG
        df_temp['dG'] = df_temp['dG'].str.replace("(","0.0)")       # trick to avoid crushing for 0.0
        df_temp['dG'] = df_temp['dG'].str.replace("0.0\)-", "-")    # not very elegant but works
        df_temp['dG'] = df_temp['dG'].str.strip(")").astype(float)
        return df_temp

    def RNAinvert(self, structure="", saveData=False, temp=None, n=5, RNAprimer="", stall="", quick=False):
        '''Returns 2*n sequences with with given structure

        :param structure: str() with secondary RNA structure
        :param saveData: boolean() default=False
        :param temp: int() default=None
        :param n: int()
        :param RNAprimer: str() sequence
        :param stall: str() sequence
        :param quick: boolean() if False uses -Fmp -f 0.01 params
        :return: DataFrame()
        '''

        # prepare folder
        directory = self.tempDir + tt.methods.timestampRandomInt() + "_RNAinvert"
        os.mkdir(directory)
        path = directory + "/"

        # folding
        if temp: self.temp = temp
        os.chdir(path)
        tt.assays.structureFile(structure=structure, stall=stall,
                      RNAprimer=RNAprimer)  # saves structure to file together with sequence constrain
        # run RNAinvert
        if quick == False:
            tt.methods.bashCommand("cat structure.txt | RNAinverse -R" + str(n) + " -Fmp -f 0.01 -T " + str(
                self.temp) + " | awk '{print $1}' > output.tab")
        elif quick == True:
            tt.methods.bashCommand("cat structure.txt | RNAinverse -R" + str(n) + " -T " + str(
                self.temp) + " | awk '{print $1}' > output.tab")

        df_temp = pd.read_csv(path + "output.tab", sep="\t", names=['sequence'])

        # clean folder
        if saveData == False:
            os.chdir(self.tempDir)
            shutil.rmtree(directory)

        # return dG
        return df_temp

    def RNAinvertStall(self, structure="", RNAprimer="", stall="A", n=200):
        '''
        Returns sequences without nt that is present in stall
        :param structure: str() with secondary RNA structure
        :param RNAprimer: str() sequence
        :param stall: str() sequence
        :param n: int()
        :return: DataFrame()
        '''

        # generates sequences according to the structure, quick parameter uses faster algorithm
        df_sequences = self.RNAinvert(structure=structure, RNAprimer=RNAprimer, n=n, quick=True)

        # filters out sequences with stall sequence outside the primer
        df_filtered = df_sequences[df_sequences['sequence'].str.strip(RNAprimer.lower()).str.find(stall) == -1]

        return df_filtered.reset_index().drop('index', axis=1)

################################################
#############        output

### Vienna ###

def analyseViennamarkGC(vienna=str(),sequence=str()):
    '''
    leaves only C and G in stem structures
    :param vienna: str() with vienna format
    :param sequence: str() sequence
    :return:
    '''

    output = []
    for i, l in enumerate(vienna):
        if l in ["(", ")"] and sequence[i] in ["C","G"]:
            output.append(sequence[i])
#         elif l in ["(", ")"] and sequence[i] not in ["C","G"]:
#             output.append(l)
        else: output.append(".")
    return "".join(output)

def markVienna(df=pd.DataFrame()):
    '''
    :param df: DataFrame()
    :return: DataFrame()
    '''
    df['marked'] = pd.Series()
    for i,row in df.iterrows():
        df.loc[i, 'marked'] = analyseViennamarkGC(row['vienna'],row['sequence'])
#     return df.apply(lambda x: analyseVienna(df['vienna'],df['sequence']), axis=1)
    return df

def selectFoldedN(data=pd.DataFrame(), n=5, pattern="(((((....)))))"):
    '''
    Takes Fold().RNAFold() df as an input. Selects n rows with a given pattern
    on the 5end and most different folding energy
    :param data: DataFrame()
    :param n: int() samples
    :param pattern: str() vienna format
    :return: DataFrame()
    '''

    n = n  # row every n
    step = len(data) / n  # step size

    df_output = pd.DataFrame()

    take = True
    try:
        for i, row in data.sort_values('dG').reset_index().iterrows():

            if take == True:
                if row['vienna'].startswith(pattern):
                    take = False
                    df_output = df_output.append(row)

            if ((i + 1) % step) == 0:
                take = True

        return df_output.drop('index', axis=1).reset_index().drop('index', axis=1)

    except:
        print("Does your input df is Fold().RNAfold() output?")

### sliding window

def join2d(df=pd.DataFrame(), use='format'):
    # uses output of nf.Fold().RNAfold(df)
    # plus strand only

    # parse names
    df_names = df['name'].str.split("_", expand=True)
    df_names.columns = ['name', 'strand', 'position', 'temp', 'window']
    # TODO minus strand
    df_names['position'] = df_names['position'].str.replace('pos', '').astype(int)
    df_names['temp'] = df_names['temp'].str.replace('temp', '').astype(int)
    df_names['window'] = df_names['window'].str.replace('win', '').astype(int)

    df = pd.concat([df_names, df.drop('name', axis=1)], axis=1)

    # prepares output - dictionary of DataFrames - one for each gene/chromosome (name)
    output = {}

    # separate by gene/chromosome
    for name, d1 in df.groupby('name'):

        # separate by strand
        for strand, d2 in d1.groupby('strand'):
            if strand == 'plus':
                l = d2.max()['position']
                df_out = pd.DataFrame(index=range(1, l))

                for i, fold in d2.iterrows():
                    stop = fold['position']
                    start = stop - fold['window']

                    df_out[i] = pd.Series(index=range(start + 1, stop + 1), data=list(fold[use]))

            elif strand == "minus":
                print("Minus strand not supported.")

        output[name] = df_out

        # output
    return output


def merge2d(df=pd.DataFrame):
    stat = df.T.describe(include='object').T
    return stat


def nascentElem(vienna="", sequence="", verbose=False):
    l = len(sequence)
    if len(vienna) != l: return False

    # extract positions of stems
    stems, mstems = tt.secondary.loopStems(vienna)

    # choose the last nt of the last substructure
    if mstems:
        last = max(max(max(stems)), max(mstems))  # in ofder: tuple, list, one of two lists
    elif stems:
        last = max(max(stems))
    else:
        return None, None

    if last in mstems:
        return sequence[mstems[0] - 1:last], l - last

    else:
        for s in stems:
            if last == s[1]:
                return sequence[s[0] - 1:s[1]], l - last


def nascentElemDataFrame(data=pd.DataFrame()):
    data = data.copy()
    indexList = []
    elemSequenceList = []
    elemDistanceList = []

    for i, df in data.iterrows():
        indexList.append(i)
        seq, dist = nascentElem(df['vienna'], df['sequence'])
        elemSequenceList.append(seq)
        elemDistanceList.append(dist)

    data['elemSequence'] = pd.Series(data=elemSequenceList, index=indexList)
    data['elemDistance'] = pd.Series(data=elemDistanceList, index=indexList)

    return data


def foldNascentElem(data=pd.DataFrame()):
    # parse names and define folding temperature
    df_names = data['name'].str.split("_", expand=True)
    df_names.columns = ['name', 'strand', 'position', 'temp', 'window']
    # TODO minus strand
    df_names['temp'] = df_names['temp'].str.replace('temp', '').astype(int)
    temps = df_names['temp'].unique()
    if len(temps) > 1:
        print("More than 1 folding temperature detected")
        return False
    elif len(temps) < 1:
        print("No folding temperature detected")
        return False
    else:
        temp = temps[0]

    # folding
    toFold = data[~data['elemSequence'].isna()]['elemSequence'].unique().tolist()

    folded = Fold().RNAfold(toFold, temp=temp)

    energies_dict = folded[['sequence', 'dG']].set_index('sequence').to_dict()['dG']
    energies_dict[None] = 0.0

    output = data.copy()
    output['elemDistance'] = output['elemDistance'].fillna(0).astype(int)
    output['elem_dG'] = output['elemSequence'].map(energies_dict)

    return output


def nascentFolding(sequence='', name='name', temp=37, window=100):
    # preparing sequence sliding window
    df = prepareNascent(sequence, temp=temp, window=window)
    # folding arbitraty windows
    df = Fold().RNAfold(df, temp=temp)
    # defining nascent element and distance
    df = nascentElemDataFrame(data=df)
    # folding nascent element
    df = foldNascentElem(df)

    # parse names
    df_names = df['name'].str.split("_", expand=True)
    df_names.columns = ['name', 'strand', 'position', 'temp', 'window']
    # TODO minus strand
    df_names['position'] = df_names['position'].str.replace('pos', '').astype(int)
    df_names['temp'] = df_names['temp'].str.replace('temp', '').astype(int)
    df_names['window'] = df_names['window'].str.replace('win', '').astype(int)

    df = pd.concat([df_names, df.drop('name', axis=1)], axis=1)

    # preparing output
    zeros = df['position'].min()
    zeros = pd.DataFrame(data={"elemDistance": [0] * (zeros - 1),
                               "elem_dG": [0.0] * (zeros - 1),
                               'position': range(1, zeros)}).set_index('position')

    df = df[['position', 'elemDistance', 'elem_dG']].set_index('position')

    return pd.concat([zeros, df])