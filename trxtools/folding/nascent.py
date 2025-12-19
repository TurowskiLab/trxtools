import pandas as pd
import trxtools as tt
import os, shutil, random

### NASCENT (nascentFolding)

################################################
#############    preparing sequence
def extendingWindow(sequence="", name="name", strand="plus", temp=30, m=7, toAdd=0):
    '''Returns DataFrame of sequences of all possible lengths between minimum (m) and length of input sequence -1.

    :param sequence: str, input sequence
    :type sequence: str
    :param name: name of the sequence, default "name"
    :type name: str
    :param strand: strand type {"plus", "minus"}, default "plus"
    :type strand: str
    :param temp: temperature, default 30
    :type temp: int
    :param m: minimum length of sequence, default 7
    :type m: int
    :param toAdd: additional length to add, default 0
    :type toAdd: int

    :return: DataFrame of sequences with names according to nascent.slidingWindow convention
    :rtype: pd.DataFrame
    '''

    seqList = []
    nameList = []
    if strand == "minus":
        sequence = tt.methods.reverse_complement(sequence)
        for i in range(m, len(sequence)):
            seqList.append(sequence[0:i])
            nameList.append(name + "_" + strand + "_pos" + str(toAdd-i+1) + "_temp" + str(temp) + "_win" + str(i))
    else:
        for i in range(m, len(sequence)):
            seqList.append(sequence[0:i])
            nameList.append(name + "_" + strand + "_pos" + str(toAdd+i) + "_temp" + str(temp) + "_win" + str(i))
    return pd.DataFrame(pd.Series(data=seqList, index=nameList, name='sequence'))

def slidingWindow(sequence="", name="name", strand="plus", temp=30, window=100):
    '''Slices sequence using sliding window.

    :param sequence: input sequence
    :type sequence: str
    :param name: name of the sequence, default "name"
    :type name: str
    :param strand: strand type {"plus", "both", "minus"}, default "plus"
    :type strand: str
    :param temp: temperature, default 30
    :type temp: int
    :param window: window size, default 100
    :type window: int

    :return: DataFrame with sliding windows
    :rtype: pd.DataFrame
    '''

    index_temp = list()
    list_temp = list()

    if strand == "plus" or strand == 'both':
        for i in range(window, len(sequence) + 1):
            string = sequence[i - window:i]
            list_temp.append(string)
            index_temp.append(name + '_plus' + '_pos' + str(i) + '_temp' + str(temp) + "_win" + str(window))

    if strand == "both" or strand == 'minus':
        sequence_revcomp = tt.methods.reverse_complement_DNA(sequence)
        for i in range(window, len(sequence_revcomp) + 1):
            string_rc = sequence_revcomp[i - window:i]
            list_temp.append(string_rc)
            index_temp.append(name + '_minus' + '_pos' + str(len(sequence_revcomp) - i + 1) + '_temp' + str(temp) + "_win" + str(window))

    df_temp = pd.DataFrame()
    df_temp['sequence'] = pd.Series(list_temp, index=index_temp)
    return df_temp

def prepareNascent(sequence="", name="name", strand="plus", temp=30, window=100):
    '''Divide long transcript into short sequences. Combines output of extendingWindow and slidingWindow.

    :param sequence: input sequence
    :type sequence: str
    :param name: name of the sequence, default "name"
    :type name: str
    :param strand: strand type {"plus", "minus"}, default "plus"
    :type strand: str
    :param temp: temperature, default 30
    :type temp: int
    :param window: window size, default 100
    :type window: int

    :return: DataFrame with sequences
    :rtype: pd.DataFrame
    '''

    if strand == "minus":
        extended = extendingWindow(sequence=sequence[-window:], name=name, strand=strand, temp=temp, toAdd=len(sequence))
    else:
        extended = extendingWindow(sequence=sequence[:window], name=name, strand=strand, temp=temp)
    slided = slidingWindow(sequence=sequence, name=name, strand=strand, temp=temp, window=window)
    return pd.concat([extended, slided])

def name2index(s1=pd.Series(dtype=str)):
    '''Extracts position from sequence name.

    :param s1: Series with names from prepareNascent function
    :type s1: pd.Series

    :return: Series with positions
    :rtype: pd.Series
    '''

    return s1.str.split("_", expand=True)[2].str.replace("pos", "").astype(int)

def handleInput(data, keepNames=True):
    '''Input data with columns: "seq" or "sequence" and "name" (optional).

    :param data: input data {str, list, Series, DataFrame}
    :type data: {str, list, pd.Series, pd.DataFrame}
    :param keepNames: if True use given names, default True
    :type keepNames: bool

    :return: DataFrame where index become name of sequence to fold
    :rtype: pd.DataFrame
    '''

    if isinstance(data, str):
        data = [data]

    if isinstance(data, list):
        data = pd.Series(data)

    if isinstance(data, pd.DataFrame):
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
        self.Na = 0.03  # for UNAfold (30 mM in yeast, 10 mM in humans)
        self.Mg = 0.008  # for UNAfold

    def bashFolding(self, method="RNA"):
        '''Runs RNA folding using bash.

        :param method: Folding method, "RNA" for ViennaRNA or "UNA" for UNAfold.
        :type method: str

        :return: None
        '''

        tt.methods.bashCommand("tab2fasta.awk temp.tab > temp.fasta")
        if method == "UNA" or method == 'both':
            tt.methods.bashCommand("hybrid-ss-min --tmin=" + str(self.temp) + " --tmax=" + str(self.temp) + " --magnesium=" + str(
                self.Mg) + " --sodium=" + str(self.Na) + " temp.fasta")
        elif method == "RNA" or method == 'both':
            tt.methods.bashCommand("RNAfold -T " + str(
                self.temp) + " -i temp.fasta | sed 's/( -/(-/g' | awk 'BEGIN{RS=\">\"}{print $1\"\t\"$2\"\t\"$3\"\t\"$4}' > RNAfold.output")

    def UNAfold(self, data, saveData=False, temp=None):
        '''Calculates dG using UNAfold.

        :param data: Input data, can be a list, Series, or DataFrame with "name" column.
        :type data: {list, pd.Series, pd.DataFrame}
        :param saveData: If True, saves the data, default is False.
        :type saveData: bool
        :param temp: Temperature for folding, default is None.
        :type temp: int

        :return: DataFrame with calculated dG values.
        :rtype: pd.DataFrame
        '''

        # check input
        data = handleInput(data)

        # prepare folder
        if temp: self.temp = temp
        directory = self.tempDir + tt.methods.timestampRandomInt() + "_" + str(self.temp) + "_UNAfold"
        os.makedirs(directory)
        path = directory + "/"

        # folding
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

    # # # # #
    def RNAfold(self, data, saveData=False, temp=None):
        '''Calculates dG using RNAfold (ViennaRNA).

        :param data: Input data, can be a list, Series, or DataFrame with "name" column.
        :type data: {list, pd.Series, pd.DataFrame}
        :param saveData: If True, saves the data, default is False.
        :type saveData: bool
        :param temp: Temperature for folding, default is None.
        :type temp: int

        :return: DataFrame with calculated dG values.
        :rtype: pd.DataFrame
        '''

        # check input
        data = handleInput(data)

        # prepare folder
        if temp: self.temp = temp
        directory = self.tempDir + tt.methods.timestampRandomInt() + "_" + str(self.temp) + "_RNAfold"
        os.makedirs(directory)
        path = directory + "/"

        # folding
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
        df_temp['dG'] = df_temp['dG'].str.replace("(","0.0)")       # trick to avoid crashing for 0.0
        df_temp['dG'] = df_temp['dG'].str.replace("0.0)-", "-")     # not very elegant but works
        df_temp['dG'] = df_temp['dG'].str.strip(")").astype(float)
        return df_temp

    def RNAinvert(self, structure="", saveData=False, temp=None, n=5, RNAprimer="", stall="", quick=False):
        '''Returns n sequences with a given structure.

        :param structure: Secondary RNA structure.
        :type structure: str
        :param saveData: If True, saves the data, default is False.
        :type saveData: bool
        :param temp: Temperature for folding, default is None.
        :type temp: int
        :param n: Number of output sequences, default is 5.
        :type n: int
        :param RNAprimer: RNA primer sequence.
        :type RNAprimer: str
        :param stall: Nucleotide sequence to avoid.
        :type stall: str
        :param quick: If False, uses -Fmp -f 0.01 params, default is False.
        :type quick: bool

        :return: DataFrame with generated sequences.
        :rtype: pd.DataFrame
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
        '''Returns sequences without the nucleotide that is present in stall.

        :param structure: Secondary RNA structure.
        :type structure: str
        :param RNAprimer: RNA primer sequence.
        :type RNAprimer: str
        :param stall: Nucleotide sequence to avoid, default is "A".
        :type stall: str
        :param n: Number of sequences to generate, default is 200.
        :type n: int

        :return: DataFrame with filtered sequences.
        :rtype: pd.DataFrame
        '''

        # generates sequences according to the structure, quick parameter uses faster algorithm
        df_sequences = self.RNAinvert(structure=structure, RNAprimer=RNAprimer, n=n, quick=True)

        # filters out sequences with stall sequence outside the primer
        df_filtered = df_sequences[df_sequences['sequence'].str.strip(RNAprimer.lower()).str.find(stall) == -1]

        return df_filtered.reset_index().drop('index', axis=1)

################################################
#############        hybrid

class Hybrid:
    def __init__(self, tempDir=None):
        if not tempDir:
            self.tempDir = os.getenv("HOME") + '/nascentFolding/temp_hybrids/'
        else:
            self.tempDir = tempDir

        self.temp = 37
        self.n = 'RNA'  # "RNA" or "DNA"
        self.Na = 0.01  # for human cell
        self.Mg = 0.01  # for human cell

    def bashHybrid(self):
        '''Runs hybrid-min using bash.
        This function executes the hybrid-min command using bash to calculate the hybridization energy between two sequences.

        :return: None
        '''
        
        tt.methods.bashCommand("tab2fasta.awk temp1.tab > temp1.fasta")
        tt.methods.bashCommand("tab2fasta.awk temp2.tab > temp2.fasta")
        tt.methods.bashCommand("hybrid-min -T " + str(self.temp)
                               + " -n " + str(self.n)
                               + " -N " + str(self.Na)
                               + " -M " + str(self.Mg)
                               + " temp1.fasta temp2.fasta")

    def RNAhybrid(self, data, saveData=False, temp=None):
        '''Calculates dG using hybrid-min.
        This function calculates the hybridization energy (dG) between RNA sequences using the hybrid-min tool.

        :param data: Input data, can be a list, Series, or DataFrame with "name" column.
        :type data: {list, pd.Series, pd.DataFrame}
        :param saveData: If True, saves the data, default is False.
        :type saveData: bool
        :param temp: Temperature for hybridization, default is None.
        :type temp: int

        :return: DataFrame with calculated dG values.
        :rtype: pd.DataFrame

        :example:

        >>> hybrid = Hybrid()
        >>> data = ["AGCUAGUCA", "CGAUCGUAG"]
        >>> hybrid.RNAhybrid(data)
        '''

        # check input
        data = handleInput(data)

        # prepare folder
        directory = self.tempDir + tt.methods.timestampRandomInt() + "_" + str(self.temp) + "_hybrid_min"
        os.makedirs(directory)
        path = directory + "/"

        # folding
        if temp: self.temp = temp
        data_temp1 = data
        data_temp1.to_csv(path + "temp1.tab", sep="\t", header=False)
        # reverse complement
        data_temp2 = pd.DataFrame(data_temp1['sequence'].apply(tt.methods.reverse_complement))
        data_temp2.to_csv(path + "temp2.tab", sep="\t", header=False)

        os.chdir(path)
        self.bashHybrid()
        df_temp = pd.read_csv(path + "temp1.fasta-temp2.fasta.dG", sep="\t", names=['temp', 'dG', 'Z'], header=0)

        # clean folder
        if saveData == False:
            os.chdir(self.tempDir)
            shutil.rmtree(directory)

        # return dG
        df_temp['dG'] = df_temp['dG'].astype(float)
        return df_temp['dG']

### sliding window

def join2d(df=pd.DataFrame(), use='format'):
    '''Joins 2D data from RNA folding results.

    :param df: DataFrame containing RNA folding results.
    :type df: pd.DataFrame
    :param use: Column to use for joining, default is 'format'.
    :type use: str

    :return: Dictionary of DataFrames, one for each gene/chromosome.
    :rtype: dict
    '''
    
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
    '''Merges 2D data from RNA folding results.

    :param df: DataFrame containing RNA folding results.
    :type df: pd.DataFrame

    :return: DataFrame with statistical description of the merged data.
    :rtype: pd.DataFrame
    '''

    stat = df.T.describe(include='object').T
    return stat

def nascentElems(vienna="", sequence=""):
    '''Describe elements of secondary structure: stems, multistems.

    :param vienna: Vienna RNA secondary structure notation.
    :type vienna: str
    :param sequence: RNA sequence.
    :type sequence: str

    :return: Tuple containing the sequence of the last nascent element and its distance from the end.
    :rtype: tuple
    '''

    l = len(sequence)
    if len(vienna) != l: return False

    # extract positions of stems
    stems, mstems = tt.secondary.loopStems(vienna)

    # choose the last nt of the last substructure
    if mstems:
        last = max(max(max(stems)), max(mstems))  # in order: tuple, list, one of two lists
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


def nascentElemsDataFrame(data=pd.DataFrame()):
    '''Apply nascentElems to DataFrame of folded sequences.

    :param data: DataFrame containing folded sequences with 'vienna' and 'sequence' columns.
    :type data: pd.DataFrame

    :return: DataFrame with additional columns for nascent element sequences and distances.
    :rtype: pd.DataFrame
    '''

    data = data.copy()
    indexList = []
    elemSequenceList = []
    elemDistanceList = []

    for i, df in data.iterrows():
        indexList.append(i)
        seq, dist = nascentElems(df['vienna'], df['sequence'])
        elemSequenceList.append(seq)
        elemDistanceList.append(dist)

    data['elemSequence'] = pd.Series(data=elemSequenceList, index=indexList)
    data['elemDistance'] = pd.Series(data=elemDistanceList, index=indexList)

    return data

def foldNascentElem(data=pd.DataFrame()):
    '''Fold the very 3' of nascent elements.

    :param data: DataFrame containing folded sequences.
    :type data: pd.DataFrame

    :return: DataFrame with additional columns for nascent element folding energies.
    :rtype: pd.DataFrame
    '''

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

def parseFoldingName(df=pd.DataFrame()):
    '''Parse folding names into separate columns.

    :param df: DataFrame containing folding names.
    :type df: pd.DataFrame

    :return: DataFrame with parsed columns.
    :rtype: pd.DataFrame
    '''

    # parse names
    df_names = df['name'].str.split("_", expand=True)
    df_names.columns = ['name', 'strand', 'position', 'temp', 'window']
    # TODO minus strand
    df_names['position'] = df_names['position'].str.replace('pos', '').astype(int)
    df_names['temp'] = df_names['temp'].str.replace('temp', '').astype(int)
    df_names['window'] = df_names['window'].str.replace('win', '').astype(int)

    return pd.concat([df_names, df.drop('name', axis=1)], axis=1)

def nascentFolding(sequence='', temp=30, window=100):
    '''Combines folding function: fold RNA, locate last nascent element and calculate dG of it.

    :param sequence: RNA sequence.
    :type sequence: str
    :param temp: Folding temperature, default is 30.
    :type temp: int
    :param window: Window size for folding, default is 100.
    :type window: int
    
    :return: DataFrame with folding results and nascent element information.
    :rtype: pd.DataFrame
    '''

    # preparing sequence sliding window
    df = prepareNascent(sequence, temp=temp, window=window)
    # folding arbitrary windows
    df = Fold().RNAfold(df, temp=temp)
    # defining nascent element and distance
    df = nascentElemsDataFrame(data=df)
    # folding nascent element
    df = foldNascentElem(df)
    # parse name
    df = parseFoldingName(df)

    # preparing output
    zeros = df['position'].min()
    zeros = pd.DataFrame(data={"elemDistance": [0] * (zeros - 1),
                               "elem_dG": [0.0] * (zeros - 1),
                               'position': range(1, zeros)}).set_index('position')

    df = df[['position', 'elemDistance', 'elem_dG']].set_index('position')

    return pd.concat([zeros, df])