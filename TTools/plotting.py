import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import trxtools.profiles as profiles
from adjustText import adjust_text
import seaborn as sns

#### PCA
def plotPCA(data=pd.DataFrame(), names=[], title="", PClimit=1,figsize=(7,7), PCval=[]):
    '''Plot PCA plot

    :param data: DataFrame
    :param names: list of names to annotate
    :param title: str
    :param PClimit: int number of PC to plot, default 1
    :param figsize: tuple, default (7,7)
    :return:

    >>> plotPCA(methods.runPCA(example_df)[0])
    '''
    nPCA = len([col for col in data.columns if 'PC' in col])
    axes = [ i +1 for i in range(nPCA)[:-1]]
    if not PCval:
        PCval = [""]*(PClimit+1)

    for nPC in axes:
        fig = plt.figure(figsize=figsize)
        #         ax = fig.add_subplot(nPCA-1,1,nPC)
        ax = fig.add_subplot(1 ,1 ,1)

        texts = []

        if 'group' in data.columns.tolist():
            for group, df_temp in data.groupby('group'):
                a = df_temp['PC' +str(nPC)].tolist()
                b = df_temp['PC' +str(nPC+1)].tolist()
                ax.scatter(x=a ,y=b ,label=group, cmap="paired")
                for x, y, s in zip(a, b, df_temp.index.tolist()):
                    if s in names:
                        texts.append(plt.text(x, y, s))

        else:
            a = data['PC' +str(nPC)].tolist()
            b = data['PC' +str(nPC+1)].tolist()
            ax.scatter(x=a ,y=b ,color='lightgray')
            for x, y, s in zip(a, b, data.index.tolist()):
                if s in names:
                    texts.append(plt.text(x, y, s))

        ax.legend()
        ax.grid(True,ls="dotted")
        plt.xlabel('PC ' +str(nPC)+" ("+str(PCval[nPC-1])+"%)")
        plt.ylabel('PC ' +str(nPC +1)+" ("+str(PCval[nPC])+"%)")
        plt.title(title)
        adjust_text(texts, only_move={'points':'y', 'texts':'y'}, arrowprops=dict(arrowstyle="->", color='black', lw=0.5))
        plt.show()

        if nPC==PClimit: break

def clusterClusterMap(df):
    '''
    Clustermap for clusters
    :param df: DataFrame
    :return:
    '''
    df_mean = pd.DataFrame()

    n = df['cluster'].max()+1

    for c in range(0, n):
        a = df[df['cluster'] == c]
        #         print("Cluster "+str(c)+"\t"+str(len(a)))
        df_mean['cluster ' + str(c) + " (" + str(len(a)) + ")"] = a.mean()

    df_mean = df_mean.T.drop('cluster', axis=1)

    sns.clustermap(df_mean, center=0, cmap="vlag",
                   linewidths=.75, figsize=(7, 7))
    plt.show()

### Peaks metaplot
def plotCumulativePeaks(ref, df2=pd.DataFrame(), local_pos=list(), dpi=150,
                        title="", start=None, stop=None, window=50, figsize=(4,3),
                        color1='green', color2="magenta", lc='red'):
    '''Plot single gene peaks metaplot.

    :param ref: str with path to csv file or DataFrame
    :param df2: DataFrame
    :param local_pos: list of features (peaks/troughs)
    :param dpi: int, default 150
    :param title: str
    :param start: int
    :param stop: int
    :param window: int, default 50
    :param figsize: tuple, default (4,3)
    :param color1: str, default "green"
    :param color2: str, default "magenta"
    :param lc: str, default "red"
    :return:
    '''
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    # extracting data for metaplot
    df_dataset1 = pd.DataFrame()
    df_dataset2 = pd.DataFrame()

    # filter local max accorgind to start and stop
    if start: local_pos = [i for i in local_pos if i > start]
    if stop: local_pos = [i for i in local_pos if i < stop]

    for i, loc in enumerate(local_pos):
        s2_dataset1 = reference['median'][loc - window:loc + window]
        df_dataset1[i] = s2_dataset1.reset_index()['median']
        s3_dataset2 = df2['median'][loc - window:loc + window]
        df_dataset2[i] = s3_dataset2.reset_index()['median']

    s_data1 = df_dataset1.mean(axis=1)
    s_data2 = df_dataset2.mean(axis=1)

    # plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    #     ax1.set_ylim(ylim)
    ax1.plot(np.arange(-window, window), s_data1, color=color1)

    # plotting dataset2
    ax1.plot(np.arange(-window, window), s_data2, color=color2)

    ax1.axvline(0, color=lc, alpha=0.5)
    ax1.legend(loc=2)
    plt.show()

### profiles
def plot_as_box_plot(df=pd.DataFrame(),title="", start=None, stop=None,
                     figsize=(7,3),ylim=(None,0.01), dpi=150, color='green',
                     h_lines=[], lc="red",offset=0):
    '''Plots figure similar to box plot: median, 2 and 3 quartiles and min-max range

    :param df: Dataframe() containing following columns:```['position'] ['mean'] ['median'] ['std']```
        optionally ```['nucleotide'] ['q1'] ['q3'] ['max'] ['min']```
    :param title: str
    :param start: int
    :param stop: int
    :param figsize: tuple, default (7,4)
    :param ylim: tuple OY axes lim. Default (None,0.01)
    :param dpi: int, default 150
    :param color: str, default "green"
    :param h_lines: list of horizontal lines
    :param lc: str color of horizontal lines, default "red"
    :param offset: int number to offset position if 5' flank was used, default 0
    :return:
    '''

    if 'nucleotide' in df.columns.values:
        df = df.drop('nucleotide', 1)
    s2 = df[start:stop]
    #plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    ax1.set_ylim(ylim)
    ax1.plot(s2.index-offset, s2['median'], color=color)
    if set(['q1','q3']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], label='range (2nd-3rd quartile)', color=color, alpha=0.2)
    if set(['min','max']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)
    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color=lc)
    ax1.legend()

def plot_to_compare(ref, df=pd.DataFrame(), color1='green', color2='black',
                    ref_label="", label="", title="", start=None, stop=None, figsize=(7,3),
                    ylim=(None,0.01), h_lines=[], lc="red", dpi=150,offset=300):
    '''Figure to compare to plots similar to box plot: median, 2 and 3 quartiles and min-max range

    :param ref: str with path to csv file or DataFrame
    :param df: DataFrame
    :param color1: str, default "green"
    :param color2: str, default "black"
    :param ref_label: str
    :param label: str
    :param title: str
    :param start: int
    :param stop: int
    :param figsize: tuple, default (7,4)
    :param ylim: tuple OY axes lim. Default (None,0.01)
    :param h_lines: list of horizontal lines
    :param lc: str color of horizontal lines, default "red"
    :param dpi: int, default 150
    :param offset: int number to offset position if 5' flank was used, default 0
    :return:
    '''

    #handling reference plot
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    dataset, s2 = df[start:stop], reference[start:stop] # prepating datasets
    #plotting
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    plt.axhline(0, color='red')

    ## plot to compare
    if len(dataset.columns) == 4: #if only two experiments
        ax1.plot(dataset.index-offset, dataset['mean'], color2, label=label)
        ax1.fill_between(dataset.index-offset, dataset['min'], dataset['max'], color=color2, alpha=0.3, label='range (min-max)')
    else: #if more than two experiments
        ax1.plot(dataset.index-offset, dataset['median'], color2, label=label)
        ax1.fill_between(dataset.index-offset, dataset['q1'], dataset['q3'], color=color2, alpha=0.2, label='range (2nd-3rd quartile)')
    ax1.set_xlabel('position')
    ax1.set_ylim(ylim)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('fraction of reads', color='black')
    for tl in ax1.get_yticklabels():
        tl.set_color('black')

    ## reference plot
    if len(s2.columns) == 4:  # if only two experiments
        ax1.plot(s2.index - offset, s2['mean'], color1, label=ref_label)
        ax1.fill_between(s2.index - offset, s2['min'], s2['max'], color=color1, alpha=0.07, label='range (min-max)')
    else:
        ax1.plot(s2.index-offset, s2['median'], color1, label=ref_label)
        ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], color=color1, alpha=0.2, label='range (q2-q3)')

    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color=lc)
    ax1.legend()

def plot_diff(ref, dataset=pd.DataFrame(), ranges='mm', label="", start=None, stop=None, plot_medians=True,
              plot_ranges=True, figsize=(7, 3), ylim=(None,0.01), h_lines=list(), offset=0):
    '''Plot given dataset and reference, differences are marked

    :param ref: str with path to csv file or DataFrame
    :param dataset: DataFrame containing following columns:```['position'] ['mean'] ['median'] ['std']```
        optionally ```['nucleotide'] ['q1'] ['q3'] ['max'] ['min']```
    :param ranges: str "mm" : min-max or "qq" : q1-q3
    :param label: str
    :param start: int
    :param stop: int
    :param plot_medians: boolean if True plot medians, default True
    :param plot_ranges: boolean if True plot ranges, default True
    :param figsize: tuple, default (7,3)
    :param ylim: tuple OY axes lim, default (None,0.01)
    :param h_lines: list of horizontal lines
    :return:
    '''

    #handling reference plot
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    ranges_dict = {'mm': 'min-max', 'qq': 'q1-q3'}

    if isinstance(reference, str):
        reference = pd.read_csv(reference, index_col=0)

    differences_df = profiles.compareMoretoRef(dataset=dataset, ranges=ranges, reference=reference)[start:stop]
    dataset, s2 = dataset[start:stop], reference[start:stop]  # prepating datasets
    # plotting
    fig, ax1 = plt.subplots(figsize=figsize)
    ax1.fill_between(differences_df.index-offset, differences_df['ear_min'], differences_df['ear_max'], color='red',
                     where=(differences_df['ear_max'] > 0), label='increased pausing (' + ranges_dict[ranges] + ')')
    ax1.fill_between(differences_df.index-offset, differences_df['rae_min'], differences_df['rae_max'], color='blue',
                     where=(differences_df['rae_max'] > 0), label='decreased pausing (' + ranges_dict[ranges] + ')')
    if plot_medians == True:
        ax1.plot(dataset.index-offset, dataset['median'], 'black', label=label)
        ax1.plot(s2.index-offset, s2['median'], 'green', label='reference RDN37-1')
    if plot_ranges == True:
        if len(dataset.columns) == 4:  # if only two experiments
            ax1.fill_between(dataset.index-offset, dataset['min'], dataset['max'], color='black', alpha=0.3, label='min-max')
        else:  # if more than two experiments
            ax1.fill_between(dataset.index-offset, dataset['q1'], dataset['q3'], color='black', alpha=0.2, label='q1-q3')
            ax1.fill_between(dataset.index-offset, dataset['min'], dataset['max'], color='black', alpha=0.07, label='min=max')
        ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], color='green', alpha=0.2, label='q1-q3')
        ax1.fill_between(s2.index-offset, s2['min'], s2['max'], color='green', alpha=0.07, label='min=max')
    ax1.set_ylim(ylim)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads ' + label, color='black')
    for i in [i for i in h_lines if i in range(start, stop)]: ax1.axvline(i, color='red')
    plt.legend()

def plot_heatmap(df=pd.DataFrame(), title='Heatmap of differences between dataset and reference plot for RDN37-1',
                 vmin=None, vmax=None, figsize=(20,10)):
    '''Plot heat map of differences, from dataframe generated by compare1toRef(dataset, heatmap=True) function

    :param df: DataFrame
    :param title: str
    :param vmin:
    :param vmax:
    :param figsize: tuple, default (20,10)
    :return:
    '''

    fig, ax = plt.subplots(figsize=figsize)
    if not vmin:
        vmin = -np.absolute(df.max().median())
    if not vmax:
        vmax = df.max().median()
    heatmap = ax.pcolormesh(df.transpose(), cmap='seismic', vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(list(df.columns.values))) + 0.5, minor=False)
    ax.set_yticklabels(list(df.columns.values), minor=False)
    fig.colorbar(heatmap)
    ax.set_title(title)

################################################
#############       STAR mapping statistics

def plotSTARstats_reads(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    reads = df['                          Number of input reads |'].astype(int)
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3

    fig, ax1 = plt.subplots(figsize=(len(labels) / 4, 3), dpi=dpi)
    plt.title('Number of reads after pre-processing')
    ax1.bar(x, reads, width)
    ax1.set_yscale('log')
    ax1.set_ylabel('Number of reads [log10]')
    ax1.grid(axis='y', c='black', ls="dotted")
    plt.xticks(x, labels, rotation=90)
    plt.show()


def plotSTARstats_readLen(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    # length of reads
    readsLen = df['                      Average input read length |'].astype(int)

    fig, ax1 = plt.subplots(figsize=(len(labels) / 4, 3), dpi=dpi)
    plt.title('Average input read length')
    ax1.bar(x, readsLen, width)
    ax1.set_ylabel('Average input read length [nt]')
    ax1.grid(axis='y', c='black', ls="dotted")
    ax1.set_ylim(0, 175)
    plt.xticks(x, labels, rotation=90)
    plt.show()


def plotSTARstats_mapping(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    # reads mapping
    mapped_uniq = df['                        Uniquely mapped reads % |'].str.strip("%").astype(float)
    mapped_multi = df['             % of reads mapped to multiple loci |'].str.strip("%").astype(float)
    mapped_tooMany = df['             % of reads mapped to too many loci |'].str.strip("%").astype(float)
    unmapped_mismatches = df['       % of reads unmapped: too many mismatches |'].str.strip("%").astype(float)
    unmapped_tooShort = df['                 % of reads unmapped: too short |'].str.strip("%").astype(float)
    unmapped_other = df['                     % of reads unmapped: other |'].str.strip("%").astype(float)

    fig, ax1 = plt.subplots(figsize=(len(labels) / 4, 3), dpi=dpi)
    plt.title('STAR mapping statistics')

    p1 = ax1.bar(x, mapped_uniq, width, color='green')
    p2 = ax1.bar(x, mapped_multi, width, bottom=mapped_uniq, color='limegreen')
    p3 = ax1.bar(x, mapped_tooMany, width, bottom=mapped_uniq + mapped_multi, color='yellow')
    p4 = ax1.bar(x, unmapped_mismatches, width, bottom=mapped_uniq + mapped_multi + mapped_tooMany, color='grey')
    p5 = ax1.bar(x, unmapped_tooShort, width, bottom=mapped_uniq + mapped_multi + mapped_tooMany + unmapped_mismatches,
                 color='orange')
    p6 = ax1.bar(x, unmapped_other, width,
                 bottom=mapped_uniq + mapped_multi + mapped_tooMany + unmapped_mismatches + unmapped_tooShort,
                 color='red')

    plt.ylabel('Cumulative mapping [%]')
    plt.xticks(x, labels, rotation=90)
    plt.yticks(np.arange(0, 101, 10))
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]), (
    'uniquely mapped', 'multimappers', 'too many mapped', 'unmapped (mistmateches)', 'unmapped (too short)',
    'unmapped (other)'))

    plt.show()


def plotSTARstats_mistmatches(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    # mistmatches and deletions
    mismatches = df['                      Mismatch rate per base, % |'].str.strip("%").astype(float)
    deletions = df['                         Deletion rate per base |'].str.strip("%").astype(float)
    insertions = df['                        Insertion rate per base |'].str.strip("%").astype(float)

    fig, ax1 = plt.subplots(figsize=(len(labels) / 4, 3), dpi=dpi)
    plt.title('Mapping errors')
    ax1.bar(x - width, mismatches, width, label='mismatches')
    ax1.bar(x, deletions, width, label='deletions')
    ax1.bar(x + width, insertions, width, label='insertions')
    ax1.grid(axis='y', c='black', ls="dotted")
    ax1.legend()
    plt.xticks(x, labels, rotation=90)
    plt.ylabel('Errors as % rate per base')
    plt.show()


def plotSTARstats_chimeric(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    reads = df['                          Number of input reads |'].astype(int)
    # Splices and chimeric reads
    splices = df['                       Number of splices: Total |'].astype(int) / reads
    chmieric = df['                            % of chimeric reads |'].str.strip("%").astype(float)

    fig, ax1 = plt.subplots(figsize=(len(labels) / 4, 3), dpi=dpi)
    plt.title('Chimeric reads')

    ax1.bar(x - width / 2, splices, width, label='splices')
    ax1.bar(x + width / 2, chmieric, width, label='chmieric')
    ax1.grid(axis='y', c='black', ls="dotted")
    ax1.legend()
    plt.xticks(x, labels, rotation=90)
    plt.ylabel('% of reads')
    plt.show()


def plotSTARstats(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    plotSTARstats_reads(df=df, dpi=dpi)
    plotSTARstats_readLen(df=df, dpi=dpi)
    plotSTARstats_mapping(df=df, dpi=dpi)
    plotSTARstats_mistmatches(df=df, dpi=dpi)
    plotSTARstats_chimeric(df=df, dpi=dpi)


def hplotSTARstats_reads(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    reads = df['                          Number of input reads |'].astype(int)
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3

    fig, ax1 = plt.subplots(figsize=(5, len(labels) / 4), dpi=dpi)
    plt.title('Number of reads after pre-processing')
    ax1.barh(x, reads, width)
    ax1.set_xscale('log')
    ax1.set_xlabel('Number of reads [log10]')
    ax1.invert_yaxis()
    ax1.grid(axis='x', c='black', ls="dotted")
    plt.yticks(x, labels, rotation=0)
    plt.show()


def hplotSTARstats_readLen(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    # length of reads
    readsLen = df['                      Average input read length |'].astype(int)

    fig, ax1 = plt.subplots(figsize=(5, len(labels) / 4), dpi=dpi)
    plt.title('Average input read length')
    ax1.barh(x, readsLen, width)
    ax1.set_xlabel('Average input read length [nt]')
    ax1.invert_yaxis()
    ax1.grid(axis='x', c='black', ls="dotted")
    ax1.set_xlim(0, 175)
    plt.yticks(x, labels, rotation=0)
    plt.show()


def hplotSTARstats_mapping(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    # reads mapping
    mapped_uniq = df['                        Uniquely mapped reads % |'].str.strip("%").astype(float)
    mapped_multi = df['             % of reads mapped to multiple loci |'].str.strip("%").astype(float)
    mapped_tooMany = df['             % of reads mapped to too many loci |'].str.strip("%").astype(float)
    unmapped_mismatches = df['       % of reads unmapped: too many mismatches |'].str.strip("%").astype(float)
    unmapped_tooShort = df['                 % of reads unmapped: too short |'].str.strip("%").astype(float)
    unmapped_other = df['                     % of reads unmapped: other |'].str.strip("%").astype(float)

    fig, ax1 = plt.subplots(figsize=(5, len(labels) / 4), dpi=dpi)
    plt.title('STAR mapping statistics')

    p1 = ax1.barh(x, mapped_uniq, width, color='green')
    p2 = ax1.barh(x, mapped_multi, width, left=mapped_uniq, color='limegreen')
    p3 = ax1.barh(x, mapped_tooMany, width, left=mapped_uniq + mapped_multi, color='yellow')
    p4 = ax1.barh(x, unmapped_mismatches, width, left=mapped_uniq + mapped_multi + mapped_tooMany, color='grey')
    p5 = ax1.barh(x, unmapped_tooShort, width, left=mapped_uniq + mapped_multi + mapped_tooMany + unmapped_mismatches,
                  color='orange')
    p6 = ax1.barh(x, unmapped_other, width,
                  left=mapped_uniq + mapped_multi + mapped_tooMany + unmapped_mismatches + unmapped_tooShort,
                  color='red')

    plt.xlabel('Cumulative mapping [%]')
    ax1.invert_yaxis()
    plt.yticks(x, labels, rotation=0)
    plt.xticks(np.arange(0, 101, 10))
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]), (
    'uniquely mapped', 'multimappers', 'too many mapped', 'unmapped (mistmateches)', 'unmapped (too short)',
    'unmapped (other)'))

    plt.show()


def hplotSTARstats_mistmatches(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    # mistmatches and deletions
    mismatches = df['                      Mismatch rate per base, % |'].str.strip("%").astype(float)
    deletions = df['                         Deletion rate per base |'].str.strip("%").astype(float)
    insertions = df['                        Insertion rate per base |'].str.strip("%").astype(float)

    fig, ax1 = plt.subplots(figsize=(5, len(labels) / 4), dpi=dpi)
    plt.title('Mapping errors')
    ax1.barh(x - width, mismatches, width, label='mismatches')
    ax1.barh(x, deletions, width, label='deletions')
    ax1.barh(x + width, insertions, width, label='insertions')
    ax1.grid(axis='x', c='black', ls="dotted")
    ax1.legend()
    plt.yticks(x, labels, rotation=0)
    plt.xlabel('Errors rate per base [%]')
    ax1.invert_yaxis()
    plt.show()


def hplotSTARstats_chimeric(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    reads = df['                          Number of input reads |'].astype(int)
    # Splices and chimeric reads
    splices = df['                       Number of splices: Total |'].astype(int) / reads
    chmieric = df['                            % of chimeric reads |'].str.strip("%").astype(float)

    fig, ax1 = plt.subplots(figsize=(5, len(labels) / 4), dpi=dpi)
    plt.title('Chimeric reads')

    ax1.barh(x - width / 2, splices, width, label='splices')
    ax1.barh(x + width / 2, chmieric, width, label='chmieric')
    ax1.grid(axis='x', c='black', ls="dotted")
    ax1.legend()
    plt.yticks(x, labels, rotation=0)
    plt.xlabel('% of reads')
    ax1.invert_yaxis()
    plt.show()


def hplotSTARstats(df=pd.DataFrame(), dpi=150):
    '''

    :param df: DataFrame
    :param dpi: int, default=150
    :return:
    '''
    hplotSTARstats_reads(df=df, dpi=dpi)
    hplotSTARstats_readLen(df=df, dpi=dpi)
    hplotSTARstats_mapping(df=df, dpi=dpi)
    hplotSTARstats_mistmatches(df=df, dpi=dpi)
    hplotSTARstats_chimeric(df=df, dpi=dpi)
