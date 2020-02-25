import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import TTools.profiles as profiles

#### PCA

def plotPCA(data=pd.DataFrame(), names=list(), title=str(), PClimit=1,figsize = (7,7)):
    nPCA = len([col for col in data.columns if 'PC' in col])
    axes = [ i +1 for i in range(nPCA)[:-1]]

    for nPC in axes:
        fig = plt.figure(figsize)
        #         ax = fig.add_subplot(nPCA-1,1,nPC)
        ax = fig.add_subplot(1 ,1 ,1)
        a = data['PC ' +str(nPC)].tolist()
        b = data['PC ' +str(nPC+1)].tolist()

        ax.scatter(x=a ,y=b ,color='lightgray')

        if names:
            for i, txt in enumerate(data.index.tolist()):
                if txt in names:
                    ax.annotate(txt, (a[i], b[i]))

        #         for mark, c in zip(l1_fake,l1_fakeColor):
        #             markTemp = data.T[mark]
        #             ax.scatter(x=markTemp['PC'+str(nPC)],y=markTemp['PC'+str(nPC+1)],marker='X', color=c)

        ax.legend()
        ax.grid(True)
        #         ax.set_xlim(-1,1)
        plt.xlabel('PC ' +str(nPC))
        plt.ylabel('PC ' +str(nPC +1))
        plt.title(title)
        plt.show()

        if nPC==PClimit: break

### Peaks metaplot

def plotCumulativePeaks(csv_path=str(), df2=pd.DataFrame(), local_pos=list(), dpi=150,
                        title=None, start=None, stop=None, window=20, figsize=(4, 3),
                        color1='green', color2="magenta", lc='red', ylim=None):
    reference = pd.read_csv(csv_path, index_col=0)

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

def plot_as_box_plot(df=pd.DataFrame(),title=None, start=None, stop=None,figsize=(7,3),ylim=(None,0.01), dpi=150, color='green', h_lines=list(), lc='red',offset=0):
    '''Plots figure similar to box plot: median, 2 and 3 quartiles and min-max range

    :param df: Dataframe() containing following columns:```['position'] ['mean'] ['median'] ['std']```
        optionally ```['nucleotide'] ['q1'] ['q3'] ['max'] ['min']```
    :param title: str()
    :param start: int()
    :param stop: int()
    :param figsize: touple(). Default = (7,4)
    :param ylim: touple() OY axes lim. Default = (None,0.01)
    :param dpi: int()
    :param color: str()
    :param h_lines: list() of horizontal lines
    :param lc: str() color of horizontal lines
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

def plot_to_compare(df=pd.DataFrame(), df2=None, color1='black', color2='darkred',
                    label=str(), title=None, start=None, stop=None, figsize=(7,3),
                    ylim=(None,0.01), h_lines=list(), dpi=150,offset=300,
                    csv_path='/home/tturowski/notebooks/RDN37_csv_path_collapsed.csv'):
    reference = pd.read_csv(csv_path, index_col=0) # reading reference
    dataset, s2 = df[start:stop], reference[start:stop] # prepating datasets
    #plotting
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    plt.axhline(0, color='red')
    if len(dataset.columns) == 4: #if only two experiments
        ax1.plot(dataset.index-offset, dataset['mean'], color1, label=label)
        ax1.fill_between(dataset.index-offset, dataset['min'], dataset['max'], color=color1, alpha=0.3, label='range (min-max)')
    else: #if more than two experiments
        ax1.plot(dataset.index-offset, dataset['median'], color1, label=label)
        ax1.fill_between(dataset.index-offset, dataset['q1'], dataset['q3'], color=color1, alpha=0.2, label='range (2nd-3rd quartile)')
#         ax1.fill_between(dataset.index-offset, dataset['min'], dataset['max'], color=color1, alpha=0.07, label='range (min-max)')
    ax1.set_xlabel('position')
    ax1.set_ylim(ylim)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('fraction of reads', color='black')
    for tl in ax1.get_yticklabels():
        tl.set_color('black')

    ax1.plot(s2.index-offset, s2['median'], 'green', label='Rpa190')
    ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], color='green', alpha=0.2, label='range (q2-q3)')
    ax1.fill_between(s2.index-offset, s2['min'], s2['max'], color='green', alpha=0.07, label='range (min-max)')
    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color='red')
    ax1.legend()

def plot_diff(dataset=pd.DataFrame(), ranges='mm', label=str(), start=None, stop=None, plot_medians=True,
              plot_ranges=True, figsize=(7, 3), ylim=(None,0.01), h_lines=list(), offset=0,
              reference='/home/tturowski/notebooks/RDN37_reference_collapsed.csv'):
    '''Plot given dataset and reference, differences are marked

    :param dataset: DataFrame() containing following columns:```['position'] ['mean'] ['median'] ['std']```
        optionally ```['nucleotide'] ['q1'] ['q3'] ['max'] ['min']```
    :param ranges: str() mm : min-max or qq : q1-q3
    :param label: str()
    :param start: int()
    :param stop: int()
    :param plot_medians: plot medians
    :param plot_ranges: plot ranges
    :param figsize: touple(). Default = (7,4)
    :param ylim: touple() OY axes lim. Default = (None,0.01)
    :param h_lines: list() of horizontal lines
    :param reference: str() with path or DataFrame() to reference data
    :return: Plot with marked sequences
    '''

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

def plot_heatmap(df=pd.DataFrame(), title='Heat map of differences between dataset and reference plot for RDN37-1', vmin=None,
                 vmax=None, figsize=(20,10)):
    '''Plot heat map of differences, from dataframe generated by compare1toRef(dataset, heatmap=True) function

    :param df: DataFrame()
    :param title: str()
    :param vmin:
    :param vmax:
    :param figsize: touple()
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
