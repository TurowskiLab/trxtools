from turtle import st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import trxtools.profiles as profiles
from adjustText import adjust_text
import seaborn as sns
import matplotlib
from matplotlib_venn import venn2, venn3, venn3_circles
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os, math

#### general functions ####

def select_colors(n=int(), name_cmap='Spectral'):
    '''_summary_

    :param n: _description_, defaults to int()
    :type n: _type_, optional
    :param name_cmap: _description_, defaults to 'Spectral'
    :type name_cmap: str, optional
    :return: _description_
    :rtype: _type_
    '''
    cmap = matplotlib.cm.get_cmap(name_cmap)
    numbers = np.linspace(0, 1, n).tolist()
    return [cmap(i) for i in numbers]

####### Statistical plots  ########

### Volcano

def enhancedVolcano(df, x_fc='log2FoldChange', y_pval='padj', pval_cutoff=0.01, fc_cutoff=2,
                         plot_n=True, labels=True, xlim=None, ylim=None,
                         figsize=(4, 4), dpi=300, title=None, save=None):
    '''Generate an enhanced volcano plot based on DataFrame values.

    :param df: Input DataFrame containing the data.
    :type df: DataFrame
    :param x_fc: The column name for the x-axis values (fold change). Defaults to 'log2FoldChange'.
    :type x_fc: str, optional
    :param y_pval: The column name for the y-axis values (adjusted p-value). Defaults to 'padj'.
    :type y_pval: str, optional
    :param pval_cutoff: The p-value cutoff for significance. Defaults to 0.01.
    :type pval_cutoff: float, optional
    :param fc_cutoff: The fold change cutoff for significance. Defaults to 2.
    :type fc_cutoff: float, optional
    :param plot_n: Whether to plot non-significant points. Defaults to True.
    :type plot_n: bool, optional
    :param labels: Whether to add labels to significant points. Defaults to True.
    :type labels: bool, optional
    :param xlim: The x-axis limits. Defaults to None.
    :type xlim: tuple, optional
    :param ylim: The y-axis limits. Defaults to None.
    :type ylim: tuple, optional
    :param figsize: The figure size. Defaults to (4, 4).
    :type figsize: tuple, optional
    :param dpi: The resolution of the saved figure. Defaults to 300.
    :type dpi: int, optional
    :param title: The title of the plot. Defaults to None.
    :type title: str, optional
    :param save: The file path to save the figure. Defaults to None.
    :type save: str, optional

    :return: None

    '''

    # Filter DataFrame into different subsets based on conditions
    df_s_fc = df[(df[y_pval] < pval_cutoff) & ((df[x_fc] <= -fc_cutoff) | (df[x_fc] >= fc_cutoff))]  # significant and fold change
    df_s_nf = df[(df[y_pval] < pval_cutoff) & ((df[x_fc] > -fc_cutoff) & (df[x_fc] < fc_cutoff))]   # significant, no fold change
    df_ns_fc = df[(df[y_pval] >= pval_cutoff) & ((df[x_fc] <= -fc_cutoff) | (df[x_fc] >= fc_cutoff))]  # non-significant, fold change
    df_ns_nf = df[(df[y_pval] >= pval_cutoff) & ((df[x_fc] > -fc_cutoff) & (df[x_fc] < fc_cutoff))]  # non-significant, no fold change
    
    # Create a new figure
    plt.figure(figsize=figsize, dpi=dpi)
    
    # Plot the scatter plots for different subsets
    a = sns.scatterplot(data=df_s_fc, x=df_s_fc[x_fc], y=-np.log10(df_s_fc[y_pval]), color='red', s=4, alpha=0.4)
    b = sns.scatterplot(data=df_s_nf, x=df_s_nf[x_fc], y=-np.log10(df_s_nf[y_pval]), color='blue', s=4, alpha=0.4)
    c = sns.scatterplot(data=df_ns_fc, x=df_ns_fc[x_fc], y=-np.log10(df_ns_fc[y_pval]), color='green', s=4, alpha=0.4)
    d = sns.scatterplot(data=df_ns_nf, x=df_ns_nf[x_fc], y=-np.log10(df_ns_nf[y_pval]), color='grey', s=4, alpha=0.4)
    
    # Set the x-axis and y-axis labels
    a.set(ylabel='-log10(Pval)')
    a.set(xlabel='log2FC')
    
    # Add vertical and horizontal lines
    a.axvline(0, alpha=0.5, c='grey')
    a.axhline(-np.log10(pval_cutoff), alpha=0.75, c='black', ls="--")
    a.axvline(-fc_cutoff, alpha=0.75, c='black', ls="--")
    a.axvline(fc_cutoff, alpha=0.75, c='black', ls="--")
    
    # Set limist
    if xlim:
        a.set(xlim=xlim)
    if ylim:
        a.set(ylim=ylim)
    
    if labels:
        texts = []
        # Add labels to significant points
        for x, y, s in zip(df_s_fc[x_fc].tolist(), (-np.log10(df_s_fc[y_pval])).tolist(), df_s_fc.index.tolist()):
            texts.append(plt.text(x, y, s, fontsize=12))
        adjust_text(texts, only_move={'points': 'y', 'texts': 'y'}, arrowprops=dict(arrowstyle="->", color='r', lw=0.5), expand_points=(2, 2))
        
        # Increase font size for titles and labels
        for item in ([a.title, a.xaxis.label, a.yaxis.label] + a.get_xticklabels() + a.get_yticklabels()):
            item.set_fontsize(15)
    
    # Create legend patches for different subsets
    patch_labels = ["p < {} and log2FC (n={})".format(pval_cutoff, len(df_s_fc)),
                    "p < {}".format(pval_cutoff),
                    "log2FC",
                    "n.s."]
    legend_patches = [mpatches.Patch(color='red', label=patch_labels[0]),
                      mpatches.Patch(color='blue', label=patch_labels[1]),
                      mpatches.Patch(color='green', label=patch_labels[2]),
                      mpatches.Patch(color='grey', label=patch_labels[3])]
    
    # Adjust the layout and add legend
    plt.subplots_adjust(top=0.85)
    plt.legend(handles=legend_patches, fontsize=8, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.12), frameon=False)
    
    # Set the title
    plt.suptitle(title)
    
    # Save the figure if a file path is provided
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    # Show the plot
    plt.show()

def vennDiagram(*dataframes, labels=['df1','df2','df3'], colors=('skyblue', 'lightgreen', 'lightpink'),
                      title=None, save=None):
    '''
    '''
    # Prepare the data for the Venn diagram
    sets = [set(df['gene_name']) for df in dataframes]
    labels = [labels[i] for i in range(len(dataframes))]

    # Determine the type of Venn diagram based on the number of DataFrames
    if len(dataframes) == 2:
        venn_func = venn2
    elif len(dataframes) == 3:
        venn_func = venn3
    else:
        raise ValueError("Venn diagram is supported for 2 or 3 DataFrames only.")

    # Create the Venn diagram
    venn = venn_func(sets, set_labels=labels)

    # Customize the Venn diagram colors
    for i, patch in enumerate(venn.patches):
        if patch: patch.set_color(colors[i % len(colors)])

    # Customize the Venn diagram labels
    for i, label in enumerate(venn.set_labels):
        label.set_fontsize(12)
        label.set_fontweight('bold')

    # Customize the Venn diagram title
    plt.title(title, fontsize=14, fontweight='bold')

    # Customize the Venn diagram circles (applicable for 3 DataFrames only)
    if len(dataframes) == 3:
        venn3_circles(sets)
    
    # Save the figure if a file path is provided
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    # Display the Venn diagram
    plt.show()

#### GO term

def GOterm(df, x='fdr',y='term.label',normalizer='pValue', cutoff=0.05, count=20, figsize=(4, 4), dpi=300, title=None, fname=None):
    
    # Select data
    df = df[(df[normalizer] < cutoff ) & (df[x] < cutoff )]
    
    if count == None:
        df = df.sort_values(normalizer,ascending=True)
    else:
        df = df.sort_values(normalizer,ascending=True)[:count]
    
    if len(df) == 0:
        print("No significant GO terms found.")
        return True

    hight = 1.1 + (len(df)/6.9) # 6.9 is a magic number to adjust the height of the plot
    figsize = (figsize[0], hight)

    y_pos = np.arange(len(df))
    # Create a new figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Color map and normalization
    data_color = [x / cutoff for x in df[normalizer]]
    my_cmap = plt.cm.get_cmap('autumn')
    colors = my_cmap(data_color)
    
    # Plot barplot
    ax.barh(y_pos, df[x], align='center', color = colors)
    ax.set_yticks(y_pos, labels=df[y])
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('FDR')
    ax.set_title(title)
    
    # color map
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,cutoff))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label(normalizer, rotation=270,labelpad=25)
    
    # output
    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

#### PCA
def PCA(data=pd.DataFrame(), names=[], title="", PClimit=1,figsize=(7,7), PCval=[]):
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

## boxplots

def boxplot1(data,labels=None,title="",figsize=(7, 6),dpi=150,log=False,lim=None,
            name_cmap='Spectral',vert=1,color=None,grid=False,fname=None):

    fig = plt.figure(figsize=figsize,dpi=dpi)
    ax = fig.add_subplot()

    # Creating axes instance
    bp = ax.boxplot(data, patch_artist = True,
                    notch ='True', vert = vert)

    ### colors ###
    if color:
        colors = [color]*len(data)
    else:
        colors = select_colors(n=len(data),name_cmap=name_cmap)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # whiskers
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 1.5,
                    linestyle =":")
    # caps
    for cap in bp['caps']:
        cap.set(color ='black', linewidth = 2)
    # median
    for median in bp['medians']:
        median.set(color ='black', linewidth = 3)
    
    # style of fliers
    for flier in bp['fliers']:
        flier.set(marker ='.',
                color ='grey',
                alpha = 0.5)

    # x-axis labels
    if labels: pass
    elif (isinstance(data, pd.DataFrame) or isinstance(data, pd.Series)):
        labels = data.index.tolist()
    else:
        labels = ["data_"+str(i) for i in range(1,len(data)+1)]
    
    if vert==0:
        ax.set_yticklabels(labels)
        if lim:
            ax.set_xlim(lim)
        if log==True:
            ax.set_xscale('log')
        if grid==True:
            ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
    elif vert==1:
        ax.set_xticklabels(labels)
        if lim:
            ax.set_ylim(lim)
        if log==True:
            ax.set_yscale('log')
        if grid==True:
            ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

    # Removing top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.title(title)
    
    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

### Peaks metaplot
def cumulativePeaks(ref, df2=pd.DataFrame(), local_pos=list(), dpi=150,
                        title="", start=None, stop=None, window=50, figsize=(4,3),equal_weight=False,
                        color1='green', color2="magenta", lc='red',fname=None,use="mean",ylim=None):
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
        s3_dataset2 = df2['median'][loc - window:loc + window]
        if equal_weight == True: #to be tested
            s2_pseudocounts = s2_dataset1.min() / 100
            s2_dataset1 = s2_dataset1.add(s2_pseudocounts) / s2_dataset1.add(s2_pseudocounts).sum()
            s3_pseudocounts = s3_dataset2.min() / 100
            s3_dataset2 = s3_dataset2.add(s3_pseudocounts) / s3_dataset2.add(s3_pseudocounts).sum()
        df_dataset1[i] = s2_dataset1.reset_index()['median'] #for multigene metaplot this value could/should be normalized
        df_dataset2[i] = s3_dataset2.reset_index()['median'] #for multigene metaplot this value could/should be normalized

    if use=='mean':
        s_data1 = df_dataset1.mean(axis=1)
        s_data2 = df_dataset2.mean(axis=1)
    elif use=='median':
        s_data1 = df_dataset1.median(axis=1)
        s_data2 = df_dataset2.median(axis=1)
    elif use=='sum':
        s_data1 = df_dataset1.sum(axis=1)
        s_data2 = df_dataset2.sum(axis=1)

    # plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel(f'fraction of reads ({use})')
    if ylim:
        ax1.set_ylim(ylim)
    ax1.plot(np.arange(-window, window), s_data1, color=color1)

    # plotting dataset2
    ax1.plot(np.arange(-window, window), s_data2, color=color2)

    ax1.axvline(0, color=lc, alpha=0.5)
    # ax1.legend(loc=2)

    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

### profiles
def plot_as_box_plot(df=pd.DataFrame(),title="", start=None, stop=None,name='median',
                     figsize=(7,3),ylim=(None,0.01), dpi=150, color='green',
                     h_lines=[], lc="red",offset=0,fname=None):
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
    :param filename: default None, string with path to file
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
    ax1.plot(s2.index-offset, s2['median'], color=color,label=name)
    if set(['q1','q3']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], label='range (2nd-3rd quartile)', color=color, alpha=0.2)
    if set(['min','max']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)
    
    if not start: start=min(s2.index-offset)
    if not stop: stop=max(s2.index-offset)
    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color=lc)
    ax1.legend()
    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()


def plotAndFolding(df=pd.DataFrame(),dG=pd.Series(), title="", start=None, stop=None,
                     figsize=(7,3),ylim=(None,0.01), dpi=150, color='green', name='median',
                     h_lines=[], lc="red",offset=0,fname=None):
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
    dG = dG[start:stop]
    #plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    ax1.set_ylim(ylim)
    ax1.plot(s2.index-offset, s2['median'], color=color,label=name)
    if set(['q1','q3']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], label='range (2nd-3rd quartile)', color=color, alpha=0.2)
    if set(['min','max']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)
    
    if not start: start=min(s2.index-offset)
    if not stop: stop=max(s2.index-offset)
    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color=lc)
    ax1.legend(loc=2)
        
    #plotting delta G
    ax2 = ax1.twinx()
    ax2.plot(dG.index-offset, dG, color="orange", label="dG_@30")
    ax2.set_ylabel(r'delta'+' G', color='orange')
#     ax2.set_yticks(np.arange(0,1,0.1))
    ax2.set_ylim(-50,20)
    for tl in ax2.get_yticklabels():
        tl.set_color('orange')
    ax2.grid()
    ax2.legend(loc=1)
    
    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()


def plot_to_compare(ref, df=pd.DataFrame(), color1='green', color2='black',
                    ref_label="", label="", title="", start=None, stop=None, figsize=(7,3),
                    ylim=(None,0.01), h_lines=[], lc="red", dpi=150,offset=300,fname=None):
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
        ax1.fill_between(dataset.index-offset, dataset['q1'], dataset['q3'], color=color2, alpha=0.2, label='range (q2-q3)')
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

    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

def plot_diff(ref, dataset=pd.DataFrame(), ranges='mm', label1="reference", label2="", 
              title="", start=None, stop=None, plot_medians=True, plot_ranges=True,
              dpi=150, figsize=(7, 3), ylim=(None,0.01), h_lines=list(),lc="red", offset=0,fname=None):
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

    differences_df = profiles.compareMoretoRef(dataset=dataset, ranges=ranges, ref=reference)[start:stop]
    dataset, s2 = dataset[start:stop], reference[start:stop]  # prepating datasets
    # plotting
    fig, ax1 = plt.subplots(figsize=figsize,dpi=dpi)
    ax1.fill_between(differences_df.index-offset, differences_df['ear_min'], differences_df['ear_max'], color='red',
                     where=(differences_df['ear_max'] > 0), label='increased occupancy (' + ranges_dict[ranges] + ')')
    ax1.fill_between(differences_df.index-offset, differences_df['rae_min'], differences_df['rae_max'], color='blue',
                     where=(differences_df['rae_max'] > 0), label='decreased occupancy (' + ranges_dict[ranges] + ')')
    if plot_medians == True:
        ax1.plot(dataset.index-offset, dataset['median'], 'black', label=label2)
        ax1.plot(s2.index-offset, s2['median'], 'green', label=label1)
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
    ax1.set_ylabel('fraction of reads', color='black')
    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color=lc)
    ax1.set_title(title)
    plt.legend()
    plt.title(title)

    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

def plot_heatmap(df=pd.DataFrame(), title='Heatmap of differences between dataset and reference plot for RDN37-1',
                 vmin=None, vmax=None, figsize=(20,10),dpi=300,fname=None):
    '''Plot heat map of differences, from dataframe generated by compare1toRef(dataset, heatmap=True) function

    :param df: DataFrame
    :param title: str
    :param vmin:
    :param vmax:
    :param figsize: tuple, default (20,10)
    :return:
    '''

    fig, ax = plt.subplots(figsize=figsize,dpi=dpi)
    if not vmin:
        vmin = -np.absolute(df.max().median())
    if not vmax:
        vmax = df.max().median()
    heatmap = ax.pcolormesh(df.transpose(), cmap='seismic', vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(list(df.columns.values))) + 0.5, minor=False)
    ax.set_yticklabels(list(df.columns.values), minor=False)
    fig.colorbar(heatmap)
    ax.set_title(title)

    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

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

def classes(df,fraction=True,dpi=300, save=None,d_color=None):
    if fraction==True:
        df = df.groupby('type').agg("sum").div(df.sum()).multiply(100)
        df = df.drop(columns='type')
    elif fraction==False:
        df = df.groupby('type').agg("sum")
              
    labels = df.columns.tolist()
    x = np.arange(len(labels))
    width = 0.3
    
    fig, ax1 = plt.subplots(figsize=(len(labels) / 4, 3), dpi=dpi)
    plt.title('Classes of mapped reads')
    
    if d_color is None:
        d_color = {'intergenic_region':"orange",
             'ncRNA'              :"yellow",
             'protein_coding'     :"blue",
             'pseudogene'         :"grey",
             'rRNA'               :"green",
             'snRNA'              :"black",
             'snoRNA'             :"pink",
             'tRNA'               :"darkred"}
        
    s0 = pd.Series(index=df.columns, data=0)

    for i,row in df.iterrows():
        ax1.bar(x, row, width, bottom=s0, label=i, color=d_color[i])
        s0 = s0+row

    plt.ylabel('Cumulative mapping [%]')
    plt.xticks(x, labels, rotation=90)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # Save the figure if a file path is provided
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
    plt.show()