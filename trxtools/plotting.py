from turtle import st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import trxtools.profiles as profiles
import trxtools.metaprofiles as meta
from adjustText import adjust_text
import seaborn as sns
import matplotlib
from matplotlib_venn import venn2, venn3, venn3_circles
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os, math

#### general functions ####

def select_colors(n=int(), name_cmap='Spectral'):
    '''Select a list of colors from a colormap.

    :param n: Number of colors to select, defaults to int()
    :type n: int, optional
    :param name_cmap: Name of the colormap to use, defaults to 'Spectral'
    :type name_cmap: str, optional

    :return: List of colors
    :rtype: list
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
    :type df: pd.DataFrame
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

    :example:

    >>> enhancedVolcano(df, x_fc='log2FoldChange', y_pval='padj', pval_cutoff=0.01, fc_cutoff=2, title='Volcano Plot')
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
    
    # Set limits
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
    '''Generate a Venn diagram for 2 or 3 DataFrames.

    :param dataframes: Input DataFrames containing the data.
    :type dataframes: list of pd.DataFrame
    :param labels: Labels for the DataFrames. Defaults to ['df1','df2','df3'].
    :type labels: list of str, optional
    :param colors: Colors for the Venn diagram. Defaults to ('skyblue', 'lightgreen', 'lightpink').
    :type colors: tuple of str, optional
    :param title: The title of the plot. Defaults to None.
    :type title: str, optional
    :param save: The file path to save the figure. Defaults to None.
    :type save: str, optional

    :return: None

    :example:

    >>> vennDiagram(df1, df2, labels=['Group 1', 'Group 2'], title='Venn Diagram')
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

def GOterm(df, x='fdr', y='term.label', normalizer='pValue', cutoff=0.05, count=20, figsize=(4, 4), dpi=300, title=None, fname=None):
    '''Generate a bar plot for GO terms based on significance.

    :param df: Input DataFrame containing GO term data.
    :type df: pd.DataFrame
    :param x: The column name for the x-axis values (e.g., FDR). Defaults to 'fdr'.
    :type x: str, optional
    :param y: The column name for the y-axis values (e.g., term label). Defaults to 'term.label'.
    :type y: str, optional
    :param normalizer: The column name used for normalization (e.g., p-value). Defaults to 'pValue'.
    :type normalizer: str, optional
    :param cutoff: The cutoff value for significance. Defaults to 0.05.
    :type cutoff: float, optional
    :param count: The number of top terms to display. Defaults to 20.
    :type count: int, optional
    :param figsize: The figure size. Defaults to (4, 4).
    :type figsize: tuple, optional
    :param dpi: The resolution of the saved figure. Defaults to 300.
    :type dpi: int, optional
    :param title: The title of the plot. Defaults to None.
    :type title: str, optional
    :param fname: The file path to save the figure. Defaults to None.
    :type fname: str, optional

    :return: None

    :example:

    >>> GOterm(df, x='fdr', y='term.label', normalizer='pValue', cutoff=0.05, count=20, title='GO Term Plot')
    '''
    # Select data
    df = df[(df[normalizer] < cutoff) & (df[x] < cutoff)]
    
    if count is None:
        df = df.sort_values(normalizer, ascending=True)
    else:
        df = df.sort_values(normalizer, ascending=True)[:count]
    
    if len(df) == 0:
        print("No significant GO terms found.")
        return True

    height = 1.1 + (len(df) / 6.9)  # 6.9 is a magic number to adjust the height of the plot
    figsize = (figsize[0], height)

    y_pos = np.arange(len(df))
    # Create a new figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Color map and normalization
    data_color = [x / cutoff for x in df[normalizer]]
    my_cmap = plt.cm.get_cmap('autumn')
    colors = my_cmap(data_color)
    
    # Plot barplot
    ax.barh(y_pos, df[x], align='center', color=colors)
    ax.set_yticks(y_pos, labels=df[y])
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('FDR')
    ax.set_title(title)
    
    # color map
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0, cutoff))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label(normalizer, rotation=270, labelpad=25)
    
    # output
    if fname:
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()

#### PCA
def PCA(data=pd.DataFrame(), names=[], title="", PClimit=1, figsize=(7, 7), PCval=[]):
    '''Plot PCA plot.

    :param data: DataFrame containing PCA results.
    :type data: pd.DataFrame
    :param names: List of names to annotate.
    :type names: list, optional
    :param title: Title of the plot.
    :type title: str, optional
    :param PClimit: Number of principal components to plot. Defaults to 1.
    :type PClimit: int, optional
    :param figsize: Size of the figure. Defaults to (7, 7).
    :type figsize: tuple, optional
    :param PCval: List of explained variance for each principal component. Defaults to an empty list.
    :type PCval: list, optional

    :return: None

    :example:

    >>> PCA(data, names=['Sample1', 'Sample2'], title='PCA Plot', PClimit=2)
    '''
    nPCA = len([col for col in data.columns if 'PC' in col])
    axes = [i + 1 for i in range(nPCA)[:-1]]
    if not PCval:
        PCval = [""] * (PClimit + 1)

    for nPC in axes:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)

        texts = []

        if 'group' in data.columns.tolist():
            for group, df_temp in data.groupby('group'):
                a = df_temp['PC' + str(nPC)].tolist()
                b = df_temp['PC' + str(nPC + 1)].tolist()
                ax.scatter(x=a, y=b, label=group, cmap="paired")
                for x, y, s in zip(a, b, df_temp.index.tolist()):
                    if s in names:
                        texts.append(plt.text(x, y, s))

        else:
            a = data['PC' + str(nPC)].tolist()
            b = data['PC' + str(nPC + 1)].tolist()
            ax.scatter(x=a, y=b, color='lightgray')
            for x, y, s in zip(a, b, data.index.tolist()):
                if s in names:
                    texts.append(plt.text(x, y, s))

        ax.legend()
        ax.grid(True, ls="dotted")
        plt.xlabel('PC ' + str(nPC) + " (" + str(PCval[nPC - 1]) + "%)")
        plt.ylabel('PC ' + str(nPC + 1) + " (" + str(PCval[nPC]) + "%)")
        plt.title(title)
        adjust_text(texts, only_move={'points': 'y', 'texts': 'y'}, arrowprops=dict(arrowstyle="->", color='black', lw=0.5))
        plt.show()

        if nPC == PClimit:
            break

def clusterClusterMap(df):
    '''Generate a clustermap for clusters.

    :param df: DataFrame containing cluster data.
    :type df: pd.DataFrame
    
    :return: None
    
    :example:
    
    >>> clusterClusterMap(df)

    '''
    df_mean = pd.DataFrame()

    n = df['cluster'].max() + 1

    for c in range(0, n):
        a = df[df['cluster'] == c]
        df_mean['cluster ' + str(c) + " (" + str(len(a)) + ")"] = a.mean()

    df_mean = df_mean.T.drop('cluster', axis=1)

    sns.clustermap(df_mean, center=0, cmap="vlag", linewidths=.75, figsize=(7, 7))
    plt.show()

## boxplots

def boxplot1(data, labels=None, title="", figsize=(7, 6), dpi=150, log=False, lim=None,
             name_cmap='Spectral', vert=1, color=None, grid=False, fname=None):
    '''Generate a box plot.

    :param data: Data to plot.
    :type data: list or pd.DataFrame or pd.Series
    :param labels: Labels for the data. Defaults to None.
    :type labels: list, optional
    :param title: Title of the plot. Defaults to an empty string.
    :type title: str, optional
    :param figsize: Size of the figure. Defaults to (7, 6).
    :type figsize: tuple, optional
    :param dpi: Resolution of the figure. Defaults to 150.
    :type dpi: int, optional
    :param log: Whether to use a logarithmic scale. Defaults to False.
    :type log: bool, optional
    :param lim: Limits for the y-axis. Defaults to None.
    :type lim: tuple, optional
    :param name_cmap: Name of the colormap to use. Defaults to 'Spectral'.
    :type name_cmap: str, optional
    :param vert: Whether to plot the boxes vertically. Defaults to 1.
    :type vert: int, optional
    :param color: Color of the boxes. Defaults to None.
    :type color: str, optional
    :param grid: Whether to display a grid. Defaults to False.
    :type grid: bool, optional
    :param fname: File path to save the figure. Defaults to None.
    :type fname: str, optional

    :return: None

    :example:

    >>> boxplot1(data, labels=['A', 'B', 'C'], title='Box Plot', log=True)
    '''
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot()

    # Creating axes instance
    bp = ax.boxplot(data, patch_artist=True, notch='True', vert=vert)

    # Colors
    if color:
        colors = [color] * len(data)
    else:
        colors = select_colors(n=len(data), name_cmap=name_cmap)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # Whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='black', linewidth=1.5, linestyle=":")
    # Caps
    for cap in bp['caps']:
        cap.set(color='black', linewidth=2)
    # Median
    for median in bp['medians']:
        median.set(color='black', linewidth=3)
    
    # Style of fliers
    for flier in bp['fliers']:
        flier.set(marker='.', color='grey', alpha=0.5)

    # X-axis labels
    if labels:
        pass
    elif isinstance(data, (pd.DataFrame, pd.Series)):
        labels = data.index.tolist()
    else:
        labels = ["data_" + str(i) for i in range(1, len(data) + 1)]
    
    if vert == 0:
        ax.set_yticklabels(labels)
        if lim:
            ax.set_xlim(lim)
        if log:
            ax.set_xscale('log')
        if grid:
            ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    elif vert == 1:
        ax.set_xticklabels(labels)
        if lim:
            ax.set_ylim(lim)
        if log:
            ax.set_yscale('log')
        if grid:
            ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    # Removing top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.title(title)
    
    if fname:
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()

### Peaks metaplot
def cumulativePeaks(ref, df2=pd.DataFrame(), local_pos=list(), dpi=150,
                        title="", start=None, stop=None, window=50, figsize=(4,3),equal_weight=False,
                        color1='green', color2="magenta", lc='red',fname=None,use="mean",ylim=None):
    '''Plot cumulative peaks metaplot for a single gene.

    :param ref: Path to CSV file or DataFrame containing reference data.
    :type ref: str or pd.DataFrame
    :param df2: DataFrame containing the second dataset.
    :type df2: pd.DataFrame
    :param local_pos: List of positions of features (peaks/troughs).
    :type local_pos: list
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional
    :param title: Title of the plot. Defaults to an empty string.
    :type title: str, optional
    :param start: Start position for filtering local positions. Defaults to None.
    :type start: int, optional
    :param stop: Stop position for filtering local positions. Defaults to None.
    :type stop: int, optional
    :param window: Window size around each feature. Defaults to 50.
    :type window: int, optional
    :param figsize: Size of the figure. Defaults to (4, 3).
    :type figsize: tuple, optional
    :param equal_weight: Whether to normalize the data to equal weight. Defaults to False.
    :type equal_weight: bool, optional
    :param color1: Color for the reference dataset plot. Defaults to 'green'.
    :type color1: str, optional
    :param color2: Color for the second dataset plot. Defaults to 'magenta'.
    :type color2: str, optional
    :param lc: Color for the vertical line at position 0. Defaults to 'red'.
    :type lc: str, optional
    :param fname: File path to save the figure. Defaults to None.
    :type fname: str, optional
    :param use: Method to aggregate data ('mean', 'median', 'sum'). Defaults to 'mean'.
    :type use: str, optional
    :param ylim: Limits for the y-axis. Defaults to None.
    :type ylim: tuple, optional

    :return: None

    :example:

    >>> cumulativePeaks('reference.csv', df2, local_pos=[100, 200, 300], title='Cumulative Peaks')
    '''
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    # extracting data for metaplot
    df_dataset1 = pd.DataFrame()
    df_dataset2 = pd.DataFrame()

    # filter local max according to start and stop
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

def cumulativeDifference(ref, df2=pd.DataFrame(), local_pos=list(), dpi=150, fill_between=True,
                        title="", start=None, stop=None, window=50, figsize=(4,3),equal_weight=False,
                        color1='green', lc='red',fname=None,use="mean",ylim=None):
    '''Plot cumulative differences for peaks metaplot.

    :param ref: Path to CSV file or DataFrame containing reference data.
    :type ref: str or pd.DataFrame
    :param df2: DataFrame containing the second dataset.
    :type df2: pd.DataFrame
    :param local_pos: List of positions of features (peaks/troughs).
    :type local_pos: list
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional
    :param fill_between: Whether to fill the area between the plot and the x-axis. Defaults to True.
    :type fill_between: bool, optional
    :param title: Title of the plot. Defaults to an empty string.
    :type title: str, optional
    :param start: Start position for filtering local positions. Defaults to None.
    :type start: int, optional
    :param stop: Stop position for filtering local positions. Defaults to None.
    :type stop: int, optional
    :param window: Window size around each feature. Defaults to 50.
    :type window: int, optional
    :param figsize: Size of the figure. Defaults to (4, 3).
    :type figsize: tuple, optional
    :param equal_weight: Whether to normalize the data to equal weight. Defaults to False.
    :type equal_weight: bool, optional
    :param color1: Color for the difference plot. Defaults to 'green'.
    :type color1: str, optional
    :param lc: Color for the vertical line at position 0. Defaults to 'red'.
    :type lc: str, optional
    :param fname: File path to save the figure. Defaults to None.
    :type fname: str, optional
    :param use: Method to aggregate data ('mean', 'median', 'sum'). Defaults to 'mean'.
    :type use: str, optional
    :param ylim: Limits for the y-axis. Defaults to None.
    :type ylim: tuple, optional

    :return: None

    :example:

    >>> cumulativeDifference('reference.csv', df2, local_pos=[100, 200, 300], title='Cumulative Differences')
    '''
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    # extracting data for metaplot
    df_diff = pd.DataFrame()

    # filter local max according to start and stop
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
        df_diff[i] = s3_dataset2.reset_index()['median'] - s2_dataset1.reset_index()['median']
        
    if use=='mean':
        s_data1 = df_diff.mean(axis=1)
    elif use=='median':
        s_data1 = df_diff.median(axis=1)
    elif use=='sum':
        s_data1 = df_diff.sum(axis=1)

    # plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel(f'cumulative difference ({use})')
    if ylim:
        ax1.set_ylim(ylim)

    if fill_between==True:
        s_zeros = pd.Series([0]*len(s_data1))
        s_positive = s_zeros+s_data1[s_data1>0]
        s_negative = s_zeros+s_data1[s_data1<0]

        ax1.fill_between(np.arange(-window, window), s_zeros, s_positive, color='red', alpha=0.2)
        ax1.fill_between(np.arange(-window, window), s_zeros, s_negative, color='blue', alpha=0.2)

    ax1.plot(np.arange(-window, window), s_data1, color=color1)

    ax1.axvline(0, color=lc, alpha=0.5)
    ax1.axhline(0, color='black', alpha=0.5)
    # ax1.legend(loc=2)

    if fname:
        plt.savefig(fname=fname,dpi=dpi,format='png',bbox_inches='tight')
    else:
        plt.show()

### metaprofile for multiple genes
def generateSubplots(dataframes, figsize=(5, 3), dpi=300, save=None):
    '''Generate subplots for multiple dataframes.

    :param dataframes: A list of dataframes to plot. Each element is a list of lists containing title, dataframe, and optional colormap.
    :type dataframes: list
    :param figsize: The size of the figure. Defaults to (5, 3).
    :type figsize: tuple, optional
    :param dpi: The resolution of the figure in dots per inch. Defaults to 300.
    :type dpi: int, optional
    :param save: The file path to save the figure. If not provided, the figure will be displayed.
    :type save: str, optional
    
    :return: None

    :example:

    >>> df1 = pd.DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})
    >>> df2 = pd.DataFrame({'x': [7, 8, 9], 'y': [10, 11, 12]})
    >>> generateSubplots([
    >>>     [['Title 1', df1, None], 
    >>>      ['Title 2', df2, None]]
    >>> ], figsize=(10, 6), dpi=150, save='plot.png')
    '''
    num_rows = len(dataframes)
    num_cols = max(len(df) for df in dataframes)
        
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(figsize[0]*num_cols, figsize[1]*num_rows))

    cmap = plt.cm.tab20

    for i, df_list in enumerate(dataframes):
        for j, p in enumerate(df_list):
            if num_rows == 1 or num_cols == 1:
                ax = axs[j]
            else:
                ax = axs[i, j]
            title = p[0]
            df = p[1]
            if p[2]: cmap = p[2]

            cs = [cmap(i) for i in range(0, len(df.columns))]
            for c, col in enumerate(df.columns):
                ax.plot(df[col], color=cs[c])
            ax.set_title(title)
            ax.set_xlabel('position')
    ax.legend(df.columns, loc="upper right")

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=dpi)
    else:
        plt.show()

def metaprofileAndHeatmap(data_metaplot, data_heatmap, subsets={"title": None}, 
                          agg_type='mean', normalize_internal=True, differential_heatmap=True,
                          figsize=(5, 3), dpi=300, save=None, bins=[50, 10, 50]):
    '''Generate metaprofile and heatmap plots for multiple subsets.

    :param data_metaplot: Data for metaprofile plot.
    :type data_metaplot: dict
    :param data_heatmap: Data for heatmap plot.
    :type data_heatmap: dict
    :param subsets: Dictionary of subsets with titles as keys and subset data as values. Defaults to {"title": None}.
    :type subsets: dict, optional
    :param agg_type: Aggregation type for metaprofile ('mean', 'median', 'sum'). Defaults to 'mean'.
    :type agg_type: str, optional
    :param normalize_internal: Whether to normalize data internally. Defaults to True.
    :type normalize_internal: bool, optional
    :param differential_heatmap: Whether to plot differential heatmap. Defaults to True.
    :type differential_heatmap: bool, optional
    :param figsize: The size of the figure. Defaults to (5, 3).
    :type figsize: tuple, optional
    :param dpi: The resolution of the figure in dots per inch. Defaults to 300.
    :type dpi: int, optional
    :param save: The file path to save the figure. If not provided, the figure will be displayed.
    :type save: str, optional
    :param bins: List of bin sizes for gene scheme. Defaults to [50, 10, 50].
    :type bins: list, optional

    :return: None

    :example:

    >>> metaprofileAndHeatmap(data_metaplot, data_heatmap, subsets={"Subset 1": subset1, "Subset 2": subset2}, 
    >>>                       agg_type='mean', normalize_internal=True, differential_heatmap=True,
    >>>                       figsize=(10, 6), dpi=150, save='plot.png', bins=[50, 10, 50])
    '''
    num_cols = len(subsets)
    heatmap_size = max(len(df) for df in data_heatmap)

    if differential_heatmap == False:
        fig, axs = plt.subplots(3, num_cols, figsize=(figsize[0]*num_cols, figsize[1]*heatmap_size/8),
                                gridspec_kw={'height_ratios': [1, 0.2, heatmap_size/8]}, sharex=True)
    elif differential_heatmap == True:
        fig, axs = plt.subplots(4, num_cols, figsize=(figsize[0]*num_cols, figsize[1]*heatmap_size/4),
                                gridspec_kw={'height_ratios': [1, 0.2, heatmap_size/8, heatmap_size/8]}, sharex=True)
        
    fig.subplots_adjust(hspace=0, wspace=0)
    
    cmap = plt.cm.tab20c
    
    heatmap_max = 0.05
    heatmap_reference = None
    
    for s, title in enumerate(subsets.keys()):
        # prepare data
        df_metaplot = meta.metaprofile(data_metaplot, agg_type=agg_type,
                                        normalize_internal=normalize_internal,
                                        subset=subsets[title]).fillna(0.0)
        df_heatmap = meta.metaprofile(data_heatmap, agg_type=agg_type,
                                        normalize_internal=normalize_internal,
                                        subset=subsets[title]).fillna(0.0)
        
        if s == 0:
            heatmap_max = df_heatmap.max().max()
        
        # preparing length of datasets
        if subsets[title] is not None:
            meta_len = max([len(df) for df in meta.selectSubsetMM(data_metaplot, subset=subsets[title]).values()])
            heat_len = max([len(df) for df in meta.selectSubsetMM(data_heatmap, subset=subsets[title]).values()])
        else:
            meta_len = max([len(df) for df in data_metaplot.values()])
            heat_len = max([len(df) for df in data_heatmap.values()])
        
        if meta_len != heat_len:
            print("Lengths of metaplot and heatmap datasets are different. n is valid for heatmap.")
        
        # preparing indexes
        if not df_metaplot.index.all() == df_heatmap.index.all():
            print("Indexes for metaplot and heatmap are different. Plotting df_heatmap index.")
        index = list(df_heatmap.index.values)
        df_metaplot = df_metaplot.reset_index().drop(columns='index')
        df_heatmap = df_heatmap.reset_index().drop(columns='index')
        
        if s == 0:
            heatmap_reference = df_heatmap

        # axs for plotting
        if num_cols == 1:
            ax1 = axs[0]
            ax_gene = axs[1]
            ax2 = axs[2]
        else:
            ax1 = axs[0, s]
            ax_gene = axs[1, s]
            ax2 = axs[2, s]
            ax3 = axs[3, s]
        
        # Plot using df_metaplot       
        cs = [cmap(i) for i in range(0, len(df_metaplot.columns))]
        for c, col in enumerate(df_metaplot.columns):
            ax1.plot(df_metaplot[col], color=cs[c])
               
        ax1.set_title(title + " (n=" + str(heat_len) + ")")
        
        # plot gene scheme
        s1_gene = pd.DataFrame(index=df_metaplot.index, data={'gene': 0})
        s1_gene[bins[0]:bins[0]+bins[1]] = 1
        h_gene = ax_gene.imshow(s1_gene.values.T, cmap='Blues', aspect='auto', vmin=0, vmax=1)
        
        # Create heatmap using df_heatmap
        df_heatmap_norm = df_heatmap.div(heatmap_max)
        heatmap = ax2.imshow(df_heatmap_norm.values.T, cmap='binary', aspect='auto', vmin=0, vmax=1)
        ax2.set_xlabel('position')
        ax2.set_xticks(np.arange(0, len(index), 10), minor=False)
        ax2.set_xticklabels(index[::10], minor=False)
        if s == 0:
            ax2.set_yticks(np.arange(len(list(df_heatmap_norm.columns.values))), minor=False)
            ax2.set_yticklabels(list(df_heatmap_norm.columns.values), minor=False)
            ax1.set_ylabel('fraction of reads')
            ax_gene.set_yticks([0], minor=False)
            ax_gene.set_yticklabels(["gene position"], minor=False)
        
        else:
            ax1.set_yticks([], minor=False)
            ax_gene.set_yticks([], minor=False)
            ax2.set_yticks([], minor=False)
            
        
        if differential_heatmap == True:
            if s == 0:
                ax3.axis('off')
            else:
                # Create differential heatmap
                df_heatmap_reference_norm = heatmap_reference.div(heatmap_max)
                df_heatmap_diff = df_heatmap_norm.add(0.1).div(df_heatmap_reference_norm.add(0.1)).applymap(np.log2)                
                heatmap2 = ax3.imshow(df_heatmap_diff.values.T, cmap='seismic', aspect='auto', vmin=-4, vmax=4)
                if s == 1:
                    ax3.set_yticks(np.arange(len(list(df_heatmap_norm.columns.values))), minor=False)
                    ax3.set_yticklabels(list(df_heatmap_norm.columns.values), minor=False)
                else:
                    ax3.set_yticks([], minor=False)
    
    # legends
    ax1.legend(df_metaplot.columns, loc='center left', bbox_to_anchor=(1, 0.5))
    ## custom ax for colorbar
    # greyscale
    fig.subplots_adjust(right=0.8)
    cbar_ax2 = fig.add_axes([0.85, 0.4, 0.02, 0.2])
    fig.colorbar(heatmap, cax=cbar_ax2)
    # seismic
    cbar_ax3 = fig.add_axes([0.85, 0.15, 0.02, 0.2])
    fig.colorbar(heatmap2, cax=cbar_ax3)

    if save:
        plt.savefig(save, dpi=dpi)
    else:
        plt.show()

def plot_as_box_plot(df=pd.DataFrame(), title="", start=None, stop=None, name='median',
                     figsize=(7, 3), ylim=(None, 0.01), dpi=150, color='green',
                     h_lines=[], lc="red", offset=0, fname=None):
    '''Plot a figure similar to a box plot: median, 2nd and 3rd quartiles, and min-max range.

    :param df: DataFrame containing the data with columns: ['position', 'mean', 'median', 'std'] and optionally ['nucleotide', 'q1', 'q3', 'max', 'min'].
    :type df: pd.DataFrame
    :param title: Title of the plot.
    :type title: str, optional
    :param start: Start position for the plot.
    :type start: int, optional
    :param stop: Stop position for the plot.
    :type stop: int, optional
    :param name: Name for the median line in the legend. Defaults to 'median'.
    :type name: str, optional
    :param figsize: Size of the figure. Defaults to (7, 3).
    :type figsize: tuple, optional
    :param ylim: Limits for the y-axis. Defaults to (None, 0.01).
    :type ylim: tuple, optional
    :param dpi: Resolution of the figure in dots per inch. Defaults to 150.
    :type dpi: int, optional
    :param color: Color for the median line. Defaults to 'green'.
    :type color: str, optional
    :param h_lines: List of horizontal lines to add to the plot.
    :type h_lines: list, optional
    :param lc: Color for the horizontal lines. Defaults to 'red'.
    :type lc: str, optional
    :param offset: Number to offset position if 5' flank was used. Defaults to 0.
    :type offset: int, optional
    :param fname: File path to save the figure. If not provided, the figure will be displayed.
    :type fname: str, optional

    :return: None

    :example:
    
    >>> df = pd.DataFrame({'position': range(100), 'median': np.random.rand(100), 'q1': np.random.rand(100), 'q3': np.random.rand(100), 'min': np.random.rand(100), 'max': np.random.rand(100)})
    >>> plot_as_box_plot(df, title="Box Plot", start=10, stop=90, name='median', figsize=(10, 5), ylim=(0, 1), dpi=200, color='blue', h_lines=[20, 40, 60], lc='black', offset=5, fname='box_plot.png')
    '''

    if 'nucleotide' in df.columns.values:
        df = df.drop('nucleotide', 1)
    s2 = df[start:stop]
    # plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    ax1.set_ylim(ylim)
    ax1.plot(s2.index-offset, s2['median'], color=color, label=name)
    if set(['q1', 'q3']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['q1'], s2['q3'], label='range (2nd-3rd quartile)', color=color, alpha=0.2)
    if set(['min', 'max']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index-offset, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)
    
    if not start: start = min(s2.index-offset)
    if not stop: stop = max(s2.index-offset)
    for i in [i for i in h_lines if i in range(start-offset, stop-offset)]: ax1.axvline(i, color=lc)
    ax1.legend()
    if fname:
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()

def plotAndFolding(df=pd.DataFrame(), dG=pd.Series(dtype="float"), title="", start=None, stop=None, legend=True,
                   figsize=(7, 3), ylim=(None, 0.01), dpi=150, color='green', name='median',
                   h_lines=[], lc="red", offset=0, fname=None):
    '''Plot a figure similar to a box plot with additional delta G values.

    :param df: DataFrame containing the data with columns: ['position', 'mean', 'median', 'std'] and optionally ['nucleotide', 'q1', 'q3', 'max', 'min'].
    :type df: pd.DataFrame
    :param dG: Series containing delta G values.
    :type dG: pd.Series
    :param title: Title of the plot.
    :type title: str, optional
    :param start: Start position for the plot.
    :type start: int, optional
    :param stop: Stop position for the plot.
    :type stop: int, optional
    :param legend: Whether to display the legend. Defaults to True.
    :type legend: bool, optional
    :param figsize: Size of the figure. Defaults to (7, 3).
    :type figsize: tuple, optional
    :param ylim: Limits for the y-axis. Defaults to (None, 0.01).
    :type ylim: tuple, optional
    :param dpi: Resolution of the figure in dots per inch. Defaults to 150.
    :type dpi: int, optional
    :param color: Color for the median line. Defaults to 'green'.
    :type color: str, optional
    :param name: Name for the median line in the legend. Defaults to 'median'.
    :type name: str, optional
    :param h_lines: List of horizontal lines to add to the plot.
    :type h_lines: list, optional
    :param lc: Color for the horizontal lines. Defaults to 'red'.
    :type lc: str, optional
    :param offset: Number to offset position if 5' flank was used. Defaults to 0.
    :type offset: int, optional
    :param fname: File path to save the figure. If not provided, the figure will be displayed.
    :type fname: str, optional

    :return: None

    :example:

    >>> plotAndFolding(df, dG, title="Plot with Delta G", start=10, stop=90, color='blue', name='median', h_lines=[20, 40, 60], lc='black', offset=5, fname='plot.png')
    '''
    if 'nucleotide' in df.columns.values:
        df = df.drop('nucleotide', 1)
    s2 = df[start:stop]
    dG = dG[start:stop]
    # plotting reference dataset
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads')
    ax1.set_ylim(ylim)
    ax1.plot(s2.index - offset, s2['median'], color=color, label=name)
    if set(['q1', 'q3']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index - offset, s2['q1'], s2['q3'], label='range (2nd-3rd quartile)', color=color, alpha=0.2)
    if set(['min', 'max']).issubset(list(s2.columns.values)):
        ax1.fill_between(s2.index - offset, s2['min'], s2['max'], label='range (min-max)', color=color, alpha=0.07)

    if not start: start = min(s2.index - offset)
    if not stop: stop = max(s2.index - offset)
    for i in [i for i in h_lines if i in range(start - offset, stop - offset)]: ax1.axvline(i, color=lc)
    if legend:
        ax1.legend(loc=2)

    # plotting delta G
    ax2 = ax1.twinx()
    ax2.plot(dG.index - offset, dG, color="orange", label="dG_@30")
    ax2.set_ylabel(r'delta' + ' G', color='orange')
    ax2.set_ylim(-50, 20)
    for tl in ax2.get_yticklabels():
        tl.set_color('orange')
    ax2.grid()

    if legend:
        ax2.legend(loc=1)

    if fname:
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()


def plot_to_compare(ref, df=pd.DataFrame(), color1='green', color2='black', legend=True,
                    ref_label="", label="", title="", start=None, stop=None, figsize=(7, 3),
                    ylim=(None, 0.01), h_lines=[], lc="red", dpi=150, offset=300, fname=None):
    '''Plot a figure to compare two datasets similar to a box plot.

    :param ref: Path to CSV file or DataFrame containing reference data.
    :type ref: str or pd.DataFrame
    :param df: DataFrame containing the data to compare.
    :type df: pd.DataFrame
    :param color1: Color for the reference dataset plot. Defaults to 'green'.
    :type color1: str, optional
    :param color2: Color for the comparison dataset plot. Defaults to 'black'.
    :type color2: str, optional
    :param ref_label: Label for the reference dataset in the legend.
    :type ref_label: str, optional
    :param label: Label for the comparison dataset in the legend.
    :type label: str, optional
    :param title: Title of the plot.
    :type title: str, optional
    :param start: Start position for the plot.
    :type start: int, optional
    :param stop: Stop position for the plot.
    :type stop: int, optional
    :param legend: Whether to display the legend. Defaults to True.
    :type legend: bool, optional
    :param figsize: Size of the figure. Defaults to (7, 3).
    :type figsize: tuple, optional
    :param ylim: Limits for the y-axis. Defaults to (None, 0.01).
    :type ylim: tuple, optional
    :param h_lines: List of horizontal lines to add to the plot.
    :type h_lines: list, optional
    :param lc: Color for the horizontal lines. Defaults to 'red'.
    :type lc: str, optional
    :param dpi: Resolution of the figure in dots per inch. Defaults to 150.
    :type dpi: int, optional
    :param offset: Number to offset position if 5' flank was used. Defaults to 300.
    :type offset: int, optional
    :param fname: File path to save the figure. If not provided, the figure will be displayed.
    :type fname: str, optional

    :return: None

    :example:

    >>> plot_to_compare('reference.csv', df, color1='blue', color2='red', ref_label='Reference', label='Comparison', title='Comparison Plot', start=10, stop=90, h_lines=[20, 40, 60], lc='black', offset=5, fname='comparison_plot.png')
    '''
    # handling reference plot
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    dataset, s2 = df[start:stop], reference[start:stop]  # preparing datasets
    # plotting
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    plt.title(title)
    plt.axhline(0, color='red')

    # plot to compare
    if len(dataset.columns) == 4:  # if only two experiments
        ax1.plot(dataset.index - offset, dataset['mean'], color2, label=label)
        ax1.fill_between(dataset.index - offset, dataset['min'], dataset['max'], color=color2, alpha=0.3, label='range (min-max)')
    else:  # if more than two experiments
        ax1.plot(dataset.index - offset, dataset['median'], color2, label=label)
        ax1.fill_between(dataset.index - offset, dataset['q1'], dataset['q3'], color=color2, alpha=0.2, label='range (q2-q3)')
    ax1.set_xlabel('position')
    ax1.set_ylim(ylim)
    ax1.set_ylabel('fraction of reads', color='black')
    for tl in ax1.get_yticklabels():
        tl.set_color('black')

    # reference plot
    if len(s2.columns) == 4:  # if only two experiments
        ax1.plot(s2.index - offset, s2['mean'], color1, label=ref_label)
        ax1.fill_between(s2.index - offset, s2['min'], s2['max'], color=color1, alpha=0.07, label='range (min-max)')
    else:
        ax1.plot(s2.index - offset, s2['median'], color1, label=ref_label)
        ax1.fill_between(s2.index - offset, s2['q1'], s2['q3'], color=color1, alpha=0.2, label='range (q2-q3)')

    for i in [i for i in h_lines if i in range(start - offset, stop - offset)]: ax1.axvline(i, color=lc)
    if legend:
        ax1.legend()

    if fname:
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()


def plot_diff(ref, dataset=pd.DataFrame(), ranges='mm', label1="reference", label2="", 
              title="", start=None, stop=None, plot_medians=True, plot_ranges=True, legend=True,
              dpi=150, figsize=(7, 3), ylim=(None, 0.01), h_lines=list(), lc="red", offset=0, fname=None):
    '''Plot given dataset and reference, highlighting differences.

    :param ref: Path to CSV file or DataFrame containing reference data.
    :type ref: str or pd.DataFrame
    :param dataset: DataFrame containing the data to compare with columns: ['position', 'mean', 'median', 'std'] and optionally ['nucleotide', 'q1', 'q3', 'max', 'min'].
    :type dataset: pd.DataFrame
    :param ranges: Type of range to highlight ('mm' for min-max or 'qq' for q1-q3). Defaults to 'mm'.
    :type ranges: str, optional
    :param label1: Label for the reference dataset in the legend. Defaults to "reference".
    :type label1: str, optional
    :param label2: Label for the comparison dataset in the legend.
    :type label2: str, optional
    :param title: Title of the plot.
    :type title: str, optional
    :param start: Start position for the plot.
    :type start: int, optional
    :param stop: Stop position for the plot.
    :type stop: int, optional
    :param plot_medians: Whether to plot medians. Defaults to True.
    :type plot_medians: bool, optional
    :param plot_ranges: Whether to plot ranges. Defaults to True.
    :type plot_ranges: bool, optional
    :param legend: Whether to display the legend. Defaults to True.
    :type legend: bool, optional
    :param dpi: Resolution of the figure in dots per inch. Defaults to 150.
    :type dpi: int, optional
    :param figsize: Size of the figure. Defaults to (7, 3).
    :type figsize: tuple, optional
    :param ylim: Limits for the y-axis. Defaults to (None, 0.01).
    :type ylim: tuple, optional
    :param h_lines: List of horizontal lines to add to the plot.
    :type h_lines: list, optional
    :param lc: Color for the horizontal lines. Defaults to 'red'.
    :type lc: str, optional
    :param offset: Number to offset position if 5' flank was used. Defaults to 0.
    :type offset: int, optional
    :param fname: File path to save the figure. If not provided, the figure will be displayed.
    :type fname: str, optional

    :return: None

    :example:

    >>> plot_diff('reference.csv', dataset, ranges='qq', label1='Reference', label2='Comparison', title='Difference Plot', start=10, stop=90, plot_medians=True, plot_ranges=True, h_lines=[20, 40, 60], lc='black', offset=5, fname='difference_plot.png')
    '''
    # handling reference plot
    if isinstance(ref, str):
        reference = pd.read_csv(ref, index_col=0)
    elif isinstance(ref, pd.DataFrame):
        reference = ref

    ranges_dict = {'mm': 'min-max', 'qq': 'q1-q3'}

    differences_df = profiles.compareMoretoRef(dataset=dataset, ranges=ranges, ref=reference)[start:stop]
    dataset, s2 = dataset[start:stop], reference[start:stop]  # preparing datasets
    # plotting
    fig, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
    for i in [i for i in h_lines if i in range(start - offset, stop - offset)]: ax1.axvline(i, color=lc)

    ax1.fill_between(differences_df.index - offset, differences_df['ear_min'], differences_df['ear_max'], color='red',
                     where=(differences_df['ear_max'] > 0), label='increased occupancy (' + ranges_dict[ranges] + ')')
    ax1.fill_between(differences_df.index - offset, differences_df['rae_min'], differences_df['rae_max'], color='blue',
                     where=(differences_df['rae_max'] > 0), label='decreased occupancy (' + ranges_dict[ranges] + ')')
    if plot_medians:
        ax1.plot(dataset.index - offset, dataset['median'], 'black', label=label2)
        ax1.plot(s2.index - offset, s2['median'], 'green', label=label1)
    if plot_ranges:
        if len(dataset.columns) == 4:  # if only two experiments
            ax1.fill_between(dataset.index - offset, dataset['min'], dataset['max'], color='black', alpha=0.3, label='min-max')
        else:  # if more than two experiments
            ax1.fill_between(dataset.index - offset, dataset['q1'], dataset['q3'], color='black', alpha=0.2, label='q1-q3')
            ax1.fill_between(dataset.index - offset, dataset['min'], dataset['max'], color='black', alpha=0.07, label='min-max')
        ax1.fill_between(s2.index - offset, s2['q1'], s2['q3'], color='green', alpha=0.2, label='q1-q3')
        ax1.fill_between(s2.index - offset, s2['min'], s2['max'], color='green', alpha=0.07, label='min-max')
    ax1.set_ylim(ylim)
    ax1.set_xlabel('position')
    ax1.set_ylabel('fraction of reads', color='black')

    ax1.set_title(title)
    if legend:
        plt.legend()
    plt.title(title)

    if fname:
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()


def plot_heatmap(df=pd.DataFrame(), title='Heatmap of differences between dataset and reference plot for RDN37-1',
                 vmin=None, vmax=None, figsize=(20, 10), dpi=300, fname=None):
    '''Plot heatmap of differences between dataset and reference.

    :param df: DataFrame containing the differences.
    :type df: pd.DataFrame
    :param title: Title of the plot. Defaults to 'Heatmap of differences between dataset and reference plot for RDN37-1'.
    :type title: str, optional
    :param vmin: Minimum value for the colormap. Defaults to None.
    :type vmin: float, optional
    :param vmax: Maximum value for the colormap. Defaults to None.
    :type vmax: float, optional
    :param figsize: Size of the figure. Defaults to (20, 10).
    :type figsize: tuple, optional
    :param dpi: Resolution of the figure in dots per inch. Defaults to 300.
    :type dpi: int, optional
    :param fname: File path to save the figure. If not provided, the figure will be displayed.
    :type fname: str, optional

    :return: None

    :example:

    >>> plot_heatmap(df, title='Difference Heatmap', vmin=-1, vmax=1, figsize=(15, 8), dpi=200, fname='heatmap.png')
    '''
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
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
        plt.savefig(fname=fname, dpi=dpi, format='png', bbox_inches='tight')
    else:
        plt.show()

################################################
#############       STAR mapping statistics

def plotSTARstats_reads(df=pd.DataFrame(), dpi=150):
    '''Plot the number of reads after pre-processing.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> plotSTARstats_reads(df)
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
    '''Plot the average input read length.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> plotSTARstats_readLen(df)
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
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
    '''Plot STAR mapping statistics.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> plotSTARstats_mapping(df)
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
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
    '''Plot mapping errors including mismatches, deletions, and insertions.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> plotSTARstats_mistmatches(df)
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
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
    '''Plot chimeric reads and splices.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> plotSTARstats_chimeric(df)
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
    reads = df['                          Number of input reads |'].astype(int)
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
    '''Plot all STAR statistics including reads, read length, mapping, mismatches, and chimeric reads.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> plotSTARstats(df)
    '''
    plotSTARstats_reads(df=df, dpi=dpi)
    plotSTARstats_readLen(df=df, dpi=dpi)
    plotSTARstats_mapping(df=df, dpi=dpi)
    plotSTARstats_mistmatches(df=df, dpi=dpi)
    plotSTARstats_chimeric(df=df, dpi=dpi)


def hplotSTARstats_reads(df=pd.DataFrame(), dpi=150):
    '''Plot the number of reads after pre-processing as a horizontal bar plot.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> hplotSTARstats_reads(df)
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
    '''Plot the average input read length as a horizontal bar plot.

    :param df: DataFrame containing STAR statistics.
    :type df: pd.DataFrame
    :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
    :type dpi: int, optional

    :return: None

    :example:

    >>> hplotSTARstats_readLen(df)
    '''
    df = df.T
    labels = df.index.tolist()
    x = np.arange(len(labels))
    width = 0.3
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
        '''Plot STAR mapping statistics as a horizontal bar plot.

        :param df: DataFrame containing STAR statistics.
        :type df: pd.DataFrame
        :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
        :type dpi: int, optional

        :return: None

        :example:

        >>> hplotSTARstats_mapping(df)
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
        '''Plot mapping errors including mismatches, deletions, and insertions as a horizontal bar plot.

        :param df: DataFrame containing STAR statistics.
        :type df: pd.DataFrame
        :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
        :type dpi: int, optional

        :return: None

        :example:

        >>> hplotSTARstats_mistmatches(df)
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
        '''Plot chimeric reads and splices as a horizontal bar plot.

        :param df: DataFrame containing STAR statistics.
        :type df: pd.DataFrame
        :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
        :type dpi: int, optional

        :return: None

        :example:

        >>> hplotSTARstats_chimeric(df)
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
        '''Plot all STAR statistics as horizontal bar plots including reads, read length, mapping, mismatches, and chimeric reads.

        :param df: DataFrame containing STAR statistics.
        :type df: pd.DataFrame
        :param dpi: Resolution of the plot in dots per inch. Defaults to 150.
        :type dpi: int, optional

        :return: None

        :example:

        >>> hplotSTARstats(df)
        '''
        hplotSTARstats_reads(df=df, dpi=dpi)
        hplotSTARstats_readLen(df=df, dpi=dpi)
        hplotSTARstats_mapping(df=df, dpi=dpi)
        hplotSTARstats_mistmatches(df=df, dpi=dpi)
        hplotSTARstats_chimeric(df=df, dpi=dpi)

    def classes(df, fraction=True, dpi=300, save=None, d_color=None):
        '''Plot classes of mapped reads as a bar plot.

        :param df: DataFrame containing read classes.
        :type df: pd.DataFrame
        :param fraction: Whether to plot as fractions (percentages). Defaults to True.
        :type fraction: bool, optional
        :param dpi: Resolution of the plot in dots per inch. Defaults to 300.
        :type dpi: int, optional
        :param save: File path to save the figure. If not provided, the figure will be displayed.
        :type save: str, optional
        :param d_color: Dictionary of colors for each class. Defaults to None.
        :type d_color: dict, optional

        :return: None

        :example:

        >>> classes(df, fraction=True, dpi=300, save='classes.png', d_color={'intergenic_region':"orange", 'ncRNA':"yellow", 'protein_coding':"blue"})
        '''
        if fraction:
            df = df.groupby('type').agg("sum").div(df.sum()).multiply(100)
            df = df.drop(columns='type')
        else:
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

        for i, row in df.iterrows():
            ax1.bar(x, row, width, bottom=s0, label=i, color=d_color[i])
            s0 = s0 + row

        plt.ylabel('Cumulative mapping [%]')
        plt.xticks(x, labels, rotation=90)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        # Save the figure if a file path is provided
        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches='tight')
            
        plt.show()