import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from trxtools.plotting import clusterClusterMap

################################################
#############       statistics and analysis

def normalize(df=pd.DataFrame, log2=False, pseudocounts=0.1, CPM=True):
    '''Normalize the DataFrame by adding pseudocounts, dividing by the sum, and optionally applying log2 transformation and CPM scaling.

    :param df: DataFrame to be normalized
    :type df: pandas.DataFrame
    :param log2: Whether to apply log2 transformation, defaults to False
    :type log2: bool, optional
    :param pseudocounts: Value to add to each element to avoid log of zero, defaults to 0.1
    :type pseudocounts: float, optional
    :param CPM: Whether to scale the data to counts per million, defaults to True
    :type CPM: bool, optional

    :return: Normalized DataFrame
    :rtype: pandas.DataFrame
    '''

    df = df.add(pseudocounts)
    df = df / df.sum()
    if CPM==True:
        df = df.multiply(1000000)
    
    if log2==True:
        return df.apply(np.log2)
    else:
        return df

def quantileCategory(s1=pd.Series(dtype=float), q=4):
    '''Quantile-based discretization function based on pandas.qcut function.

    :param s1: Series to be discretized
    :type s1: pandas.Series
    :param q: Number of quantiles, defaults to 4
    :type q: int, optional

    :return: Series with quantile categories
    :rtype: pandas.Series

    :example:

    >>> quantileCategory(pd.Series([1, 2, 3, 4, 5]), q=2)
    0       0
    1       0
    2       1
    3       1
    4       1
    dtype: int64
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

    :param data: DataFrame containing the input data
    :type data: pandas.DataFrame
    :param n_components: Number of principal components to compute, defaults to 2
    :type n_components: int, optional

    :return: Tuple consisting of DataFrame with PCA results and a list of explained variance ratios
    :rtype: tuple

    :example:

    >>> runPCA(pd.DataFrame([[1, 2], [3, 4], [5, 6]]), n_components=2)
    (          PC1       PC2
    0 -4.242641  0.000000
    1  0.000000  0.000000
    2  4.242641  0.000000, [100.0, 0.0])
    '''

    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(data)
    principalDf = pd.DataFrame(data=principalComponents,
                               columns=['PC' + str(i + 1) for i in np.arange(0, n_components)])

    values = pca.explained_variance_ratio_
    principalDf.index = data.index

    return principalDf, [round(i*100,2) for i in values.tolist()]

def addCluster(df=pd.DataFrame(), n=10):
    '''Assigns n clusters to the data using KMeans algorithm

    :param df: DataFrame containing the data to be clustered
    :type df: pandas.DataFrame
    :param n: Number of clusters to form, defaults to 10
    :type n: int, optional

    :return: DataFrame with an additional 'cluster' column indicating cluster assignment
    :rtype: pandas.DataFrame

    :example:

    >>> df = pd.DataFrame({'x': [1, 2, 3, 4, 5], 'y': [5, 4, 3, 2, 1]})
    >>> addCluster(df, n=2)
       x  y  cluster
    0  1  5        1
    1  2  4        1
    2  3  3        0
    3  4  2        0
    4  5  1        0
    '''

    if 'cluster' in df.columns:
        df = df.drop('cluster', 1)
    kmeans = KMeans(n_clusters=n, random_state=0).fit(df)
    df['cluster'] = kmeans.labels_
    df = df.sort_values('cluster', ascending=True)

    # summary
    clusterClusterMap(df)

    return df
