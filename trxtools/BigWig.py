import pyBigWig
import pandas as pd
import numpy as np
import warnings

# BigWig functions for extracting data from BigWig files. 
# DEPRECATION WARNING: These functions will be removed in a future version.

def strip_BigWig_names(files=list()):
    '''Remove _fwd.bw and _rev.bw from the list of BigWig files

    :param files: list of filenames, defaults to list()
    :type files: list, required
    :return: list of unique names
    :rtype: list

    >>> strip_BigWig_names(["file1_fwd.bw", "file1_rev.bw", "file2_fwd.bw", "file2_rev.bw"])
    ['file2', 'file1']
    '''

    warnings.warn("The strip_BigWig_names function is deprecated and will be removed in a future version.", DeprecationWarning)
    #returns uniq names for *fwd.bw and *rev.bw files
    return list(set([f.replace("_fwd.bw","").replace("_rev.bw","") for f in files]))

def getSeqData(gene_name, data_path, name, gtf, ranges=0):
    '''Load BigWig data for a gene/chromosome region and return as pandas Series.

    :param gene_name: Name of the gene/chromosome to retrieve data for.
    :type gene_name: str
    :param data_path: Path to the directory containing BigWig files.
    :type data_path: str
    :param name: Base name of the BigWig files (without strand suffix and extension). Can be obtained using strip_BigWig_names().
    :type name: str
    :param gtf: GTF object that provides gene information such as strand, chromosome, and coordinates.
    :type gtf: object
    :param ranges: Number of bases to extend the region on both sides, defaults to 0.

    :return: BigWig data for the specified gene or region as a pandas Series.
    :rtype: pandas.Series

    :example:

    >>> gtf = GTF('/path/to/gtf_file.gtf')
    >>> data = getSeqData('BRCA1', '/path/to/bigwig/', 'sample', gtf, ranges=100)
    >>> print(data)
    
    '''
    warnings.warn("The getSeqData function is deprecated and will be removed in a future version.", DeprecationWarning)

    strand, chromosome, coordinates = gtf.strand(gene_name), gtf.chromosome(gene_name), gtf.chromosomeCoordinates(gene_name)
    if strand == "+":
        bw = pyBigWig.open(data_path+name+"_fwd.bw")
        return pd.Series(bw.values(chromosome,min(coordinates)-ranges,max(coordinates)+ranges))
    if strand == "-":
        bw = pyBigWig.open(data_path+name+"_rev.bw")
        return pd.Series(bw.values(chromosome,min(coordinates)-ranges,max(coordinates)+ranges)[::-1])

def geneFromBigWig(gene_name, data_path, data_files, gtf, ranges=0,verbose=False):
    '''Extracts nucleotide sequences and associated data for a given gene from BigWig files.
    
    :param gene_name: Name of the gene to extract data for.
    :type gene_name: str
    :param data_path: Path to the directory containing BigWig files.
    :type data_path: str
    :param data_files: List of BigWig file names to extract data from.
    :type data_files: list
    :param gtf: GTF file object containing gene annotations.
    :type gtf: GTF
    :param ranges: Range of nucleotides to extract, defaults to 0.
    :param verbose: If True, prints the name of each BigWig file being processed, defaults to False.
    
    :return: DataFrame containing nucleotide sequences and associated data.
    :rtype: pd.DataFrame

    :example:
    
    >>> gtf = GTF("/path/to/annotations.gtf")
    >>> df = geneFromBigWig("BRCA1", "/path/to/bigwig/", ["file1.bw", "file2.bw"], gtf, ranges=1000, verbose=True)
    >>> print(df.head())
    '''    
    
    warnings.warn("The geneFromBigWig function is deprecated and will be removed in a future version.", DeprecationWarning)
    df_t1 = pd.DataFrame()
    df_t1["nucleotide"] = "_".join(gtf.genomicSequence(gene_name,ranges=ranges)).split("_")
    for name in strip_BigWig_names(data_files):
        if verbose==True: print(name)
        df_t1[name] = getSeqData(gene_name, data_path, name, gtf, ranges=ranges)
    return df_t1

def foldingFromBigWig(gene_name, data_path, data_files, gtf, ranges=0,range5end=0,offset=15):
    '''Extracts folding energy for a given gene from BigWig files, 
    and applies a folding offset to the data.

    :param gene_name: Name of the gene to extract data for.
    :type gene_name: str
    :param data_path: Path to the directory containing BigWig files.
    :type data_path: str
    :param data_files: List of BigWig file names to extract data from.
    :type data_files: list
    :param gtf: GTF file object containing gene annotations.
    :type gtf: GTF
    :param ranges: Range of nucleotides to extract, defaults to 0.
    :type ranges: int, optional
    :param range5end: Number of bases to exclude from the 5' end, defaults to 0.
    :type range5end: int, optional
    :param offset: Number of bases to shift the data by, defaults to 15.
    :type offset: int, optional
    
    :return: DataFrame containing nucleotide sequences and associated data with applied folding offset.
    :rtype: pd.DataFrame

    :example:

    >>> gtf = GTF("/path/to/annotations.gtf")
    >>> df = foldingFromBigWig("BRCA1", "/path/to/bigwig/", ["file1.bw", "file2.bw"], gtf, ranges=1000, range5end=100, offset=15)
    >>> print(df.head())
    '''
    warnings.warn("The foldingFromBigWig function is deprecated and will be removed in a future version.", DeprecationWarning)
    df_t1 = pd.DataFrame()
    df_t1["nucleotide"] = "_".join(gtf.genomicSequence(gene_name,ranges=ranges)).split("_")
    for name in strip_BigWig_names(data_files):
        df_t1[name] = getSeqData(gene_name, data_path, name, gtf, ranges=ranges)
        df_t1[name][0:ranges-range5end+int(name.replace("w","").replace("_dG",""))] = np.nan #remove folding for the 5'end
        df_t1[name+"_add"+str(offset)+"nt"] = pd.Series((offset*[np.nan])+df_t1[name].tolist()) #shifting to the extrustion point
    return df_t1