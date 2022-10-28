import pyBigWig
import pandas as pd
import numpy as np

def strip_BigWig_names(files=list()):
    '''_summary_

    :param files: _description_, defaults to list()
    :type files: _type_, optional
    :return: _description_
    :rtype: _type_
    '''
    #returns uniq names for *fwd.bw and *rev.bw files
    return list(set([f.replace("_fwd.bw","").replace("_rev.bw","") for f in files]))

def getSeqData(gene_name, data_path, name, gtf, ranges=0):
    '''_summary_

    :param gene_name: _description_
    :type gene_name: _type_
    :param data_path: _description_
    :type data_path: _type_
    :param name: _description_
    :type name: _type_
    :param gtf: _description_
    :type gtf: _type_
    :param ranges: _description_, defaults to 0
    :type ranges: int, optional
    :return: _description_
    :rtype: _type_
    '''
    strand, chromosome, coordinates = gtf.strand(gene_name), gtf.chromosome(gene_name), gtf.chromosomeCoordinates(gene_name)
    if strand == "+":
        bw = pyBigWig.open(data_path+name+"_fwd.bw")
        return pd.Series(bw.values(chromosome,min(coordinates)-ranges,max(coordinates)+ranges))
    if strand == "-":
        bw = pyBigWig.open(data_path+name+"_rev.bw")
        return pd.Series(bw.values(chromosome,min(coordinates)-ranges,max(coordinates)+ranges)[::-1])

def geneFromBigWig(gene_name, data_path, data_files, gtf, ranges=0,verbose=False):
    '''_summary_

    :param gene_name: _description_
    :type gene_name: _type_
    :param data_path: _description_
    :type data_path: _type_
    :param data_files: _description_
    :type data_files: _type_
    :param gtf: _description_
    :type gtf: _type_
    :param ranges: _description_, defaults to 0
    :type ranges: int, optional
    :param verbose: _description_, defaults to False
    :type verbose: bool, optional
    :return: _description_
    :rtype: _type_
    '''
    df_t1 = pd.DataFrame()
    df_t1["nucleotide"] = "_".join(gtf.genomicSequence(gene_name,ranges=ranges)).split("_")
    for name in strip_BigWig_names(data_files):
        if verbose==True: print(name)
        df_t1[name] = getSeqData(gene_name, data_path, name, gtf, ranges=ranges)
    return df_t1

def foldingFromBigWig(gene_name, data_path, data_files, gtf, ranges=0,range5end=0,offset=15):
    '''_summary_

    :param gene_name: _description_
    :type gene_name: _type_
    :param data_path: _description_
    :type data_path: _type_
    :param data_files: _description_
    :type data_files: _type_
    :param gtf: _description_
    :type gtf: _type_
    :param ranges: _description_, defaults to 0
    :type ranges: int, optional
    :param range5end: _description_, defaults to 0
    :type range5end: int, optional
    :param offset: _description_, defaults to 15
    :type offset: int, optional
    :return: _description_
    :rtype: _type_
    '''
    df_t1 = pd.DataFrame()
    df_t1["nucleotide"] = "_".join(gtf.genomicSequence(gene_name,ranges=ranges)).split("_")
    for name in strip_BigWig_names(data_files):
        df_t1[name] = getSeqData(gene_name, data_path, name, gtf, ranges=ranges)
        df_t1[name][0:ranges-range5end+int(name.replace("w","").replace("_dG",""))] = np.nan #remove folding for the 5'end
        df_t1[name+"_add"+str(offset)+"nt"] = pd.Series((offset*[np.nan])+df_t1[name].tolist()) #shifting to the extrustion point
    return df_t1