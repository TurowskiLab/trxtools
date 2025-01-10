import pandas as pd
import trxtools.methods as ttm

### Functions used in desing of in vitro assays i.e. elongation, termination, etc. for RNA polymerases

################################################
#############        input

def findPrimer(seq=str(), query=str()):
    '''Find query sequence within given sequence

    :param seq: sequence
    :type seq: str
    :param query: query within searched sequences
    :type query: str

    :return: (start,stop) in python indexing
    :rtype: tuple

    :example:

    >>> findPrimer("ATCGATCG", "ATCG")
    (0, 4)
    '''

    if not query:
        primerStart = 0
    elif query:
        primerStart = seq.find(query)

    primerStop = primerStart + len(query)

    return (primerStart, primerStop)

################################################
#############        selection

def sequenceConstrain(structure=str(), stall="", RNAprimer=""):
    '''Returns sequence constrained with the 5prime RNA primer

    :param structure: sequence structure
    :type structure: str
    :param stall: single letter representing the stall sequence
    :type stall: str
    :param RNAprimer: RNA primer sequence
    :type RNAprimer: str

    :return: sequence constraints for folding algorithm
    :rtype: str
    '''

    stallComplement = {"A": "B",  # use nucleotides that are not present in stall sequence
                       "T": "V",  # IUPAC nucleotide code
                       "C": "D",
                       "G": "H",
                       "": "N"}

    length = len(structure)
    length2fill = length - len(RNAprimer)
    nt = stallComplement[stall]
    if RNAprimer:
        return RNAprimer.lower() + length2fill * nt  # lower provide strict constrain
    else:
        return length2fill * nt

def structureFile(structure=str(), stall="", RNAprimer=""):
    '''Saves structure file in current directory

    :param structure: sequence structure
    :type structure: str
    :param stall: single letter representing the stall sequence
    :type stall: str
    :param RNAprimer: RNA primer sequence
    :type RNAprimer: str

    :return: True
    :rtype: bool
    '''

    f = open("structure.txt", "w+")
    f.write(structure + "\n")  # saves structure format
    f.write(sequenceConstrain(structure=structure, stall=stall,
                              RNAprimer=RNAprimer) + "\n")  # saves nt constrains, N's if no RNAprimer
    f.close()
    return True


def stalled(seq=str(), stall="AAA", primer=str()):
    '''Finds and returns stalled sequence

    :param seq: sequence
    :type seq: str
    :param stall: stall sequence, default "AAA"
    :type stall: str
    :param primer: primer sequence
    :type primer: str

    :return: stalled sequence
    :rtype: str
    '''

    (primerStart, primerStop) = findPrimer(seq, primer)
    stall = seq.find(stall, primerStop)
    return seq[:stall]


def extruded(seq=str(), buried=13):
    '''Returns the extruded sequence by removing the buried nucleotides from the end

    :param seq: sequence
    :type seq: str
    :param buried: length of sequence buried within RNAP
    :type buried: int

    :return: extruded sequence
    :rtype: str
    '''
    return seq[:-buried]


################################################
#############        ordering

def templateDNA(seq=str(), overhang5end=str()):
    '''Returns the DNA sequence of the template strand. Takes reverse complement of RNA sequence and the 5'end DNA overhang as input.
    
    :param seq: sequence of RNA
    :type seq: str
    :param overhang5end: sequence of the 5'end DNA overhang
    :type overhang5end: str

    :return: DNA sequence of template strand
    :rtype: str
    '''
    return ttm.reverse_complement_RNA(seq) + overhang5end[::-1]


def nonTemplateDNA(seq="", primer=""):
    '''Returns the DNA sequence of the non-template strand.

    :param seq: sequence
    :type seq: str
    :param primer: primer sequence
    :type primer: str

    :return: DNA sequence of non-template strand
    :rtype: str
    '''
    seq = seq.upper()
    primer = primer.upper()

    (primerStart, primerStop) = findPrimer(seq, primer)

    return seq[primerStop + 1:]


def bulkInputIDT(data=pd.DataFrame()):
    '''Transforms dataframe with columns "template" and "non-template" and prepares table to be used as bulk input.
    Discards oligos longer than 200 nt.

    :param data: Table with columns "template" and "non-template"
    :type data: DataFrame

    :return: Table with sequences to order
    :rtype: DataFrame
    '''

    # copy data
    working = pd.DataFrame()
    working['template'] = data['template']
    working['non-template'] = data['non-template']

    # selecting DNA oligos that can be ordered as DNA oligo
    print(str(len(working[working[
                              'non-template'].str.len() > 200])) + " constructs has been discarded due to length requirement (200 nt)")
    working = working[working['non-template'].str.len() <= 200]

    # table to order
    output = pd.DataFrame()
    for i, row in working.iterrows():
        output = output.append(pd.Series(index=['sequence'], data=row["template"], name=str(i) + "_template"))
        output = output.append(pd.Series(index=['sequence'], data=row["non-template"], name=str(i) + "_non-temp"))

    # synthesis scale
    output['scale'] = "25nm"
    output['scale'][output['sequence'].str.len() > 60] = '100nm'
    output['scale'][output['sequence'].str.len() > 90] = '250nm'
    output['scale'][output['sequence'].str.len() > 100] = '4nmU'

    # purification
    output['purification'] = "STD"
    return output


def testScaffold(data=pd.DataFrame(), overhang5end="", RNA_primer="", ):
    '''Prints scaffold to test before ordering

    :param data: DataFrame with columns "template" and "non-template"
    :type data: DataFrame
    :param overhang5end: sequence of the 5'end DNA overhang
    :type overhang5end: str
    :param RNA_primer: RNA primer sequence
    :type RNA_primer: str

    :return: None
    '''
    for i, row in data.iterrows():
        template = row['template']
        nonTemplate = row['non-template']
        shift = " " * len(overhang5end)

        print(template[::-1])
        print(shift + RNA_primer + " " + nonTemplate)
        print()


def toOrder(data=pd.DataFrame(), buried="GUCUGUUUUGUGG", stallFull="AAA", afterStall="TGATCGGTAC",
            overhang5end="TGA", RNA_primer="AGGCCGAAA", bulkInput=True, test=True, lengthMax=200):
    '''Prepares sequences for ordering by applying folding functions to DNA sequences

    :param data: DataFrame with columns "sequence" and "name"
    :type data: DataFrame
    :param buried: sequence buried within RNAP, default "GUCUGUUUUGUGG"
    :type buried: str
    :param stallFull: stall sequence, default "AAA"
    :type stallFull: str
    :param afterStall: sequence after stall, default "TGATCGGTAC"
    :type afterStall: str
    :param overhang5end: sequence of the 5'end DNA overhang, default "TGA"
    :type overhang5end: str
    :param RNA_primer: RNA primer sequence, default "AGGCCGAAA"
    :type RNA_primer: str
    :param bulkInput: whether to prepare bulk input for ordering, default True
    :type bulkInput: bool
    :param test: whether to print scaffold for testing before ordering, default True
    :type test: bool
    :param lengthMax: maximum length of sequences to order, default 200
    :type lengthMax: int
    
    :return: DataFrame of sequences to order
    :rtype: DataFrame
    '''

    buried = buried.replace("T", "U")
    stallFull = stallFull.replace("T", "U")
    afterStall = afterStall.replace("T", "U")
    RNA_primer = RNA_primer.replace("T", "U")

    data = data.drop('name', axis=1)
    data['length_extruded'] = data['sequence'].str.len()

    # preparing output data
    data['stalled'] = data['sequence'].str.upper().str.replace("T", "U") + buried
    data['runoff'] = data['stalled'] + stallFull + afterStall
    data['non-template'] = data['runoff'].apply(nonTemplateDNA, primer=RNA_primer).str.replace("U", "T")
    data['template'] = data['runoff'].apply(templateDNA, overhang5end=overhang5end).str.replace("U", "T")

    if test == True:
        testScaffold(data, overhang5end=overhang5end, RNA_primer=RNA_primer)

    if bulkInput == False:
        return data

    else:
        return bulkInputIDT(data)