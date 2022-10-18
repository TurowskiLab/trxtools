import pandas as pd
import trxtools.methods as ttm

### ASSAYS (ElongationAssay)

################################################
#############        input

def findPrimer(seq=str(), query=str()):
    '''Find query sequence within given sequence

    :param seq: str() containing sequence
    :param query: str() within searched sequences
    :return: (start,stop)
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

    :param structure: str()
    :param stall: str() single letter
    :param RNAprimer: str()
    :return: str() with sequence constrains for folding algorithm
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

    :param structure: str()
    :param stall: str()
    :param RNAprimer: str()
    :return: True
    '''

    f = open("structure.txt", "w+")
    f.write(structure + "\n")  # saves structure format
    f.write(sequenceConstrain(structure=structure, stall=stall,
                              RNAprimer=RNAprimer) + "\n")  # saves nt constrains, N's if no RNAprimer
    f.close()
    return True


def stalled(seq=str(), stall="AAA", primer=str()):
    '''Finds and returns stalled sequence

    :param seq: str()
    :param stall: str() default "AAA"
    :param primer: str()
    :return: stalled sequence, str
    '''

    (primerStart, primerStop) = findPrimer(seq, primer)
    stall = seq.find(stall, primerStop)
    return seq[:stall]


def extruded(seq=str(), buried=13):
    '''
    :param seq: str()
    :param buried: int() length of sequence buried within RNAP
    :return: Extruded sequence, str
    '''
    return seq[:-buried]


################################################
#############        ordering

def templateDNA(seq=str(), overhang5end=str()):
    '''
    :param seq: str() sequence of RNA
    :param overhang5end: str() sequence of the 5'end DNA ovethang
    :return: DNA sequence of template strand, str
    '''
    return ttm.reverse_complement_RNA(seq) + overhang5end[::-1]


def nonTemplateDNA(seq="", primer=""):
    '''
    :param seq: str() sequence
    :param primer: str() sequence
    :return: DNA sequence of non-template strand, str
    '''
    seq = seq.upper()
    primer = primer.upper()

    (primerStart, primerStop) = findPrimer(seq, primer)

    return seq[primerStop + 1:]


def bulkInputIDT(data=pd.DataFrame()):
    '''Tranforms dataframe with colums "template" and "non-template" and prepares table to be used as bulk input
    discards oligos longer than 200 nt.

    :param data: DataFrame()
    :return: DataFrame with sequences to order
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

    :param data: DataFrame()
    :param overhang5end: str() overhang to shift sequences
    :param RNA_primer: str() sequence
    :return:
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
    '''Transfers folding function to DNA sequences

    :param data: DataFrame()
    :param buried: str() sequence
    :param stallFull: str() sequence
    :param afterStall: str() sequence
    :param overhang5end: str() sequence
    :param RNA_primer: str() sequence
    :param bulkInput: boolean() default True
    :param test: boolean() default True
    :param lengthMax: int() default 200
    :return: DataFrame of sequences to order
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