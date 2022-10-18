import re
import pandas as pd

### SECONDARY (SecondarySctructure)
# functions to deal with RNA secondary structure files (default it Vienna)

################################################
#############        secondary structure trxtools
### all tools written using genomic positions cointer from 1 (not from 0 as python default)

def checkVienna(sequence="", vienna=""):
    '''Validates integrity of vienna file

    :param sequence: str
    :param vienna: str
    :return: True if pass
    '''
    if len(sequence) != len(vienna): return False

    a = vienna.count("(")
    b = vienna.count(")")
    c = vienna.count(".")

    if a != b: return False
    if (a + b + c) != len(sequence): return False
    return True


def Lstem(vienna=""):
    '''Returns list of positions where "(" is found using coordinates {1,inf}

    :param vienna: str
    :return: list
    '''
    return [i.start() + 1 for i in re.finditer("\(", vienna)]


def Rstem(vienna=""):
    '''Returns list of positions where ")" is found using coordinates {1,inf}

    :param vienna: str
    :return: list
    '''
    return [i.start() + 1 for i in re.finditer("\)", vienna)]


def loops(vienna=""):
    '''
    Returns first positions outside the loop i.e. ".((....))." returns [(3,8)]

    :param vienna: vienna
    :return: list of tuples
    '''
    # TO DO: check does loop is correct
    loops = [l.span() for l in re.finditer(r"\((\.+)\)", vienna)]  # list of loops (position of limiting base pairing)
    # genomic coordinates {1,inf}
    return [(l + 1, k) for (l, k) in loops]


def test(vienna="",sequence="", loops=None, stems=None, multistems=None, linkers=None):
    '''Prints vienna with given features

    :param vienna: str
    :param loops: list of tuples (option)
    :param stems: list of tuples (option)
    :param multistems: list (option)
    :param linkers: list (option)
    :return: None
    '''
    # genomic coordinates {1,inf}
    scale = [" "] * len(vienna)
    scale[0] = "1"
    for i in range(1, len(vienna) - 2):
        if i % 10 == 0:
            scale[i - 1] = str(int(i / 10))
            scale[i] = str(0)
    print("".join(scale))

    print(vienna)
    if sequence:
        print(sequence)

    if loops:
        printMarks = ["_"] * len(vienna)
        for i, j in loops:
            printMarks[i - 1] = "O"
            printMarks[j - 1] = "O"
    if stems:
        for i, j in stems:
            printMarks[i - 1] = "S"
            printMarks[j - 1] = "S"
    if multistems:
        for i in multistems:
            printMarks[i - 1] = "M"
    if linkers:
        for i, j in linkers:
            printMarks[i - 1] = "L"
            printMarks[j - 1] = "L"

    print("".join(printMarks))
    return None


def loopStems(vienna="", sequence="", loopsList=None, testPrint=False):
    '''Returns postions of stem of single hairpins and multiloop stems. Use coordinates {1:inf}.
    Warninig: tested with single multiloop stems only

    :param vienna: str
    :param loopsList: list (option)
    :param testPrint: boolean to default=False
    :return: list, list (stems: list of tuples; multistems: list)
    '''
    # genomic coordinates {1,inf}
    if not loopsList:
        loopsList = loops(vienna)

    stems = []  # output
    multistems = []

    stemL = Lstem(vienna)  # all left halves of stems
    stemR = Rstem(vienna)  # all right halves of stems

    for loop in loopsList:
        start, stop = loop[0], loop[1] - 1  # loop

        stemLpotential = [i for i in stemL if i <= start]  # stem elements upstream the loop
        stemLalt = [i for i in stemR if i <= start]  # alternative stems upstream the loop
        stemRpotential = [i for i in stemR if i >= stop]  # stem elements downstream the loop
        stemRalt = [i for i in stemL if i >= stop]  # alternative stems downstream the loop

        # last stem element upstream
        if not stemLalt:
            start = min(stemL)
        else:
            stemLpotential = [i for i in stemLpotential if i >= max(stemLalt)]
            start = min(stemLpotential)

        # last stem element downstream
        if not stemRalt:
            stop = max(stemR)
        else:
            stemRpotential = [i for i in stemRpotential if i <= min(stemRalt)]
            stop = max(stemRpotential)

        # detect and exclude potential multiloop stems

        ##WARNING works only with single multiloop structures (optional TODO)
        stemUp = vienna[start - 1:loop[0]]
        stemDown = vienna[loop[1] - 1:stop]
        if stemUp.count("(") != stemDown.count(")"):
            if stemUp.count("(") > stemDown.count(")"):
                multistems.append(start)

                while vienna[start - 1:loop[0]].count("(") != vienna[loop[1] - 1:stop].count(")"):
                    start += 1
                multistems.append(start - 1)

                while len(vienna[start - 1:loop[0]]) != len(vienna[start - 1:loop[0]].strip(".")):
                    start += 1

            if stemUp.count("(") < stemDown.count(")"):
                last = stop

                while vienna[start - 1:loop[0]].count("(") != vienna[loop[1] - 1:stop].count(")"):
                    stop -= 1
                multistems.append(stop + 1)
                multistems.append(last)

                while len(vienna[loop[1]:stop]) != len(vienna[loop[1]:stop].strip(".")):
                    stop -= 1

        stems.append((start, stop))

    # print to check
    if testPrint == True:
        test(vienna, sequence, loopsList, stems, multistems)

    return stems, multistems


def vienna2format(vienna="", sequence="", loopsList=None, stemsList=None, multistemsList=None, testPrint=False):
    #TODO find name to this "format"
    '''Converts vienna format to letters: O - loop, S - stem, M - multiloop stem and L - linker

    :param vienna: str
    :param loopsList: list (optional)
    :param stemsList: list (optional)
    :param multistemsList: list (optional)
    :param testPrint: defauls=False
    :return: str in "format"
    '''
    # genomic coordinates {1,inf}

    if not loopsList:
        loopsList = loops(vienna)

    if not stemsList or not multistemsList:
        stemsList, multistemsList = loopStems(vienna)

    output = list(vienna)

    # overwrite stems
    for (s_start, s_stop) in stemsList:
        s_len = s_stop - s_start + 1
        output[s_start - 1:s_start + s_len - 1] = ['S'] * s_len

    # overwrite loops
    for (l_start, l_stop) in loopsList:
        l_len = l_stop - l_start - 1
        output[l_start:l_start + l_len] = ['O'] * l_len

    # overwrite multistems
    if multistemsList:
        if len(multistemsList) != 4:
            print("More than one multiloop detected")
            pass
        elif len(multistemsList) == 4:
            output[multistemsList[0] - 1:multistemsList[1]] = ['M'] * (multistemsList[1] - multistemsList[0] + 1)
            output[multistemsList[2] - 1:multistemsList[3]] = ['M'] * (multistemsList[3] - multistemsList[2] + 1)

    output = "".join(output)

    if testPrint == True:
        scale = [" "] * len(vienna)
        scale[0] = "1"
        scale2 = [" "] * len(vienna)
        scale2[0] = "|"
        for i in range(1, len(vienna) - 2):
            if i % 10 == 0:
                scale[i - 1] = str(int(i / 10))
                scale[i] = str(0)
                scale2[i - 1] = "|"
        print("".join(scale))
        print("".join(scale2))

        print(vienna)
        if sequence:
            print(sequence)
        print("".join(scale2))
        print(output.replace(".", "L"))

    return output.replace(".", "L")


def substructures(vienna="", sequence=""):
    '''list sub-structures of the given structure

    :param vienna: str
    :param sequence: str
    :return: Series
    '''
    if len(vienna) != len(sequence): return None

    stems, mstems = loopStems(vienna)

    seqList = []
    nameList = []

    for (s1, s2) in stems:
        seqList.append(sequence[s1 - 1:s2])
        nameList.append("stem_" + str(s1))

    if mstems:
        seqList.append(sequence[mstems[0] - 1:mstems[3]])
        nameList.append("mstem_" + str(mstems[0]))

    return pd.Series(data=seqList, index=nameList)