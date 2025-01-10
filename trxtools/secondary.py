import re
import pandas as pd

### SECONDARY (SecondarySctructure)
# functions to deal with RNA secondary structure files (default it Vienna)

################################################
#############        secondary structure trxtools
### all tools written using genomic positions cointer from 1 (not from 0 as python default)

def checkVienna(sequence="", vienna=""):
    '''Validates the integrity of a Vienna file by checking the sequence and structure lengths and the balance of parentheses.

    :param sequence: RNA sequence
    :type sequence: str
    :param vienna: Vienna format secondary structure
    :type vienna: str

    :return: True if the Vienna file is valid, False otherwise
    :rtype: bool

    :example:

    >>> checkVienna("GCAU", "((..))")
    False
    '''
    if len(sequence) != len(vienna): return False

    a = vienna.count("(")
    b = vienna.count(")")
    c = vienna.count(".")

    if a != b: return False
    if (a + b + c) != len(sequence): return False
    return True


def Lstem(vienna=""):
    '''Returns a list of positions where "(" is found using coordinates {1,inf}.

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :return: List of positions with "("
    :rtype: list

    Example:
    >>> Lstem("((..))")
    [1, 2]
    '''
    return [i.start() + 1 for i in re.finditer("\(", vienna)]


def Rstem(vienna=""):
    '''Returns a list of positions where ")" is found using coordinates {1,inf}.

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :return: List of positions with ")"
    :rtype: list

    Example:
    >>> Rstem("((..))")
    [5, 6]
    '''
    return [i.start() + 1 for i in re.finditer("\)", vienna)]


def loops(vienna=""):
    '''Returns the first positions outside the loop, e.g., ".((....))." returns [(3,8)].

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :return: List of tuples with loop positions
    :rtype: list

    Example:
    >>> loops(".((....)).")
    [(3, 8)]
    '''
    loops = [l.span() for l in re.finditer(r"\((\.+)\)", vienna)]
    return [(l + 1, k) for (l, k) in loops]


def test(vienna="", sequence="", loops=None, stems=None, multistems=None, linkers=None):
    '''Prints the Vienna structure with given features marked.

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :param sequence: RNA sequence (optional)
    :type sequence: str
    :param loops: List of tuples with loop positions (optional)
    :type loops: list
    :param stems: List of tuples with stem positions (optional)
    :type stems: list
    :param multistems: List of multistem positions (optional)
    :type multistems: list
    :param linkers: List of linker positions (optional)
    :type linkers: list
    :return: None

    Example:
    >>> test("((..))", "GCAU", loops=[(3, 4)])
    1         
    ((..))
    GCAU
    __O_O
    '''
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
    '''Returns positions of stems of single hairpins and multiloop stems using coordinates {1:inf}.

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :param sequence: RNA sequence (optional)
    :type sequence: str
    :param loopsList: List of loop positions (optional)
    :type loopsList: list
    :param testPrint: Boolean to print test output, default=False
    :type testPrint: bool

    :return: List of stem positions and list of multistem positions
    :rtype: tuple

    :example:

    >>> loopStems("((..))", "GCAU")
    ([(1, 6)], [])
    '''
    if not loopsList:
        loopsList = loops(vienna)

    stems = []
    multistems = []

    stemL = Lstem(vienna)
    stemR = Rstem(vienna)

    for loop in loopsList:
        start, stop = loop[0], loop[1] - 1

        stemLpotential = [i for i in stemL if i <= start]
        stemLalt = [i for i in stemR if i <= start]
        stemRpotential = [i for i in stemR if i >= stop]
        stemRalt = [i for i in stemL if i >= stop]

        if not stemLalt:
            start = min(stemL)
        else:
            stemLpotential = [i for i in stemLpotential if i >= max(stemLalt)]
            start = min(stemLpotential)

        if not stemRalt:
            stop = max(stemR)
        else:
            stemRpotential = [i for i in stemRpotential if i <= min(stemRalt)]
            stop = max(stemRpotential)

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

    if testPrint:
        test(vienna, sequence, loopsList, stems, multistems)

    return stems, multistems


def vienna2format(vienna="", sequence="", loopsList=None, stemsList=None, multistemsList=None, testPrint=False):
    '''Converts Vienna format to a custom format with letters: O - loop, S - stem, M - multiloop stem, and L - linker.

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :param sequence: RNA sequence (optional)
    :type sequence: str
    :param loopsList: List of loop positions (optional)
    :type loopsList: list
    :param stemsList: List of stem positions (optional)
    :type stemsList: list
    :param multistemsList: List of multistem positions (optional)
    :type multistemsList: list
    :param testPrint: Boolean to print test output, default=False
    :type testPrint: bool

    :return: Custom format string
    :rtype: str

    :example:

    >>> vienna2format("((..))", "GCAU")
    'SSLLSS'
    '''
    if not loopsList:
        loopsList = loops(vienna)

    if not stemsList or not multistemsList:
        stemsList, multistemsList = loopStems(vienna)

    output = list(vienna)

    for (s_start, s_stop) in stemsList:
        s_len = s_stop - s_start + 1
        output[s_start - 1:s_start + s_len - 1] = ['S'] * s_len

    for (l_start, l_stop) in loopsList:
        l_len = l_stop - l_start - 1
        output[l_start:l_start + l_len] = ['O'] * l_len

    if multistemsList:
        if len(multistemsList) != 4:
            print("More than one multiloop detected")
            pass
        elif len(multistemsList) == 4:
            output[multistemsList[0] - 1:multistemsList[1]] = ['M'] * (multistemsList[1] - multistemsList[0] + 1)
            output[multistemsList[2] - 1:multistemsList[3]] = ['M'] * (multistemsList[3] - multistemsList[2] + 1)

    output = "".join(output)

    if testPrint:
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
    '''Lists sub-structures of the given structure.

    :param vienna: Vienna format secondary structure
    :type vienna: str
    :param sequence: RNA sequence
    :type sequence: str

    :return: Series of sub-structures
    :rtype: pd.Series

    :example:
    
    >>> substructures("((..))", "GCAU")
    stem_1    GCAU
    dtype: object
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