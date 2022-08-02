import pandas as pd
import numpy as np
import os, collections, shutil, re
import TTools as tt

### support functions ###
def groupCIGAR(cigar_string=""):
    '''Split CIGAR string to list of tuples.

    :param cigar_string: CIGAR string ""
    :type cigar_string: str
    :return: list of tuples
    :rtype: list
    '''

    return re.findall(r'(\d+)([A-Z]{1})', cigar_string)


def tostripCIGARfive(match=[]):
    '''Nucleotides without alignment at the 5' end of the read

    :param match: output from groupCIGAR (list of tuples)
    :type match: list
    :return: number of substituted nucleotides
    :rtype: int
    '''
    if "S" in match[0]:
        return int(match[0][0])
    else:
        return None


def tostripCIGARthree(match=[]):
    '''Nucleotides without alignment at the 3' end of the read

    :param match: output from groupCIGAR (list of tuples)
    :type match: list
    :return: number of substituted nucleotides
    :rtype: int
    '''
    if "S" in match[-1]:
        return int(match[-1][0])
    else:
        return None


def stripSubstitutions(match):
    '''Strip substutiotns on both ends

    :param match: output from groupCIGAR (list of tuples)
    :type match: list
    :return: list of tuples
    :rtype: list
    '''
    if "S" in match[0]:
        match = match[1:]
    if "S" in match[-1]:
        match = match[:-1]
    return match


def countRead(i=tuple()):
    '''Takes tuple (position,CIGARstring) and returns list of mapped positions

    :param i: tuple (first position of the read ,CIGAR string)
    :type i: tuple
    :return: list of mapped positions
    :rtype: np.array
    '''
    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    strip_shift = tostripCIGARfive(match)
    if not strip_shift: strip_shift = 0
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    output = np.arange(length, dtype=int)
    return output + position + strip_shift


def countMiddle(i=tuple(), expand=0):
    '''Takes tuple (position,CIGARstring) and returns list with the middle of mapped read

    :param i: tuple (first position of the read ,CIGAR string)
    :type i: tuple
    :param expand: number of nucleotides to expand for each side, defaults to 0
    :type expand: int, optional
    :return: list of mapped positions
    :rtype: np.array
    '''
    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    strip_shift = tostripCIGARfive(match)
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    mid = round(length / 2) + strip_shift
    if expand > 0:
        output = np.arange(mid - expand, mid + expand)
        return output + position + strip_shift
    else:
        return np.arange(mid, mid + 1) + position + strip_shift


def countDeletion(i=tuple(), expand=0):
    '''Takes tuple (position,CIGARstring) and returns list of mapped deletions. Each deletion is counted once, but expand parameter can merge some longer deletions (common expansion).

    :param i: tuple (first position of the read ,CIGAR string)
    :type i: tuple
    :param expand: number of nucleotides to expand for each side, defaults to 0
    :type expand: int, optional
    :return: list of mapped positions
    :rtype: np.array
    '''
    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    output = np.arange(length, dtype=int)
    outputMask = np.zeros(length, dtype=int)
    cumsum = 0
    for i, x in match:
        i = int(i)
        if x == 'D' and expand == 0:
            outputMask[cumsum:cumsum + i] = 1
        elif x == 'D' and expand > 0:
            outputMask[cumsum - expand:cumsum + i + expand] = 1
        cumsum += i
    output = (output + position) * outputMask
    return output[output != 0]


### main functions ###
def transcript2profile(l=[], length=int()):
    '''Takes list of tuples position,CIGARstring) and generates profile. Works with entire reads, for both strands.'''
    list_hits = [countRead(i) for i in l]  # list of lists
    hits = [item for sublist in list_hits for item in sublist]  # transforms to flat list
    # returning the profile
    profile = pd.Series(collections.Counter(hits))  # faster that using zip and numpy
    return profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts


def transcript2profileDeletions(l=[], expand=0, length=int()):
    '''Takes list of tuples position,CIGARstring) and generates profile. Not tested for MINUS strand.'''
    list_hits = [countDeletion(i, expand=expand) for i in l]  # list of lists
    hits = [item for sublist in list_hits for item in sublist]  # transforms to flat list
    # returning the profile
    profile = pd.Series(collections.Counter(hits))  # faster that using zip and numpy
    return profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts



def reads2profile(name=str(), dirPath=str(), df_details=pd.DataFrame()):
    '''Takes list of aligned reads and transform to profile for each transctipt. Tested only for PLUS strand.
    '''
    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']
    df_input = pd.read_csv(dirPath + "/" + name + ".tab", sep="\t", names=cols)

    output_df = pd.DataFrame()
    log = []

    for n, df in df_input.groupby('name'):
        try:
            details = df_details.loc[n]
            length = details['transcript_length']
            l = list(zip(df['position'], df['CIGAR']))
            s1 = transcript2profile(l, length=length)
            output_df = output_df.append(s1.rename(n))

            log.append(n + " - profile generated successfully")
        except:
            log.append(n + " - profile FAILED")

    return output_df, log


def reads2profileDeletions(name=str(), dirPath=str(), df_details=pd.DataFrame(), expand=5):
    '''Tested only for PLUS strand.'''
    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']
    df_input = pd.read_csv(dirPath + "/" + name + ".tab", sep="\t", names=cols)

    output_df = pd.DataFrame()
    outputDel_df = pd.DataFrame()
    outputDelExt_df = pd.DataFrame()
    log = []

    for n, df in df_input.groupby('name'):
        try:
            details = df_details.loc[n]
            length = details['transcript_length']
            l = list(zip(df['position'], df['CIGAR']))

            # reads
            s1 = transcript2profile(l, length=length)
            output_df = output_df.append(s1.rename(n))

            # deletions
            s2 = transcript2profileDeletions(l, length=length)
            outputDel_df = outputDel_df.append(s2.rename(n))

            # deletions ext
            s3 = transcript2profileDeletions(l, length=length, expand=expand)
            outputDelExt_df = outputDelExt_df.append(s3.rename(n))

            log.append(n + " - profile generated successfully")
        except:
            log.append(n + " - profile FAILED")

    return output_df, outputDel_df, outputDelExt_df, log


def sam2profiles(filename="", path='', geneList=[], toClear='', df_details=pd.DataFrame(),
                 deletions=False, expand=5, pickle=False,chunks=0):
    '''Tested only for PLUS strand.'''
    # making working directory
    name = filename.replace(".sam", "")
    if toClear:
        name = name.replace(toClear, "")
    timestamp = tt.methods.timestampRandomInt()
    dirPath = path + name + timestamp
    os.makedirs(dirPath)

    # tempfiles
    ##list
    if chunks > 0:
        chunkList = [i for i in range(0,len(geneList)+1,chunks)]
    else:
        chunkList = [0]
        chunks = len(geneList)+1

    genes = pd.DataFrame(pd.Series(geneList))
    for i in chunkList:
        geneListFileName = "/geneList_" + str(i) + ".tab"
        genes[i:i+chunks].to_csv(dirPath + geneListFileName, sep='\t', index=False, header=False)

    ##reads
    os.chdir(dirPath)
    for i in chunkList:
        geneListFileName = "geneList_" + str(i) + ".tab"
        #FLAG 0 and 256 for single end reads, aligned, forward match
        command = "grep -v ^@ ../" + filename + " | grep -f "+geneListFileName+\
                  " | awk -F'\t' 'BEGIN{OFS = FS} $2==0||$2==256{print $2,$3, $4, $6, $10, $12}' > " + name+"_" + str(i) + ".tab"
        tt.methods.bashCommand(command)
    tt.methods.bashCommand("cat "+name+"* > " + name+".tab")

    print("Reads selected.")

    # get profiles
    if deletions == False:
        df_profiles, log = reads2profile(name=name, dirPath=dirPath, df_details=df_details)
        #save output
        if pickle==True:
            df_profiles.to_pickle(path + name + "_PROFILES_reads.pcl")
        elif pickle==False:
            df_profiles.to_csv(path + name + "_PROFILES_reads.csv")
        #save log
        with open(path + name + "_PROFILES_reads.log", "w") as log_file:
            for row in log:
                log_file.write(str(row) + '\n')

    elif deletions == True:
        output_df, outputDel_df, outputDelExt_df, log = reads2profileDeletions(name=name, dirPath=dirPath,
                                                                               df_details=df_details, expand=expand)
        # save output
        if pickle==True:
            output_df.to_pickle(path + name + "_PROFILES_reads.pcl")
            outputDel_df.to_pickle(path + name + "_PROFILES_deletions_expand0.pcl")
            outputDelExt_df.to_pickle(path + name + "_PROFILES_deletions_expand" + str(expand) + ".pcl")
        elif pickle==False:
            output_df.to_csv(path + name + "_PROFILES_reads.csv")
            outputDel_df.to_csv(path + name + "_PROFILES_deletions_expand0.csv")
            outputDelExt_df.to_csv(path + name + "_PROFILES_deletions_expand" + str(expand) + ".csv")
        #save log
        with open(path + name + "_PROFILES_reads.log", "w") as log_file:
            for row in log:
                log_file.write(str(row) + '\n')

    # clean
    os.chdir(path)
    shutil.rmtree(name + timestamp)

    print("Done.")

def reads2genome(name=str(), dirPath=str(), df_details=pd.DataFrame()):
    '''Works for both strands, extensive testing needed.'''

    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']
    df_input_fwd = pd.read_csv(dirPath + "/" + name + "_fwd.tab", sep="\t", names=cols)
    df_input_rev = pd.read_csv(dirPath + "/" + name + "_rev.tab", sep="\t", names=cols)

    output_df_fwd = pd.DataFrame()
    output_df_rev = pd.DataFrame()
    log = []

    for n, df in df_input_fwd.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1 = transcript2profile(l, length=length)
            output_df_fwd = output_df_fwd.append(s1.rename(n))

            log.append(n + " - FWD profile generated successfully")
        except:
            log.append(n + " - FWD profile FAILED")
            
    for n, df in df_input_rev.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1 = transcript2profile(l, length=length)
            output_df_rev = output_df_rev.append(s1.rename(n))

            log.append(n + " - REV profile generated successfully")
        except:
            log.append(n + " - REV profile FAILED")

    return output_df_fwd, output_df_rev, log


def sam2genome(filename="", path='', geneList=[], toClear='', df_details=pd.DataFrame(),
                pickle=False,chunks=0):
    '''Works for both strands, extensive testing needed.'''
    # making working directory
    name = filename.replace(".sam", "")
    if toClear:
        name = name.replace(toClear, "")
    timestamp = tt.methods.timestampRandomInt()
    dirPath = path + name + timestamp
    os.makedirs(dirPath)

    #geneList = chromosome list from the @SQ header of SAM file'
    command = "grep @SQ " + filename + " | cut -f2 | sed 's/SN\://' > " + dirPath + "/" + name + "_chr.list"
    tt.methods.bashCommand(command)
    geneList = tt.methods.read_list(dirPath + "/" + name + "_chr.list")
    #chromosome lengths from SAM header
    command = "grep @SQ " + filename + " | cut -f4 | sed 's/LN\://' > " + dirPath + "/" + name + "_chr.len"
    tt.methods.bashCommand(command)
    geneLen = tt.methods.read_list(dirPath + "/" + name + "_chr.len")
    
    df_details = pd.DataFrame(data={"chrName" : geneList, "length" : geneLen}).set_index('chrName')
    df_details['length'] = df_details['length'].astype(int)
    
    # tempfiles
    ##list of files if red as chunks
    if chunks > 0:
        chunkList = [i for i in range(0,len(geneList)+1,chunks)]
    else:
        chunkList = [0]
        chunks = len(geneList)+1

    ##genes means chromosomes for this function
    genes = pd.DataFrame(pd.Series(geneList))
    for i in chunkList:
        geneListFileName = "/geneList_" + str(i) + ".tab"
        genes[i:i+chunks].to_csv(dirPath + geneListFileName, sep='\t', header=False,index=False)
        
    ##Select FWD and REV reads from SAM file
    os.chdir(dirPath)
    
    for i in chunkList:
        geneListFileName = "geneList_" + str(i) + ".tab"
        
        print
        #FLAG 0 and 256 for single end reads, aligned, forward match
        command = "grep -v ^@ ../" + filename + " | grep -f " + geneListFileName +\
                  " | awk -F'\t' 'BEGIN{OFS = FS} $2==0||$2==256{print $2,$3, $4, $6, $10, $12}' > "+\
                  name + "_" + str(i) + "_fwd.tab"
        tt.methods.bashCommand(command)
        
        #FLAG 16 and 272 for single end reads, aligned, forward match
        command = "grep -v ^@ ../" + filename + " | grep -f " + geneListFileName +\
                  " | awk -F'\t' 'BEGIN{OFS = FS} $2==16||$2==272{print $2,$3, $4, $6, $10, $12}' > "+\
                  name + "_" + str(i) + "_rev.tab"
        tt.methods.bashCommand(command)
    
    tt.methods.bashCommand("cat "+name+"*fwd.tab > " + name +"_fwd.tab")
    tt.methods.bashCommand("cat "+name+"*rev.tab > " + name +"_rev.tab")

    print("Reads selected.")

    #Reads to profiles
    df_profiles_fwd, df_profiles_rev, log = reads2genome(name=name, dirPath=dirPath,df_details=df_details)
    
    
    #save output
    if pickle==True:
        df_profiles_fwd.to_pickle(path + name + "_PROFILES_reads_fwd.pcl")
        df_profiles_rev.to_pickle(path + name + "_PROFILES_reads_rev.pcl")
    elif pickle==False:
        df_profiles_fwd.to_csv(path + name + "_PROFILES_reads_fwd.csv")
        df_profiles_rev.to_csv(path + name + "_PROFILES_reads_rev.csv")
    #save log
    with open(path + name + "_PROFILES_reads.log", "w") as log_file:
        for row in log:
            log_file.write(str(row) + '\n')

    # clean
    os.chdir(path)
#     shutil.rmtree(name + timestamp)

    print("Done.")