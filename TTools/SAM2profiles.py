import pandas as pd
import numpy as np
import os, collections, shutil, re
import TTools as tt

### support functions ###
def groupCIGAR(cigar_string=""):
    return re.findall(r'(\d+)([A-Z]{1})', cigar_string)


def tostripCIGARfive(match=[]):
    '''Nucleotides without alignment at the 5 end of the read'''
    if "S" in match[0]:
        return int(match[0][0])
    else:
        return None


def tostripCIGARthree(match=[]):
    '''Nucleotides without alignment at the 3 end of the read'''
    if "S" in match[-1]:
        return int(match[-1][0])
    else:
        return None


def stripSubstitutions(match):
    '''Strip substutiotns on both ends'''
    if "S" in match[0]:
        match = match[1:]
    if "S" in match[-1]:
        match = match[:-1]
    return match


def countRead(i=tuple()):
    ''' Takes tuple (position,CIGARstring) and returns list of mapped positions'''
    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    output = np.arange(length, dtype=int)
    return output + position


def countMiddle(i=tuple(), expand=0):
    ''' Takes tuple (position,CIGARstring) and returns list with the middle of mapped read'''
    (position, cigar_string) = i
    match = groupCIGAR(cigar_string)
    # stripping substitutions
    match = stripSubstitutions(match)

    # outputs
    length = sum([int(i) for i, x in match])
    mid = round(length / 2)
    if expand > 0:
        output = np.arange(mid - expand, mid + expand)
        return output + position
    else:
        return np.arange(mid, mid + 1) + position


def countDeletion(i=tuple(), expand=0):
    ''' Takes tuple (position,CIGARstring) and returns list of mapped deletions'''
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
    '''Takes list of tuples position,CIGARstring) and generates profile'''
    list_hits = [countRead(i) for i in l]  # list of lists
    hits = [item for sublist in list_hits for item in sublist]  # transforms to flat list
    # returning the profile
    profile = pd.Series(collections.Counter(hits))  # faster that using zip and numpy
    return profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts


def transcript2profileDeletions(l=[], expand=0, length=int()):
    '''Takes list of tuples position,CIGARstring) and generates profile'''
    list_hits = [countDeletion(i, expand=expand) for i in l]  # list of lists
    hits = [item for sublist in list_hits for item in sublist]  # transforms to flat list
    # returning the profile
    profile = pd.Series(collections.Counter(hits))  # faster that using zip and numpy
    return profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts


def reads2profile(name=str(), dirPath=str(), df_details=pd.DataFrame()):
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


def sam2profiles(filename="", path='', geneList=[], toClear='', df_details=pd.DataFrame(), deletions=False, expand=5):
    # making working directory
    name = filename.replace(".sam", "")
    if toClear:
        name = name.replace(toClear, "")
    timestamp = tt.methods.timestampRandomInt()
    dirPath = path + name + timestamp
    os.makedirs(dirPath)

    # tempfiles
    ##list
    genes = pd.DataFrame(pd.Series(geneList))
    genes.to_csv(dirPath + "/geneList.tab", sep='\t', index=False, header=False)

    ##reads
    os.chdir(dirPath)
    command = "grep -v ^@ ../" + filename + " | grep -f geneList.tab | awk -F'\t' 'BEGIN{OFS = FS} $2==0||$2==256{print $2,$3, $4, $6, $10, $12}' > " + name + ".tab"
    tt.methods.bashCommand(command)
    print("Reads selected.")

    # get profiles
    if deletions == False:
        df_profiles, log = reads2profile(name=name, dirPath=dirPath, df_details=df_details)
        df_profiles.to_csv(path + name + "_PROFILES_reads.csv")
        with open(path + name + "_PROFILES_reads.log", "w") as log_file:
            for row in log:
                log_file.write(str(row) + '\n')

    elif deletions == True:
        output_df, outputDel_df, outputDelExt_df, log = reads2profileDeletions(name=name, dirPath=dirPath,
                                                                               df_details=df_details, expand=expand)
        output_df.to_csv(path + name + "_PROFILES_reads.csv")
        outputDel_df.to_csv(path + name + "_PROFILES_deletions_expand0.csv")
        outputDelExt_df.to_csv(path + name + "_PROFILES_deletions_expand" + str(expand) + ".csv")
        with open(path + name + "_PROFILES_reads.log", "w") as log_file:
            for row in log:
                log_file.write(str(row) + '\n')

    # clean
    os.chdir(path)
    shutil.rmtree(name + timestamp)

    print("Done.")