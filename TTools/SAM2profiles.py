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


# def stripCIGAR(match=[]):
#     stripFive = tostripCIGARfive(match)
#     print(stripFive)
#     if stripFive:
#         match = match[1:]

#     stripThree = tostripCIGARthree(match)
#     print(stripThree)
#     if stripThree:
#         match = match[-1:]
#     return match

def countRead(cigar_string="", position=int(), how="read", expand=0):
    '''
    '''
    if not how in ['read', 'middle', 'deletion', 'del']:
        exit("How not recognized")

    match = groupCIGAR(cigar_string)
    # stripping substitutions
    if "S" in match[0]:
        match = match[1:]
    if "S" in match[-1]:
        match = match[:-1]

    # outputs
    length = sum([int(i) for i, x in match])

    if how == 'read':
        output = np.arange(length, dtype=int)
        return output + position

    elif how == 'middle':
        mid = round(length / 2)
        if expand > 0:
            output = np.arange(mid - expand, mid + expand)
            return output + position
        else:
            return np.arange(mid, mid + 1) + position

    elif how in ['deletion', 'del']:
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
def transcript2profile(n="", df=pd.DataFrame(), seq=pd.Series(), how='read', expand=0):
    if not how in ['read', 'middle', 'deletion', 'del']:
        exit("How not recognized")

    output_df = pd.DataFrame()
    output_df['sequence'] = list(seq[0])  # sequence is python-style numbered

    hits = np.array([])
    for i, row in df.iterrows():
        pos = int(row['position'])
        match = groupCIGAR(row['CIGAR'])
        shiftFive = tostripCIGARfive(match)
        length = sum([int(i) for i, x in match])
        read = row['sequence']

        # testing
        #         print(read)
        #         print(output_df[pos-1:pos+length-shiftFive-1])
        #         print(countRead(row['CIGAR'],position=pos))
        #         print(countRead(row['CIGAR'],position=pos,how='del'))
        hits = np.append(hits, countRead(row['CIGAR'], position=pos, how=how, expand=expand))

    # returning the profile
    profile = pd.Series(collections.Counter(hits))  # faster that using zip and numpy
    return profile.reindex(pd.RangeIndex(profile.index.max() + 1)).fillna(0)  # fills spaces with 0 counts


def reads2profile(name=str(), dirPath=str(),
                  df_details=pd.DataFrame(),
                  df_sequences=pd.DataFrame(),
                  how='read', expand=0):
    if not how in ['read', 'middle', 'deletion', 'del']:
        exit("How not recognized")
    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']

    df_input = pd.read_csv(dirPath + "/" + name + ".tab", sep="\t", names=cols)

    output_df = pd.DataFrame()

    for n, df in df_input.groupby('name'):
        details = df_details.loc[n]
        seq = df_sequences.loc[n]

        try:
            s1 = transcript2profile(n=n, df=df, seq=seq, how=how, expand=expand)
        except:
            print(n)

        output_df = output_df.append(s1.rename(n))

    return output_df


def sam2profiles(filename="", path='', geneList=[], toClear='',
                 df_details=pd.DataFrame(), df_sequences=pd.DataFrame(),
                 how='read', expand=0):
    if not how in ['read', 'middle', 'deletion', 'del']:
        exit("How not recognized")

    # making working directory
    name = filename.strip(toClear + ".sam")
    timestamp = tt.methods.timestampRandomInt()
    dirPath = path + name + timestamp
    os.makedirs(dirPath)

    # tempfiles
    ##list
    genes = pd.Series(geneList)
    genes = pd.DataFrame("_" + genes + "_").to_csv(dirPath + "/geneList.tab", sep='\t', index=False, header=False)

    ##reads
    os.chdir(dirPath)
    command = "grep -v ^@ ../" + filename + " | grep -f geneList.tab | awk -F'\t' 'BEGIN{OFS = FS} $2==0||$2==256{print $2,$3, $4, $6, $10, $12}' > " + name + ".tab"
    tt.methods.bashCommand(command)
    print("Reads selected.")

    # get profiles
    df_profiles = reads2profile(name=name, dirPath=dirPath, df_details=df_details, df_sequences=df_sequences, how=how,
                                expand=expand)
    df_profiles.to_csv(path + name + "_PROFILES_" + how + "_expand" + str(expand) + ".csv")

    # clean
    os.chdir(path)
    shutil.rmtree(name + timestamp)

    print("Done.")