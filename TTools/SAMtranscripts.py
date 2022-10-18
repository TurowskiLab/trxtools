from trxtools.SAMgeneral import *

####################################################
#   intermediate functions level -2 (from final)   #
####################################################

def transcript2profile(l=[], length=int()):
    '''Takes list of tuples position,CIGARstring) and generates profile. Works with entire reads, for both strands.'''
    list_hits = [countRead(i) for i in l]  # list of lists
    hits = [item for sublist in list_hits for item in sublist]  # transforms to flat list
    # returning the profile
    profile = pd.Series(collections.Counter(hits)).sort_index().astype(float)  # faster that using zip and numpy
    # return profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts
    return profile

def transcript2profileDeletions(l=[], expand=0, length=int()):
    '''Takes list of tuples position,CIGARstring) and generates profile. Not tested for MINUS strand.'''
    list_hits = [countDeletion(i, expand=expand) for i in l]  # list of lists
    hits = [item for sublist in list_hits for item in sublist]  # transforms to flat list
    # returning the profile
    profile = pd.Series(collections.Counter(hits)).sort_index().astype(float)  # faster that using zip and numpy
    # return profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts
    return profile

####################################################
#   intermediate functions level -1 (from final)   #
####################################################

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

####################################################
#              final functions (level 0)           #
####################################################

def sam2profiles(filename="", path='', geneList=[], toClear='', df_details=pd.DataFrame(),
                 deletions=False, expand=5, pickle=False,chunks=0):
    '''Function handling SAM files and generating profiles. Executed using wrapping script SAM2profiles.py.

    :param filename: 
    :type filename: str
    :param path: 
    :type path: str
    :param geneList: list of transcript to be selected
    :type geneList: list
    :param toClear: element of filename to be removed, defaults to ''
    :type toClear: str, optional
    :param df_details: Details of transcripts
    :type df_details: DataFrame
    :param deletions: Generate profile of deletions, defaults to False
    :type deletions: bool, optional
    :param expand: Expand deletions, defaults to 5
    :type expand: int, optional
    :param pickle: save output in pickle format, defaults to False
    :type pickle: bool, optional
    :param chunks: Read SAM file in chunks, defaults to 0
    :type chunks: int, optional

    '''
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