from trxtools.SAMgeneral import *
from trxtools.SAMtranscripts import transcript2profile
import time, shutil

####################################################
#   intermediate functions level -2 (from final)   #
####################################################

def chromosome2profile3end(l=[], length=int(), strand='FWD'):
    '''Generates profile for the 3' ends of reads and saves position of non-coded end

    :param l: list of triple tuples (position, cigar_string, sequence), defaults to []
    :type l: list
    :param length: length of chromosome
    :type length: int
    :param strand: ``'FWD'`` or ``'REV'``, defaults to 'FWD'
    :type strand: str, optional
    :return: profile, noncoded
    :rtype: np.array, list of tuples

    >>> chromosome2profile3end(l=[(10,"3S15M1D9M2S","TTTGCGCAGTCGTGCGGGGCGCAGCGCCC")],length=50,strand="FWD")
    (0     0.0
    1     0.0
    ...
    34    0.0
    35    1.0
    36    0.0
    ...
    50    0.0
    dtype: float64,
    [(35, 'CC')])
    >>> chromosome2profile3end(l=[(40,"3S15M1D9M2S","TTTGCGCAGTCGTGCGGGGCGCAGCGCCC")],length=50,strand="REV")
    (0     0.0
    1     0.0
    ...
    39    0.0
    40    1.0
    41    0.0
    ...
    50    0.0
    dtype: float64,
    [(40, 'AAA')])
    '''

    hits = []
    noncoded = []
    for (position, cigar_string, sequence) in l:
        match = groupCIGAR(cigar_string)
        match = stripCIGAR(match) #removes H

        # outputs
        if strand=="FWD": #strand "+"
            nonMatchThree = tostripCIGARthree(match)
            match = stripSubstitutions(match) #removes S
            readLength = sum([int(i) for i, x in match])
            readEnd = position+readLength-1
            hits.append(readEnd) #appending
            if nonMatchThree: 
                noncoded.append((readEnd,sequence[-nonMatchThree:]))
                
        elif strand=="REV": #strand "-"
            hits.append(position) #appending
            nonMatchFive = tostripCIGARfive(match)
            if nonMatchFive:
                revcomp = tt.methods.reverse_complement_DNA(sequence[:nonMatchFive])
                noncoded.append((position,revcomp))
      
    # returning the profile
    profile = pd.Series(collections.Counter(hits))  # faster that using zip and numpy
    profile = profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts
    return profile, noncoded

####################################################
#   intermediate functions level -1 (from final)   #
####################################################

def reads2genome(name=str(), dirPath=str(), df_details=pd.DataFrame()):
    '''Function used by sam2genome. Works for both strands.

    :param name: name of experiment
    :type name: str
    :param dirPath:
    :type dirPath: str
    :param df_details: lengths of chromosomes
    :type df_details: DataFrame
    :return: output_df_fwd, output_df_rev, log
    :rtype: DataFrame, DataFrame, list
    '''

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
            # output_df_fwd = output_df_fwd.append(s1.rename(n))
            output_df_fwd = pd.concat([output_df_fwd,s1_profile.rename(n)],axis=0, join='outer')

            log.append(n + " - FWD profile generated successfully")
        except:
            log.append(n + " - FWD profile FAILED")
            
    for n, df in df_input_rev.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1 = transcript2profile(l, length=length)
            # output_df_rev = output_df_rev.append(s1.rename(n))
            output_df_rev = pd.concat([output_df_rev,s1_profile.rename(n)],axis=0, join='outer')

            log.append(n + " - REV profile generated successfully")
        except:
            log.append(n + " - REV profile FAILED")

    return output_df_fwd, output_df_rev, log


def reads2genome3end(name=str(), dirPath=str(), df_details=pd.DataFrame(),noncoded=True):
    '''Function used by sam2genome3end. Works for both strands.

    :param name: name of experiment
    :type name: str
    :param dirPath:
    :type dirPath: str
    :param df_details: lengths of chromosomes
    :type df_details: DataFrame
    :param noncoded: If True then will parse and save non-coded ends, defaults to True
    :type noncoded: bool, optional
    :return: output_df_fwd, output_df_rev, log, noncoded_fwd, noncoded_rev
    :rtype: DataFrame, DataFrame, list, DataFrame, DataFrame
    '''

    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']
    df_input_fwd = pd.read_csv(dirPath + "/" + name + "_fwd.tab", sep="\t", names=cols)
    df_input_rev = pd.read_csv(dirPath + "/" + name + "_rev.tab", sep="\t", names=cols)

    output_df_fwd = pd.DataFrame()
    output_df_rev = pd.DataFrame()
    noncoded_fwd = {}
    noncoded_rev = {}
    log = []

    # strand "+"
    for n, df in df_input_fwd.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='FWD')
            # output_df_fwd = output_df_fwd.append(s1_profile.rename(n)) #method is deprecated
            output_df_fwd = pd.concat([output_df_fwd,s1_profile.rename(n)],axis=0, join='outer')
            if noncoded==True:
                noncoded_fwd[n] = l1_noncoded

            log.append(n + " - FWD profile and/or list of noncoded ends generated successfully")
        except:
            log.append(n + " - FWD profile and/or list of noncoded ends FAILED")
    
    # strand "-"
    for n, df in df_input_rev.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='REV')
            # output_df_rev = output_df_rev.append(s1_profile.rename(n)) #method is deprecated
            output_df_rev = pd.concat([output_df_rev,s1_profile.rename(n)],axis=0, join='outer')
            if noncoded==True:
                noncoded_rev[n] = l1_noncoded

            log.append(n + " - REV profile and/or list of noncoded ends generated successfully")
        except:
            log.append(n + " - REV profile and/or list of noncoded ends FAILED")

    return output_df_fwd, output_df_rev, log, noncoded_fwd, noncoded_rev

def reads2genome5end(name=str(), dirPath=str(), df_details=pd.DataFrame()):
    '''Function used by sam2genome5end. Works for both strands.

    :param name: name of experiment
    :type name: str
    :param dirPath:
    :type dirPath: str
    :param df_details: lengths of chromosomes
    :type df_details: DataFrame
    :return: output_df_fwd, output_df_rev, log
    :rtype: DataFrame, DataFrame, list
    '''

    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']
    df_input_fwd = pd.read_csv(dirPath + "/" + name + "_fwd.tab", sep="\t", names=cols)
    df_input_rev = pd.read_csv(dirPath + "/" + name + "_rev.tab", sep="\t", names=cols)

    output_df_fwd = pd.DataFrame()
    output_df_rev = pd.DataFrame()
    log = []

    # strand "+"
    for n, df in df_input_fwd.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='REV') #shortcut: swapping FWD and REV strands
            # output_df_fwd = output_df_fwd.append(s1_profile.rename(n))
            output_df_fwd = pd.concat([output_df_fwd,s1_profile.rename(n)],axis=0, join='outer')

            log.append(n + " - FWD profile and/or list of noncoded ends generated successfully")
        except:
            log.append(n + " - FWD profile and/or list of noncoded ends FAILED")
    
    # strand "-"
    for n, df in df_input_rev.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='FWD') #shortcut: swapping FWD and REV strands
            # output_df_rev = output_df_rev.append(s1_profile.rename(n))
            output_df_rev = pd.concat([output_df_rev,s1_profile.rename(n)],axis=0, join='outer')

            log.append(n + " - REV profile and/or list of noncoded ends generated successfully")
        except:
            log.append(n + " - REV profile and/or list of noncoded ends FAILED")

    return output_df_fwd, output_df_rev, log

def parseHeader(filename,name,dirPath):
    #geneList = chromosome list from the @SQ header of SAM file'
    command = "grep @SQ " + filename + " | cut -f2 | sed 's/SN\://' > " + dirPath + "/" + name + "_chr.list"
    tt.methods.bashCommand(command)
    geneList = tt.methods.read_list(dirPath + "/" + name + "_chr.list")
    #chromosome lengths from SAM header
    #-f4 for sam as direct STAR output, but -f3 for BAM output and samtools view -h
    command = "grep @SQ " + filename + " | cut -f3 | sed 's/LN\://' > " + dirPath + "/" + name + "_chr.len"
    tt.methods.bashCommand(command)
    geneLen = tt.methods.read_list(dirPath + "/" + name + "_chr.len")
    
    df_details = pd.DataFrame(data={"chrName" : geneList, "length" : geneLen}).set_index('chrName')
    df_details['length'] = df_details['length'].astype(int)

    return df_details

####################################################
#              final functions (level 0)           #
####################################################

def sam2genome(filename="", path='', toClear='', pickle=False,chunks=0):
    '''Function handling SAM files and generating profiles. Executed using wrapping script SAM2profilesGenomic.py.

    :param filename: 
    :type filename: str
    :param path: 
    :type path: str
    :param toClear: element of filename to be removed, defaults to ''
    :type toClear: str, optional
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

    df_details = parseHeader(filename=filename,name=name,dirPath=dirPath)
    geneList = df_details.index.tolist()

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
    shutil.rmtree(name + timestamp, ignore_errors=True)

    print("Done.")

def sam2genome3end(filename="", path='', toClear='', pickle=False,chunks=0,noncoded_pA=True,noncoded_raw=False):
    '''Function handling SAM files and generating profiles for the 3' end of reads. Executed using wrapping script SAM2profilesGenomic.py.

    :param filename:
    :type filename: str
    :param path: 
    :type path: str
    :param toClear: element of filename to be removed, defaults to ''
    :type toClear: str, optional
    :param pickle: save output in pickle format, defaults to False
    :type pickle: bool, optional
    :param chunks: Read SAM file in chunks, defaults to 0
    :type chunks: int, optional
    :param noncoded_pA: Save non-coded polyA ends, defaults to True
    :type noncoded_pA: bool, optional
    :param noncoded_raw: Save all non-coded ends, defaults to False
    :type noncoded_raw: bool, optional
    '''

    # making working directory
    name = filename.replace(".sam", "")
    if toClear:
        name = name.replace(toClear, "")
    timestamp = tt.methods.timestampRandomInt()
    dirPath = path + name + timestamp
    os.makedirs(dirPath)

    df_details = parseHeader(filename=filename,name=name,dirPath=dirPath)
    geneList = df_details.index.tolist()

    # tempfiles
    ##list of files if red as chunks
    if chunks > 0:
        chunkList = [i for i in range(0,len(geneList)+1,chunks)]
    else:
        chunkList = [0]
        chunks = len(geneList)+1

    ##genes means chromosomes for this function
    genes = pd.DataFrame(pd.Series(geneList,dtype=str))
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

    noncoded=False
    if noncoded_pA==True: noncoded=True
    elif noncoded_raw==True: noncoded=True

    #Reads to profiles
    df_profiles_fwd, df_profiles_rev, log, noncoded_fwd, noncoded_rev = reads2genome3end(name=name, dirPath=dirPath,df_details=df_details,noncoded=noncoded)

    #parse raw non-coded ends, keep only >=3 nt long
    noncoded_fwd = parseNoncoded(noncoded_fwd, minLen=3)
    noncoded_rev = parseNoncoded(noncoded_rev, minLen=3)

    #saving all non-coded ends
    if noncoded_raw == True:
        noncoded_fwd.to_pickle(path + name + "_noncoded_raw_3end_fwd.pcl")
        noncoded_rev.to_pickle(path + name + "_noncoded_raw_3end_rev.pcl")

    #select only polyA reads (>=3 and 75% of A content) and prepare profile
    try:
        noncoded_profile_fwd = noncoded2profile(selectPolyA(noncoded_fwd),df_details=df_details)
        noncoded_profile_fwd = noncoded_profile_fwd.T
        name_nc_fwd = name
    except:
        print("polyA reads not found for FWD strand")
        name_nc_fwd = name+"_EMPTY"
        noncoded_profile_fwd = pd.DataFrame()

    try:
        noncoded_profile_rev = noncoded2profile(selectPolyA(noncoded_rev),df_details=df_details)
        noncoded_profile_rev = noncoded_profile_rev.T
        name_nc_rev = name
    except:
        print("polyA reads not found for REV strand")
        name_nc_rev = name+"_EMPTY"
        noncoded_profile_rev = pd.DataFrame()

    #save output
    if pickle==True:
        df_profiles_fwd.to_pickle(path + name + "_PROFILES_3end_fwd.pcl")
        df_profiles_rev.to_pickle(path + name + "_PROFILES_3end_rev.pcl")
        if noncoded_pA==True:
            noncoded_profile_fwd.to_pickle(path + name_nc_fwd + "_noncoded_PROFILES_3end_fwd.pcl")
            noncoded_profile_rev.to_pickle(path + name_nc_rev + "_noncoded_PROFILES_3end_rev.pcl")
    
    elif pickle==False:
        df_profiles_fwd.to_csv(path + name + "_PROFILES_3end_fwd.csv")
        df_profiles_rev.to_csv(path + name + "_PROFILES_3end_rev.csv")
        if noncoded_pA==True:
            noncoded_profile_fwd.to_csv(path + name_nc_fwd + "_noncoded_PROFILES_3end_fwd.csv")
            noncoded_profile_rev.to_csv(path + name_nc_rev + "_noncoded_PROFILES_3end_rev.csv")
    
    #save log
    with open(path + name + "_PROFILES_3end.log", "w") as log_file:
        for row in log:
            log_file.write(str(row) + '\n')

    # clean
    os.chdir(path)
    shutil.rmtree(name + timestamp,ignore_errors=True)

    print("Done.")


def sam2genome5end(filename="", path='', toClear='', pickle=False,chunks=0):
    '''Function handling SAM files and generating profiles for the 3' end of reads. Executed using wrapping script SAM2profilesGenomic.py.

    :param filename:
    :type filename: str
    :param path: 
    :type path: str
    :param toClear: element of filename to be removed, defaults to ''
    :type toClear: str, optional
    :param pickle: save output in pickle format, defaults to False
    :type pickle: bool, optional
    :param chunks: Read SAM file in chunks, defaults to 0
    :type chunks: int, optional
    :param noncoded_pA: Save non-coded polyA ends, defaults to True
    :type noncoded_pA: bool, optional
    :param noncoded_raw: Save all non-coded ends, defaults to False
    :type noncoded_raw: bool, optional
    '''

    # making working directory
    name = filename.replace(".sam", "")
    if toClear:
        name = name.replace(toClear, "")
    timestamp = tt.methods.timestampRandomInt()
    dirPath = path + name + timestamp
    os.makedirs(dirPath)

    df_details = parseHeader(filename=filename,name=name,dirPath=dirPath)
    geneList = df_details.index.tolist()
    
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
    df_profiles_fwd, df_profiles_rev, log = reads2genome5end(name=name, dirPath=dirPath,df_details=df_details)

    #save output
    if pickle==True:
        df_profiles_fwd.to_pickle(path + name + "_PROFILES_5end_fwd.pcl")
        df_profiles_rev.to_pickle(path + name + "_PROFILES_5end_rev.pcl")
    
    elif pickle==False:
        df_profiles_fwd.to_csv(path + name + "_PROFILES_5end_fwd.csv")
        df_profiles_rev.to_csv(path + name + "_PROFILES_5end_rev.csv")
    
    #save log
    with open(path + name + "_PROFILES_5end.log", "w") as log_file:
        for row in log:
            log_file.write(str(row) + '\n')

    # clean
    os.chdir(path)
    shutil.rmtree(name + timestamp, ignore_errors=True)

    print("Done.")