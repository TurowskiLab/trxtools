from TTools.SAMgeneral import *
from TTools.SAMtranscripts import transcript2profile

####################################################
#   intermediate functions level -2 (from final)   #
####################################################

def chromosome2profile3end(l=[], length=int(), strand='FWD'):
    '''Takes list of tuples (position,CIGARstring,sequence) and generates profile'''
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
            readEnd = position+readLength
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


def reads2genome3end(name=str(), dirPath=str(), df_details=pd.DataFrame(),noncoded=True):
    cols = ['score', 'name', 'position', 'CIGAR', 'sequence', 'details']
    df_input_fwd = pd.read_csv(dirPath + "/" + name + "_fwd.tab", sep="\t", names=cols)
    df_input_rev = pd.read_csv(dirPath + "/" + name + "_rev.tab", sep="\t", names=cols)

    output_df_fwd = pd.DataFrame()
    output_df_rev = pd.DataFrame()
    noncoded_fwd = {}
    noncoded_rev = {}
    log = []
    #strand "+"
    for n, df in df_input_fwd.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='FWD')
            output_df_fwd = output_df_fwd.append(s1_profile.rename(n))
            if noncoded==True:
                noncoded_fwd[n] = l1_noncoded

            log.append(n + " - FWD profile and/or list of noncoded ends generated successfully")
        except:
            log.append(n + " - FWD profile and/or list of noncoded ends FAILED")
    
    #strand "-"
    for n, df in df_input_rev.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='REV')
            output_df_rev = output_df_rev.append(s1_profile.rename(n))
            if noncoded==True:
                noncoded_rev[n] = l1_noncoded

            log.append(n + " - REV profile and/or list of noncoded ends generated successfully")
        except:
            log.append(n + " - REV profile and/or list of noncoded ends FAILED")

    return output_df_fwd, output_df_rev, log, noncoded_fwd, noncoded_rev


####################################################
#              final functions (level 0)           #
####################################################

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
    shutil.rmtree(name + timestamp)

    print("Done.")

def sam2genome3end(filename="", path='', geneList=[], toClear='', df_details=pd.DataFrame(),
                pickle=False,chunks=0,noncoded=True):
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
    df_profiles_fwd, df_profiles_rev, log, noncoded_fwd, noncoded_rev = reads2genome3end(name=name, dirPath=dirPath,df_details=df_details)
    
    noncoded_fwd = selectPolyA(parseNoncoded(noncoded_fwd))
    noncoded_rev = selectPolyA(parseNoncoded(noncoded_rev))
    
    #save output
    if pickle==True:
        df_profiles_fwd.to_pickle(path + name + "_PROFILES_3end_fwd.pcl")
        df_profiles_rev.to_pickle(path + name + "_PROFILES_3end_rev.pcl")
    elif pickle==False:
        df_profiles_fwd.to_csv(path + name + "_PROFILES_3end_fwd.csv")
        df_profiles_rev.to_csv(path + name + "_PROFILES_3end_rev.csv")
    if noncoded==True:
        noncoded_fwd.to_csv(path + name + "_noncoded_3end_fwd.csv")
        noncoded_rev.to_csv(path + name + "_noncoded_3end_rev.csv")
    #save log
    with open(path + name + "_PROFILES_3end.log", "w") as log_file:
        for row in log:
            log_file.write(str(row) + '\n')

    # clean
    os.chdir(path)
    shutil.rmtree(name + timestamp)

    print("Done.")

    return noncoded_rev