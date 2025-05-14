from trxtools.sam.SAMgeneral import *
from trxtools.sam.SAMtranscripts import transcript2profile, transcript2profileDeletions
from trxtools.methods import timestamp
import time, shutil
import pyBigWig

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
    profile = pd.Series(collections.Counter(hits)).sort_index().astype(float)  # faster that using zip and numpy
    # profile = profile.reindex(pd.RangeIndex(length + 1)).fillna(0)  # fills spaces with 0 counts
    return profile, noncoded

####################################################
#   intermediate functions level -1 (from final)   #
####################################################

def reads2genome(name=str(), dirPath=str(), df_details=pd.DataFrame(),use="read",logName=str()):
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

    #open log file
    log_file = open(logName, "a")
    log_file.write("Processing SAM to profiles: -u reads" + '\n')
    temp_paths = {} #stores paths for each temp file

    # strand "+" 
    for n, df in df_input_fwd.groupby('name'):
        n_name = n+"_"+use
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1_profile = transcript2profile(l, length=length)
            to_save = pd.DataFrame(s1_profile.rename(n))
            fileName = dirPath+"/temp_"+n_name+"_fwd.pcl.gz"
            to_save.to_pickle(path=fileName,compression='gzip')
            
            temp_paths[n_name+"_fwd"] = fileName
            log_file.write(timestamp()+"\t"+n_name + " - FWD profile generated successfully" + '\n')
        except:
            temp_paths[n_name+"_fwd"] = None
            log_file.write(timestamp()+"\t"+n_name + " - FWD profile FAILED" + '\n')
            
    # strand "-"    
    for n, df in df_input_rev.groupby('name'):
        n_name = n+"_"+use
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1_profile = transcript2profile(l, length=length)
            to_save = pd.DataFrame(s1_profile.rename(n))
            fileName = dirPath+"/temp_"+n_name+"_rev.pcl.gz"
            to_save.to_pickle(path=fileName,compression='gzip')
            
            temp_paths[n_name+"_rev"] = fileName
            log_file.write(timestamp()+"\t"+n_name+ " - REV profile generated successfully" + '\n')
        except:
            
            temp_paths[n_name+"_rev"] = None
            log_file.write(timestamp()+"\t"+n_name+" - REV profile FAILED" + '\n')
    
    log_file.close()
    return temp_paths

def reads2genomeDeletions(name=str(), dirPath=str(), df_details=pd.DataFrame(),use="del", expand=0, logName=str()):
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

    #open log file
    log_file = open(logName, "a")
    log_file.write("Processing SAM to profiles: -u reads" + '\n')
    temp_paths = {} #stores paths for each temp file

    # strand "+" 
    for n, df in df_input_fwd.groupby('name'):
        n_name = n+"_"+use+str(expand)
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1_profile = transcript2profileDeletions(l,expand=expand, length=length)
            to_save = pd.DataFrame(s1_profile.rename(n))
            fileName = dirPath+"/temp_"+n_name+"_fwd.pcl.gz"
            to_save.to_pickle(path=fileName,compression='gzip')
            
            temp_paths[n_name+"_fwd"] = fileName
            log_file.write(timestamp()+"\t"+n_name + " - FWD profile generated successfully" + '\n')
        except:
            temp_paths[n_name+"_fwd"] = None
            log_file.write(timestamp()+"\t"+n_name + " - FWD profile FAILED" + '\n')
            
    # strand "-"    
    for n, df in df_input_rev.groupby('name'):
        n_name = n+"_"+use+str(expand)
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR']))
            s1_profile = transcript2profileDeletions(l,expand=expand, length=length)
            to_save = pd.DataFrame(s1_profile.rename(n))
            fileName = dirPath+"/temp_"+n_name+"_rev.pcl.gz"
            to_save.to_pickle(path=fileName,compression='gzip')
            
            temp_paths[n_name+"_rev"] = fileName
            log_file.write(timestamp()+"\t"+n_name+ " - REV profile generated successfully" + '\n')
        except:
            
            temp_paths[n_name+"_rev"] = None
            log_file.write(timestamp()+"\t"+n_name+" - REV profile FAILED" + '\n')
    
    log_file.close()
    return temp_paths

def reads2genome3end(name=str(), dirPath=str(), df_details=pd.DataFrame(),use="3end",
                noncoded=True,ends='polyA',logName=str(),minLen=3):
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
    #name is reference name (chr name)
    df_input_fwd = pd.read_csv(dirPath + "/" + name + "_fwd.tab", sep="\t", names=cols)
    df_input_rev = pd.read_csv(dirPath + "/" + name + "_rev.tab", sep="\t", names=cols)

    #open log file
    log_file = open(logName, "a")
    log_file.write(timestamp()+"\t"+"Processing SAM to profiles: -u "+use+'\n')
    temp_paths = {} #stores paths for each temp file
    
    # strand "+"
    for n, df in df_input_fwd.groupby('name'):
        n_name = n+"_"+use
        # try:
        length = df_details.loc[n]['length']
        l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
        #processing reads
        s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='FWD')
        to_save = pd.DataFrame(s1_profile.rename(n))
        fileName = dirPath+"/temp_"+n_name+"_fwd.pcl.gz"
        to_save.to_pickle(path=fileName,compression='gzip')
        
        temp_paths[n_name+"_fwd"] = fileName
        log_file.write(timestamp()+"\t"+n_name + " - FWD profile generated successfully" + '\n')

        if noncoded==True:
            try:
                noncoded_profile = selectNoncodedAndProfile(l=l1_noncoded,minLen=3,tail="AAA",letter="A",content=0.75)
                
                # l1_noncoded = parseNoncodedList(l1_noncoded, minLen=minLen)
                # #above file could be saved as raw noncoded ends
                # noncoded_profile = noncoded2profile1(selectEnds(l1_noncoded,ends=ends),length=df_details.loc[n]['length'])
                to_save = pd.DataFrame(noncoded_profile.rename(n))
                fileName = dirPath+"/temp_"+n_name+"_"+ends+"_fwd.pcl.gz"
                to_save.to_pickle(path=fileName,compression="gzip")
                temp_paths[n_name+"_"+ends+"_fwd"] = fileName
                log_file.write(timestamp()+"\t"+n_name + " - FWD "+ends+" profile generated successfully" + '\n')
            except:
                # temp_paths[n+"_fwd"] = None
                log_file.write(timestamp()+"\t"+n_name + " - FWD profile/noncode profile FAILED" + '\n')

    # strand "-"
    for n, df in df_input_rev.groupby('name'):
        n_name = n+"_"+use
        # try:
        length = df_details.loc[n]['length']
        l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
        #processing reads
        s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='REV')
        to_save = pd.DataFrame(s1_profile.rename(n))
        fileName = dirPath+"/temp_"+n_name+"_rev.pcl.gz"
        to_save.to_pickle(path=fileName,compression='gzip')
        
        temp_paths[n_name+"_rev"] = fileName
        log_file.write(timestamp()+"\t"+n_name + " - REV profile generated successfully" + '\n')

        if noncoded==True:
            try:
                noncoded_profile = selectNoncodedAndProfile(l=l1_noncoded,minLen=3,tail="AAA",letter="A",content=0.75)
                # l1_noncoded = parseNoncodedList(l1_noncoded, minLen=minLen)
                # #above file could be saved as raw noncoded ends
                # noncoded_profile = noncoded2profile1(selectEnds(l1_noncoded,ends=ends),length=df_details.loc[n]['length'])
                to_save = pd.DataFrame(noncoded_profile.rename(n))
                fileName = dirPath+"/temp_"+n_name+"_"+ends+"_rev.pcl.gz"
                to_save.to_pickle(path=fileName,compression="gzip")
                temp_paths[n_name+"_"+ends+"_rev"] = fileName
                log_file.write(timestamp()+"\t"+n_name + " - REV "+ends+" profile generated successfully" + '\n')
            except:
                # temp_paths[n+"_fwd"] = None
                log_file.write(timestamp()+"\t"+n_name + " - REV profile/noncode profile FAILED" + '\n')

    return temp_paths

def reads2genome5end(name=str(), dirPath=str(), df_details=pd.DataFrame(),use="5end",logName=str()):
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

    #open log file
    log_file = open(logName, "a")
    log_file.write(timestamp()+"\t"+"Processing SAM to profiles: -u 5end" + '\n')
    temp_paths = {} #stores paths for each temp file

    # strand "+"
    for n, df in df_input_fwd.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='REV') #shortcut: swapping FWD and REV strands
            to_save = pd.DataFrame(s1_profile.rename(n))
            fileName = dirPath+"/temp_"+n+"_fwd.pcl.gz"
            to_save.to_pickle(path=fileName,compression='gzip')
            
            n = n+"_"+use
            temp_paths[n+"_fwd"] = fileName
            log_file.write(timestamp()+"\t"+n + " - FWD profile generated successfully" + '\n')
        except:
            n = n+"_"+use
            temp_paths[n+"_fwd"] = None
            log_file.write(timestamp()+"\t"+n + " - FWD profile FAILED" + '\n')
    
    # strand "-"
    for n, df in df_input_rev.groupby('name'):
        try:
            length = df_details.loc[n]['length']
            l = list(zip(df['position'].astype(int), df['CIGAR'], df['sequence']))         
            #processing reads
            s1_profile, l1_noncoded = chromosome2profile3end(l, length=length, strand='FWD') #shortcut: swapping FWD and REV strands
            to_save = pd.DataFrame(s1_profile.rename(n))
            fileName = dirPath+"/temp_"+n+"_rev.pcl.gz"
            to_save.to_pickle(path=fileName,compression='gzip')
            
            n = n+"_"+use
            temp_paths[n+"_rev"] = fileName
            log_file.write(timestamp()+"\t"+n + " - REV profile generated successfully" + '\n')
        except:
            n = n+"_"+use
            temp_paths[n+"_rev"] = None
            log_file.write(timestamp()+"\t"+n + " - REV profile FAILED" + '\n')

    return temp_paths

# def parseHeader(filename,name,dirPath):
#     ### old version of the function
#     '''Parses the header of a SAM file to extract chromosome names and lengths.

#     :param filename: Path to the SAM file
#     :type filename: str
#     :param name: Name of the experiment
#     :type name: str
#     :param dirPath: Directory path to save intermediate files
#     :type dirPath: str

#     :return: DataFrame with chromosome names and lengths
#     :rtype: pd.DataFrame

#     :example:

#     >>> parseHeader("example.sam", "experiment1", "/path/to/dir")
#     chrName  length
#     chr1     248956422
#     chr2     242193529
#     ...
#     '''
#     #geneList = chromosome list from the @SQ header of SAM file'
#     command = "grep @SQ " + filename + " | cut -f2 | sed 's/SN\://' > " + dirPath + "/" + name + "_chr.list"
#     tt.methods.bashCommand(command)
#     geneList = tt.methods.read_list(dirPath + "/" + name + "_chr.list")
#     #chromosome lengths from SAM header
#     #-f4 for sam as direct STAR output, but -f3 for BAM output and samtools view -h
#     command = "grep @SQ " + filename + " | cut -f3 | sed 's/LN\://' > " + dirPath + "/" + name + "_chr.len"
#     tt.methods.bashCommand(command)
#     geneLen = tt.methods.read_list(dirPath + "/" + name + "_chr.len")
    
#     df_details = pd.DataFrame(data={"chrName" : geneList, "length" : geneLen}).set_index('chrName')
#     df_details['length'] = df_details['length'].astype(int)

#     return df_details

def parseHeader(filename, name, dirPath):
    '''Parses the header of a SAM file to extract chromosome names and lengths.

    :param filename: Path to the SAM file
    :type filename: str
    :param name: Name of the experiment
    :type name: str
    :param dirPath: Directory path to save intermediate files
    :type dirPath: str

    :return: DataFrame with chromosome names and lengths
    :rtype: pd.DataFrame

    :example:

    >>> parseHeader("example.sam", "experiment1", "/path/to/dir")
    chrName  length
    chr1     248956422
    chr2     242193529
    ...
    '''
    geneList = []
    geneLen = []

    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('@'):
                break
            if line.startswith('@SQ'):
                match = re.search(r'SN:(\S+)\s+(?:AS:(\S+)\s+)?LN:(\d+)', line)
                if match:
                    chr_name = match.group(1)
                    chr_len = match.group(3)
                geneList.append(chr_name)
                geneLen.append(int(chr_len))

    df_details = pd.DataFrame(data={"chrName": geneList, "length": geneLen}).set_index('chrName')
    return df_details

####################################################
#              final functions (level 0)           #
####################################################

def sam2genome(filename="", path='', toClear='',chunks=0,use="3end",expand=0,noncoded=True,noncoded_suffix="polyA"):
    '''Function handling SAM files and generating profiles. Executed using wrapping script SAM2profilesGenomic.py.

    :param filename: 
    :type filename: str
    :param path: 
    :type path: str
    :param toClear: element of filename to be removed, defaults to ''
    :type toClear: str, optional
    :param chunks: Read SAM file in chunks, defaults to 0
    :type chunks: int, optional
    :param noncoded_pA: Save non-coded polyA ends, defaults to True
    :type noncoded_pA: bool, optional
    :param noncoded_raw: Save all non-coded ends, defaults to False
    :type noncoded_raw: bool, optional
    '''
    # checking use
    if use not in ["read","3end","5end",'del']: exit("Wrong -u parameter selected")
    if use=="3end" and noncoded==True: 
        if noncoded_suffix not in ["polyA"]: exit("Wrong -s parameter selected")
    
    # making working directory
    name = filename.replace(".sam", "")
    if toClear:
        name = name.replace(toClear, "")
    timestampRnd = tt.methods.timestampRandomInt()
    dirPath = path + name + timestampRnd
    os.makedirs(dirPath)

    #initialize log file
    logName=path + name + "_" + use + ".log"
    log_file = open(logName, "w")
    log_file.write(timestamp()+"\t"+"Start"+"\n")
    log_file.write(timestamp()+"\t"+"Initializing"+"\n")

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
    log_file.write(timestamp()+"\t"+"Selecting reads from SAM files"+"\n")
    os.chdir(dirPath)
    
    for i in chunkList:
        geneListFileName = "geneList_" + str(i) + ".tab"

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

    #Reads to profiles
    log_file.write(timestamp()+"\t"+"Converting reads to profiles"+"\n")
    log_file.close()

    if use=="read":
        temp_paths = reads2genome(name=name, dirPath=dirPath,df_details=df_details, logName=logName)
    elif use=="3end":
        temp_paths = reads2genome3end(name=name, dirPath=dirPath, df_details=df_details,
                                        noncoded=noncoded,ends=noncoded_suffix,logName=logName,minLen=3)
    elif use=="5end":
        temp_paths = reads2genome5end(name=name, dirPath=dirPath,df_details=df_details, logName=logName)
    elif use=="del":
        temp_paths = reads2genomeDeletions(name=name, dirPath=dirPath, df_details=df_details, expand=expand, logName=logName)
    
    
    df_details.to_pickle("details.pcl")

    #save output
    log_file = open(logName, "a")
    log_file.write(timestamp()+"\t"+"Saving output"+"\n")

    chroms = list(df_details['length'].sort_values(ascending=True).to_dict().items()) #sorted for chrom length

    if use=="del":
        use = "del"+str(expand)

    ### fwd strand ###
    log_file.write(timestamp()+"\t"+"Saving FWD strand"+"\n")
    suffix = "_"+use+"_fwd"
    paths = selectSortPaths(paths=temp_paths,chroms=chroms,suffix=suffix)
    bw_name = path + name + "_PROFILE_"+use+"_fwd.bw"
    l = saveBigWig(paths=paths,suffix=suffix,bw_name=bw_name,chroms=chroms)
    log_file.write(timestamp()+"\t"+l+"\n")
    
    ### rev strand ###
    log_file.write(timestamp()+"\t"+"Saving REV strand"+"\n")
    suffix = "_"+use+"_rev"
    paths = selectSortPaths(paths=temp_paths,chroms=chroms,suffix=suffix)
    bw_name = path + name + "_PROFILE_"+use+"_rev.bw"
    l = saveBigWig(paths=paths,suffix=suffix,bw_name=bw_name,chroms=chroms)
    log_file.write(timestamp()+"\t"+l+"\n")

    if use=="3end" and noncoded==True:
        log_file.write(timestamp()+"\t"+"Saving output for noncoded ends"+"\n")

        ### fwd strand ###
        log_file.write(timestamp()+"\t"+"Saving FWD strand (noncoded)"+"\n")
        suffix = "_"+use+"_"+noncoded_suffix+"_fwd"
        paths = selectSortPaths(paths=temp_paths,chroms=chroms,suffix=suffix)
        bw_name = path + name + "_PROFILE_"+use+"_"+noncoded_suffix+"_fwd.bw"
        l = saveBigWig(paths=paths,suffix=suffix,bw_name=bw_name,chroms=chroms)
        log_file.write(timestamp()+"\t"+l+"\n")

        ### rev strand ###
        log_file.write(timestamp()+"\t"+"Saving REV strand (noncoded)"+"\n")
        suffix = "_"+use+"_"+noncoded_suffix+"_rev"
        paths = selectSortPaths(paths=temp_paths,chroms=chroms,suffix=suffix)
        bw_name = path + name + "_PROFILE_"+use+"_"+noncoded_suffix+"_rev.bw"
        l = saveBigWig(paths=paths,suffix=suffix,bw_name=bw_name,chroms=chroms)
        log_file.write(timestamp()+"\t"+l+"\n")

    # clean
    log_file.write(timestamp()+"\t"+"Cleaninig"+"\n")
    os.chdir(path)
    shutil.rmtree(name + timestampRnd, ignore_errors=True)

    log_file.write(timestamp()+"\t"+"Done"+"\n")
    log_file.close()
