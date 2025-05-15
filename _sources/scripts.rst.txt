Scripts
=======

Main scripts
------------

SAM2profilesGenomic.py
^^^^^^^^^^^^^^^^^^^^^^

This script generates genomic profiles from SAM files. It provides various options to customize the generation of profiles based on the input SAM file.

Usage:
    SAM2profilesGenomic.py -f file.sam -u 3end -n -s polyA

Options:
    -f FILE, --sam_file FILE
        SAM file to be processed.

    -u {read,3end,5end,del}, --use {read,3end,5end,del}
        Generate profiles using reads (read) or the 3' ends (3end), 5' ends (5end), or deletions (del).

    -n, --noncoded
        Save non-coded ends. Can be used ONLY with: -u 3end.

    -s {polyA}, --noncoded_suffix {polyA}
        Select non-coded ends. Can be used ONLY with: -u 3end.

    -e EXPAND, --expand EXPAND
        Will expand position of each deletion by +/- the value. Works ONLY with: -u del.

    -c TOCLEAR, --toClear TOCLEAR
        String of signs to be cleared from the name of SAM file.

    --chunks CHUNKS
        Divide list of genes into chunks of given size.

Example:
    SAM2profilesGenomic.py -f file.sam -u 3end -n -s polyA

SAM2profilesTranscripts.py
^^^^^^^^^^^^^^^^^^^^^^^^^^

This script generates profiles from SAM files for a list of given transcripts.
It reads a SAM file, a list of genes, and a file with transcriptome details,
and processes the data to create profiles. The output can be saved as a pickle
DataFrame if specified.

Usage:
    SAM2profilesTranscripts.py -f file.sam -l list_of_genes.txt -d transcriptome_details.tab

Options:
    -f FILE, --sam_file FILE
        SAM file to be processed.
    -l FILE, --list_file FILE
        File containing a list of genes.
    -d FILE, --details_file FILE
        File with transcriptome details. Default is '/homes/tollervey/COVID19/databases/transcriptome_details_biomart_tRNA_rRNA_UPDATEJan2021.tab'.
    --del
        Generate additional profiles for deletions.
    -p, --pickle
        Save output as a pickle DataFrame.
    -e INT, --expand INT
        For deletions, position can be expanded by the value on each side (e=5 gives 10 nt long). Default is 5.
    -c STRING, --toClear STRING
        String of signs to be cleared from the name of SAM file. Default is '_comp_flexbar_STARAligned.out'.
    --chunks INT
        Divide list of genes into chunks of given size. Default is 0.

Example:
    SAM2profilesTranscripts.py -f file.sam -l list_of_genes.txt -d transcriptome_details.tab


genomeNascentFolding.py
^^^^^^^^^^^^^^^^^^^^^^^

This script processes a given FASTA file to generate sliding windows of sequences and their reverse complements. 
It then folds these sequences using RNAfold and parses the folding output to calculate the free energy (dG) of the folded sequences.
The results are saved in WIG and BigWig formats for visualization.

Usage:
    python genomeNascentFolding.py -f genome.fasta -w 65 -t 37 -s both

Arguments:

    -f, --file
      Path to the input FASTA file (required).

    -g, --gzip
      Enable gzip compression for the input file.

    -w, --window
      Size of the sliding window (default: 65).

    -t, --temp
      Temperature for RNA folding (default: 37).
    
    -s, --strand
      Strand to process ('plus', 'minus', or 'both'; default: 'both').

Output:
    The script generates sliding window sequences and their reverse complements, folds them using RNAfold, 
    and saves the folding free energy (dG) values in WIG and BigWig formats. The output files are saved in the same 
    directory as the input file, with appropriate naming conventions.



Supporting scripts
------------------
genome2WindowsGTF.py
^^^^^^^^^^^^^^^^^^^^

This script generates a GTF file for a given genome with fixed-size windows (default 100bp).

Usage:
    python genome2WindowsGTF.py -f genome.len -w 1000

Options:
    -f FILE      File with genome length
    -w WINDOW    Window size (default: 100)

Description:
    The script reads a genome length file and generates a GTF file with fixed-size windows.
    Each window is represented as an exon feature in the GTF file. The script creates windows
    for both the positive and negative strands of the genome.

Output:
    The script outputs a GTF file with the specified window size for the given genome.

csv2pickle.py
^^^^^^^^^^^^^
This script reads a CSV file and saves it as a pickle file. It supports optional gzip compression.

Usage:
    python csv2pickle.py -f <filename.csv> --gzip

Arguments:

    -f FILE, --csv_file FILE
        Path to the input CSV file.

    --gzip
        Enable gzip compression for the output pickle file.

Notes:
    The script reads the CSV file from the current working directory and saves the pickle file in the same directory.

mergeSalmon.py
^^^^^^^^^^^^^^

This script merges multiple Salmon quantification files into a single tab-delimited file. It takes the name of the Salmon output directory and the column to be used (TPM or NumReads). The script can also filter the files to be loaded, clear some strings from the experiment name, and add some strings to the experiment name. The output is saved in the current directory with the name 'merged.tab'.

Usage:
    python mergeSalmon.py -i 'Salmon_output' -u TPM -f 'sample' -a 'RNAseq'

Options:
    -i
      String of signs to be found in Salmon output directory

    -o
      Name of output file (default: 'merged')

    -u 
      Column to be used {TPM, NumReads} (default: 'NumReads')

    -f 
      String of signs to be found in Salmon output directory. Optional as additional filter (default: None)

    -c  
      String of signs to be found in Salmon output directory, will be cleared (default: None)

    -a 
      String of signs to be added to experiment name (default: None)


fasta2slidingWindows.py
^^^^^^^^^^^^^^^^^^^^^^^
This script converts a FASTA file to sliding windows as sequences and their reverse complements (RC). 
The output is saved in a folder named with the prefix 'window' followed by the window size. Each chromosome 
is saved in a separate file. If the tab_output option is selected, the output is saved as a tab-separated file; 
otherwise, it is saved as FASTA files.

Usage:
    python fasta2slidingWindows.py -f genome.fasta -w 65 -t 37 -s both --tab

Arguments:

    -f, --file       
        Path to the input FASTA file (required).

    -g, --gzip       
        Enable gzip compression for the input file.

    -w, --window     
        Sliding window size (default: 65).

    -t, --temp       
        Temperature (default: 37).

    -s, --strand     
        Strand to process, options are "plus", "minus", or "both" (default: "both").
        
    --tab            
        Save output as tab-separated file.

polyAfraction.py
^^^^^^^^^^^^^^^^^^^^^^

This script generates fraction (eg. polyA) for genomic profiles from BigWig files. Bot fwd and rev files will be processed.

Usage:
    polyAfraction.py -f file_polyA_fwd.bw

Options:
    -f FILE, --sam_file FILE
        BigWig file to be processed. Requires file with suffix "_fwd.bw" or "_rev.bw". Either can be used.
        The script will automatically detect the suffix and process the file accordingly.
        The script will also look for the corresponding file with the other suffix (if not provided).
        For example, if you provide "file_polyA_fwd.bw", it will look for "file_polyA_rev.bw" in the same directory.
    -s SUFFIX, --suffix SUFFIX
        Suffix for the input file (default: polyA).
    -h, --help
        Show this help message and exit.

Example:
    polyAfraction.py -f file_polyA_fwd.bw