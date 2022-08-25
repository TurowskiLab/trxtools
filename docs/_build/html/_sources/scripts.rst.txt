Scripts
=======

Main scripts
------------

SAM2profilesGenomic.py
^^^^^^^^^^^^^^^^^^^^^^

usage: Generates genomic profiles from SAM files.

optional arguments:
  -h, --help           show this help message and exit

Options for input files:
  -f FILE              SAM file (default: None)
  -u {read,3end,5end}  Generate profiles using reads (read) or the 3' ends (end) (default: 3end)
  -n                   Save non-coded polyA ends. Can be used ONLY with: -u 3end (default: False)
  -r                   Save all non-coded ends. Can be used ONLY with: -u 3end (default: False)
  --csv                Save output as csv. Default is pickle (default: False)
  -c TOCLEAR           String of signs to be cleared from the name of SAM file (default: _comp_flexbar_STARAligned.out)
  --chunks CHUNKS      Divide list of genes into a chunks of given size (default: 0)


SAM2profiles.py
^^^^^^^^^^^^^^^

usage: Generates profiles from SAM files for a list of given transcripts.

optional arguments:
  -h, --help       show this help message and exit

Options for input files:
  -f FILE          SAM file (default: None)
  -l FILE          list of genes file (default: None)
  -d DETAILS_FILE  file with transcriptome details. (default:
                   /homes/tollervey/COVID19/databases/transcriptome_details_biomart_tRNA_rRNA_UPDATEJan2021.tab)
  --del            Generate additional profiles for deletions (default: False)
  -p               SAve output as pickle DataFrame (default: False)
  -e EXPAND        For deletions position can be expanded by the value on each side (e=5 gives 10 nt long) (default: 5)
  -c TOCLEAR       String of signs to be cleared from the name of SAM file (default: _comp_flexbar_STARAligned.out)
  --chunks CHUNKS  Divide list of genes into a chunks of given size (default: 0)


Supporting scripts
------------------

csv2pickle.py
^^^^^^^^^^^^^
usage: Read CSV file and save as a pickle

optional arguments:
  -h, --help  show this help message and exit

Options for input files:
  -f FILE     SAM file (default: None)
  --gzip      gzip compression (default: False)


mergeSalmon.py
^^^^^^^^^^^^^^
usage: merges Salmon quantification files to one tab file

optional arguments:
  -h, --help  show this help message and exit

Options for input files:
  -i INPUT    String of signs to be found in Salmon output directory (default: None)
  -o OUTPUT   Name of output file (default: merged)
  -u USE      column to be used {TPM,NumReads} (default: NumReads)
  -f FILTER   String of signs to be found in Salmon output directory. Optional as additional filter (default: None)
  -c CLEAR    String of signs to be found in Salmon output directory, will be cleared (default: None)
  -a ADD      String of signs to be added to experiment name (default: None)

