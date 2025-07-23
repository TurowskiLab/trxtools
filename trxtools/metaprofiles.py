import pandas as pd
import numpy as np
from scipy.sparse import coo_array, lil_array, vstack, save_npz, load_npz
import pyBigWig
import warnings
import trxtools as tt
from pybedtools import BedTool
import pickle
import os


### ALIAS FUNCTIONS
### these will be removed in a later release!

def read_bed(bed_path):
    """Simple BED file parser

    :param bed_path: Path to BED file
    :type bed_path: str
    :return: Contents of BED file in DataFrame form
    :rtype: pandas.DataFrame

    .. warning:: This function is deprecated and will be removed in future versions.
    """
    warnings.warn(
        'read_bed() will be renamed to readBED() in a future release. Update your code to silence this warning.',
        FutureWarning
        )
    return readBED(bed_path)

def matrix_from_bigwig(bw_path, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, align_3end=False):
    '''WARNINIG: This function is replaced by matrixFromBigWig() \n 
    Get matrix with BigWig scores for all regions in a bed df from a single BigWig file.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.
    Flanks are strand-aware.

    :param bw_path: Path to target BigWig file
    :type bw_path: str
    :param bed_df: DataFrame with BED data (use read_bed())
    :type bed_df: pandas.DataFrame
    :param flank_5: length of 5'flank to extend BED regions by, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true replace NaN values with 0 (pybigwig returns positions with 0 coverage as NaN), defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to retrieved scores, defaults to None
    :type pseudocounts: float, optional
    :param align_3end: If true, position 0 in the resulting matrix will be set at the target region's 3'end instead of 5'end, defaults to False
    :type align_3end: bool, optional

    :return: DataFrame with the result score matrix
    :rtype: pandas.DataFrame

    .. warning:: This function is deprecated and will be removed in future versions.
    '''
    warnings.warn(
        'matrix_from_bigwig() will be renamed to matrixFromBigWig() in a future release. Update your code to silence this warning.',
        FutureWarning
    )
    return matrixFromBigWig(
        bw_path=bw_path,
        bed_df=bed_df,
        flank_5=flank_5,
        flank_3=flank_3,
        fill_na=fill_na,
        pseudocounts=pseudocounts,
        align_3end=align_3end
        )

def get_multiple_matrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True, align_3end=False):
    '''WARNINIG: This function is replaced by getMultipleMatrices() \n
    Get score matrices for positions in given regions (with optional flanks) from multiple BigWig files.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.

    :param bw_paths_plus: list of paths to BigWig files (+ strand)
    :type bw_paths_plus: list
    :param bw_paths_minus: list of paths to BigWig files (- strand)
    :type bw_paths_minus: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
    :type flank_5: int, optional
    :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
    :type flank_3: int, optional
    :param fill_na: _description_, defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to datapoints, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
    :type normalize_libsize: bool, optional
    :param align_3end: if True, align profiles at the 3'end of features. Otherwise align at 5'end, defaults to False
    :type align_3end: bool, optional

    :return:  A dictionary containing score matrices for individual BigWig files. Dictionary keys are BigWig file names.
    :rtype: dict

    .. warning:: This function is deprecated and will be removed in future versions.
    '''
    
    warnings.warn(
        'get_multiple_matrices() will be renamed to getMultipleMatrices() in a future release. Update your code to silence this warning.',
        FutureWarning
    )
    return getMultipleMatrices(
        bw_paths_plus=bw_paths_plus,
        bw_paths_minus=bw_paths_minus,
        bed_df=bed_df,
        flank_5=flank_5,
        flank_3=flank_3,
        fill_na=fill_na,
        pseudocounts=pseudocounts,
        normalize_libsize=normalize_libsize,
        align_3end=align_3end
    )


### level -3
def readBED(bed_path):
    '''Simple BED file parser

    :param bed_path: Path to BED file
    :type bed_path: str
    :return: Contents of BED file in DataFrame form
    :rtype: pandas.DataFrame
    '''
    bed_df = pd.read_csv(bed_path, sep='\t', header=None)
    if len(bed_df.columns) > 6:
        bed_df.drop([6,7,8,9,10,11], axis=1, inplace=True)
    # Set score column to 0
    bed_df[4] = 0
    return bed_df

### level -2
def trim_bed(row, align_3end, len_cutoff):
            region_length = row[2] - row[1]
            if align_3end:
                # go len_cutoff bases upstream from 3'end
                if region_length > len_cutoff:
                    if row[5] == '+':
                        row[1] = row[2] - len_cutoff
                    elif row[5] == '-':
                        row[2] = row[1] + len_cutoff
            else:  # go len_cutoff bases downstream from 5'end
                if region_length > len_cutoff:
                    if row[5] == '+':
                        row[2] = row[1] + len_cutoff
                    elif row[5] == '-':
                        row[1] = row[2] - len_cutoff
            return row

def get_chrom_lens(bw):
    '''Get chromosome lengths from a BigWig file.

    :param bw: A pyBigWig object to retrieve chromosome lengths from.
    :type bw: pyBigWig
    :return: A dictionary with chromosome names as keys and their lengths as values.
    :rtype: dict
    '''
    # Get chromosome lengths from a BigWig file.
    chrom_lens = {chrom: bw.chroms()[chrom] for chrom in bw.chroms()}
    return chrom_lens

def get_bw_data(bed_row, bw, flank_5=0, flank_3=0, align_3end=False):
    '''Retrieve BigWig scores for positions in a given region, optionally including flanks of given length.

    :param bed_row: A row from a BED file, containing chromosome, start, end, and strand information.
    :type bed_row: list or pandas.Series
    :param bw: A pyBigWig object to retrieve values from.
    :type bw: pyBigWig
    :param flank_5: Length of the 5' flank to include, defaults to 0.
    :type flank_5: int, optional
    :param flank_3: Length of the 3' flank to include, defaults to 0.
    :type flank_3: int, optional
    :param align_3end: Whether to align the series to the 3' end, defaults to False.
    :type align_3end: bool, optional
    :return: A pandas Series containing the BigWig scores for the specified region.
    :rtype: pandas.Series
    '''

    # Retrieve BigWig scores for positions in a given region, optionally including flanks of given length.
    if bed_row[5] == '+':
        outseries = pd.Series(bw.values(bed_row[0], int(bed_row[1]-flank_5), int(bed_row[2]+flank_3)))
    if bed_row[5] == '-':
        # flip orientation on minus strand
        # keep the index so pandas doesn't flip it again
        outseries = pd.Series(bw.values(bed_row[0], bed_row[1]-flank_3, bed_row[2]+flank_5), dtype='float32')
        rev_index = outseries.index
        outseries = outseries.iloc[::-1]
        outseries.index = rev_index
    if align_3end:
        outseries.index = outseries.index - flank_5 - (bed_row[2] - bed_row[1])
    else:
        outseries.index = outseries.index - flank_5
    return outseries

### level -1
def bed_split_strands(bed_df):
    '''Splits a BED dataframe into two separate dataframes based on strand information.

    :param bed_df: A dataframe containing BED format data. The strand information is expected to be in the 6th column (index 5).
    :type bed_df: pandas.DataFrame
    :return: A tuple containing two dataframes:
        - bed_plus (pandas.DataFrame): Dataframe containing entries with the '+' strand.
        - bed_minus (pandas.DataFrame): Dataframe containing entries with the '-' strand.
    :rtype: tuple
    '''
    # Split dataframe (from read_bw) by strand
    bed_plus = bed_df[bed_df[5] == "+"]
    bed_minus = bed_df[bed_df[5] == "-"]
    return bed_plus, bed_minus

def matrixFromBigWig(
        bw_path, bed_df, flank_5=0, flank_3=0,
        fill_na=True, pseudocounts=None,
        align_3end=False, len_cutoff=None, skip_zeroes=False,
        chunk_baselimit=None, verbose=False):
    '''Get matrix with BigWig scores for all regions in a bed df from a single BigWig file.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.
    Flanks are strand-aware.

    :param bw_path: Path to target BigWig file
    :type bw_path: str
    :param bed_df: DataFrame with BED data (use read_bed())
    :type bed_df: pandas.DataFrame
    :param flank_5: length of 5'flank to extend BED regions by, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true replace NaN values with 0 (pybigwig returns positions with 0 coverage as NaN), defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to retrieved scores, defaults to None
    :type pseudocounts: float, optional
    :param align_3end: If true, position 0 in the resulting matrix will be set at the target region's 3'end instead of 5'end, defaults to False
    :type align_3end: bool, optional
    :return: Dictionary with 'regions' (numpy array), 'matrix' (numpy array), and 'lengths' (numpy array) keys
    :rtype: dict
    '''


    # Sanity checks for len_cutoff
    if len_cutoff is not None:
        if not isinstance(len_cutoff, int):
            raise ValueError("len_cutoff must be an integer")
        else:
            if not align_3end and flank_3 > 0:
                raise ValueError("len_cutoff is not compatible with flank_3 when align_3end=False.")
            if align_3end and flank_5 > 0:
                raise ValueError("len_cutoff is not compatible with flank_5 when align_3end=True.")
            
    if bed_df.empty:
        return {'regions': np.array([]), 'matrix': np.array([[]], dtype=np.float32), 'lengths': np.array([])}

    # Trim regions in BED when using len_cutoff
    if len_cutoff is not None:
        bed_df = bed_df.apply(trim_bed, axis=1, align_3end=align_3end, len_cutoff=len_cutoff)
        # No need to check lengths:
        max_length = flank_5 + flank_3 + len_cutoff
    else:
    # First pass: determine maximum region length for padding
        max_length = 0
        for _, row in bed_df.iterrows():
            region_length = (row[2] - row[1]) + flank_5 + flank_3
            max_length = max(max_length, region_length)
    
    if max_length == 0:
        return {'regions': np.array([]), 'matrix': np.array([[]], dtype=np.float32), 'lengths': np.array([])}

    bw = pyBigWig.open(bw_path)
    if verbose:
        print(f"Opened BigWig file: {bw_path}")
        print(f"Max region length (including flanks): {max_length}")
    # Get chromosome lengths
    chrom_lens = get_chrom_lens(bw)  
    try:
        data_arrays = []
        region_names = []
        region_lengths = []
        curr_chunk_size = 0
        chunk_num = 0
        dumped_data = False
        chunk_files = []
        
        for _, row in bed_df.iterrows():
            region_length = row[2] - row[1]
            if verbose:
                print(f"Current region length: {region_length}", end="")
                print("\r", end="")

            curr_chunk_size += region_length
            
            # Get BigWig values as numpy array directly
            if row[5] == '+':
                left = int(row[1]-flank_5)
                right = int(row[2]+flank_3)
                # Trim region if it extends beyond chromosome length
                if left < 0:
                    left = 0
                if right > chrom_lens[row[0]]:
                    right = chrom_lens[row[0]]
                values = bw.values(row[0], left, right)
            else:  # minus strand
                left = int(row[1]-flank_3)
                right = int(row[2]+flank_5)
                # Trim region if it extends beyond chromosome length
                if left < 0:
                    left = 0
                if right > chrom_lens[row[0]]:
                    right = chrom_lens[row[0]]
                values = bw.values(row[0], left, right)
                values = values[::-1] # reverse for minus strand
            
            # Convert to numpy array with memory-efficient dtype
            arr = np.array(values, dtype=np.float32)
            actual_length = len(arr)
            
            # # Reverrse data for - strand
            # if row[5] == '-':
            #     arr = arr[::-1]
            
            # Handle NaN values
            if fill_na:
                arr = np.nan_to_num(arr, nan=0.0)
            
            # Add pseudocounts
            if pseudocounts is not None:
                arr += pseudocounts
            
            # Skip features with no coverage to reduce mamory load
            if skip_zeroes and np.sum(arr) == 0:
                continue
            
            # Pad array to max_length with zeros
            if actual_length < max_length:
                padded_arr = np.zeros(max_length, dtype=np.float32)
                # Pad on 5'end side if aligning to 3'end
                if align_3end:
                    padded_arr[-actual_length:] = arr
                else:
                    padded_arr[:actual_length] = arr
                arr = padded_arr
            elif actual_length > max_length:
                # Truncate if somehow longer than expected
                arr = arr[:max_length]
                actual_length = max_length
            
            # Store data
            data_arrays.append(arr)
            region_names.append(row[3])
            region_lengths.append(actual_length)
            
            # Check if we need to dump data to disk
            if isinstance(chunk_baselimit, int) and curr_chunk_size > chunk_baselimit:
                if len(data_arrays) > 0:
                    chunk_file = f'temp_chunk_{chunk_num}.npz'
                    # Stack arrays and save with numpy
                    chunk_matrix = np.stack(data_arrays, axis=0)
                    chunk_regions = np.array(region_names)
                    chunk_lengths = np.array(region_lengths)
                    np.savez_compressed(chunk_file, 
                                      matrix=chunk_matrix, 
                                      regions=chunk_regions,
                                      lengths=chunk_lengths)
                    chunk_files.append(chunk_file)
                    
                    # Clear memory
                    data_arrays = []
                    region_names = []
                    region_lengths = []
                    chunk_num += 1
                    dumped_data = True
                    curr_chunk_size = region_length
        
        # Combine all data
    
    finally:
        bw.close()
        if verbose:
            print(f"\nClosed BigWig")
        if dumped_data:
            if verbose:
                print(f"Loading {len(chunk_files)} chunks from disk...")
            # Load and combine all chunks
            all_matrices = []
            all_regions = []
            all_lengths = []
            
            # Load dumped chunks
            for chunk_file in chunk_files:
                chunk_data = np.load(chunk_file)
                all_matrices.append(chunk_data['matrix'])
                all_regions.append(chunk_data['regions'])
                all_lengths.append(chunk_data['lengths'])
            if verbose:
                print("OK")
            
            # Add current data if any
            if verbose:
                print(f"Adding current data...")
            if data_arrays:
                current_matrix = np.stack(data_arrays, axis=0)
                current_regions = np.array(region_names)
                current_lengths = np.array(region_lengths)
                all_matrices.append(current_matrix)
                all_regions.append(current_regions)
                all_lengths.append(current_lengths)
            
            # Combine all
            if verbose:
                print(f"Combining all chunks...")
            if all_matrices:
                final_matrix = np.vstack(all_matrices)
                final_regions = np.concatenate(all_regions)
                final_lengths = np.concatenate(all_lengths)
            else:
                final_matrix = np.array([[]], dtype=np.float32)
                final_regions = np.array([])
                final_lengths = np.array([])
            if verbose:
                print("OK")
        else:
            # No chunking needed
            if verbose:
                print(f"Combining all chunks...")
            if data_arrays:
                final_matrix = np.stack(data_arrays, axis=0)
                final_regions = np.array(region_names)
                final_lengths = np.array(region_lengths)
            else:
                final_matrix = np.array([[]], dtype=np.float32)
                final_regions = np.array([])
                final_lengths = np.array([])
            if verbose:
                print("OK")
        # Clean up temporary files
        if verbose:
            print(f"Cleaning up {len(chunk_files)} temporary files...")
        for chunk_file in chunk_files:
            if os.path.exists(chunk_file):
                os.remove(chunk_file)
        if verbose:
            print("OK")
    if verbose:
            print(f"Finished processing BigWig file: {bw_path}")
    return {'regions': final_regions, 'matrix': final_matrix, 'lengths': final_lengths}

def matrixFromBigWigSparse(
        bw_path, bed_df, flank_5=0, flank_3=0,
        fill_na=True, pseudocounts=None,
        align_3end=False, len_cutoff=None, skip_zeroes=False,
        chunk_baselimit=None, verbose=False):
    '''
    matrixFromBigWig optimized for sparse data (features with mostly 0 coverage).
    Use when your bait protein binds in well defined sites (e.g TFs) as opposed to e.g. processive helicases or polymerases.
    Get matrix with BigWig scores for all regions in a bed df from a single BigWig file.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.
    Flanks are strand-aware.
    

    :param bw_path: Path to target BigWig file
    :type bw_path: str
    :param bed_df: DataFrame with BED data (use read_bed())
    :type bed_df: pandas.DataFrame
    :param flank_5: length of 5'flank to extend BED regions by, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true replace NaN values with 0 (pybigwig returns positions with 0 coverage as NaN), defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to retrieved scores, defaults to None
    :type pseudocounts: float, optional
    :param align_3end: If true, position 0 in the resulting matrix will be set at the target region's 3'end instead of 5'end, defaults to False
    :type align_3end: bool, optional
    :return: Dictionary with 'regions' (numpy array), 'matrix' (numpy array), and 'lengths' (numpy array) keys
    :rtype: dict
    '''

    # Sanity checks for len_cutoff
    if len_cutoff is not None:
        if not isinstance(len_cutoff, int):
            raise ValueError("len_cutoff must be an integer")
        else:
            if not align_3end and flank_3 > 0:
                raise ValueError("len_cutoff is not compatible with flank_3 when align_3end=False.")
            if align_3end and flank_5 > 0:
                raise ValueError("len_cutoff is not compatible with flank_5 when align_3end=True.")
            
    if bed_df.empty:
        return {'regions': np.array([]), 'matrix': np.array([[]], dtype=np.float32), 'lengths': np.array([])}

    # Trim regions in BED when using len_cutoff
    if len_cutoff is not None:
        bed_df = bed_df.apply(trim_bed, axis=1, align_3end=align_3end, len_cutoff=len_cutoff)
        # No need to check lengths:
        max_length = flank_5 + flank_3 + len_cutoff
    else:
    # First pass: determine maximum region length for padding
        max_length = 0
        for _, row in bed_df.iterrows():
            region_length = (row[2] - row[1]) + flank_5 + flank_3
            max_length = max(max_length, region_length)
    
    if max_length == 0:
        return {'regions': np.array([]), 'matrix': np.array([[]], dtype=np.float32), 'lengths': np.array([])}

    bw = pyBigWig.open(bw_path)
    if verbose:
        print(f"Opened BigWig file: {bw_path}")
        print(f"Max region length (including flanks): {max_length}") 
    # Get chromosome lengths
    chrom_lens = get_chrom_lens(bw)       
    try:
        data_arrays = []
        region_names = []
        region_lengths = []
        curr_chunk_size = 0
        chunk_num = 0
        dumped_data = False
        chunk_files = []
        name_files = []
        
        for _, row in bed_df.iterrows():
            region_length = row[2] - row[1]
            if verbose:
                print(f"Current region length: {region_length}", end="")
                print("\r", end="")

            curr_chunk_size += region_length
            
            # Get BigWig values as numpy array directly
            if row[5] == '+':
                left = int(row[1]-flank_5)
                right = int(row[2]+flank_3)
                # Trim region if it extends beyond chromosome length
                if left < 0:
                    left = 0
                if right > chrom_lens[row[0]]:
                    right = chrom_lens[row[0]]
                values = bw.values(row[0], left, right)
            else:  # minus strand
                left = int(row[1]-flank_3)
                right = int(row[2]+flank_5)
                # Trim region if it extends beyond chromosome length
                if left < 0:
                    left = 0
                if right > chrom_lens[row[0]]:
                    right = chrom_lens[row[0]]
                values = bw.values(row[0], left, right)
                values = values[::-1] # reverse for minus strand
            
            # Convert to numpy array with memory-efficient dtype
            arr = np.array(values, dtype=np.float32)
            actual_length = len(arr)
            
            # # Reverrse data for - strand
            # if row[5] == '-':
            #     arr = arr[::-1]
            
            # Handle NaN values
            if fill_na:
                arr = np.nan_to_num(arr, nan=0.0)
            
            # Add pseudocounts
            if pseudocounts is not None:
                arr += pseudocounts
            
            # Skip features with no coverage to reduce mamory load
            if skip_zeroes and np.sum(arr) == 0:
                continue
            
            # Pad array to max_length with zeros
            if actual_length < max_length:
                padded_arr = np.zeros(max_length, dtype=np.float32)
                # Pad on 5'end side if aligning to 3'end
                if align_3end:
                    padded_arr[-actual_length:] = arr
                else:
                    padded_arr[:actual_length] = arr
                arr = padded_arr
            elif actual_length > max_length:
                # Truncate if somehow longer than expected
                arr = arr[:max_length]
                actual_length = max_length
            
            # Store data
            data_arrays.append(arr)
            region_names.append(row[3])
            region_lengths.append(actual_length)
            
            # Check if we need to dump data to disk
            if isinstance(chunk_baselimit, int) and curr_chunk_size > chunk_baselimit:
                if len(data_arrays) > 0:
                    chunk_file = f'temp_chunk_{chunk_num}.npz'
                    name_file = f'temp_chunk_{chunk_num}_names.npz'
                    # Stack arrays and save with numpy
                    chunk_matrix = np.stack(data_arrays, axis=0)
                    chunk_matrix = coo_array(chunk_matrix)
                    chunk_regions = np.array(region_names)
                    # chunk_lengths = np.array(region_lengths)
                    np.savez_compressed(name_file, regions=chunk_regions)
                    save_npz(chunk_file, chunk_matrix)
                    chunk_files.append(chunk_file)
                    name_files.append(name_file)
                    
                    # Clear memory
                    data_arrays = []
                    region_names = []
                    region_lengths = []
                    chunk_num += 1
                    dumped_data = True
                    curr_chunk_size = region_length
        
        # Combine all data
    
    finally:
        bw.close()
        if verbose:
            print(f"\nClosed BigWig")
        if dumped_data:
            if verbose:
                print(f"Loading {len(chunk_files)} chunks from disk...", end="")
            # Load and combine all chunks
            all_matrices = []
            all_regions = []
            all_lengths = []
            
            # Load dumped chunks
            for chunk_file, name_file in zip(chunk_files, name_files):
                chunk_data = lil_array(load_npz(chunk_file))
                chunk_regions_loaded = np.load(name_file, allow_pickle=True)['regions']
                all_matrices.append(chunk_data)
                all_regions.append(chunk_regions_loaded)
                # all_lengths.append(chunk_data['lengths'])
            if verbose:
                print("OK")
            
            # Add current data if any
            if verbose:
                print(f"Adding current data...")
            if data_arrays:
                current_matrix = np.stack(data_arrays, axis=0)
                current_matrix = lil_array(current_matrix)
                current_regions = np.array(region_names)
                # current_lengths = np.array(region_lengths)
                all_matrices.append(current_matrix)
                all_regions.append(current_regions)
                # all_lengths.append(current_lengths)
            
            # Combine all
            if verbose:
                print(f"Combining all chunks...", end="")
            if all_matrices:
                final_matrix = vstack(all_matrices)
                final_regions = np.concatenate(all_regions)
                # final_lengths = np.concatenate(all_lengths)
            else:
                final_matrix = np.array([[]], dtype=np.float32)
                final_regions = np.array([])
                # final_lengths = np.array([])
            if verbose:
                print("OK")
        else:
            # No chunking needed
            if verbose:
                print(f"Combining all chunks...", end="")
            if data_arrays:
                final_matrix = np.stack(data_arrays, axis=0)
                final_matrix = lil_array(final_matrix)
                final_regions = np.array(region_names)
                # final_lengths = np.array(region_lengths)
            else:
                final_matrix = np.array([[]], dtype=np.float32)
                final_regions = np.array([])
                # final_lengths = np.array([])
            if verbose:
                print("OK")
        # Clean up temporary files
        if verbose:
            print(f"Cleaning up {len(chunk_files)} temporary files...", end="")
        for chunk_file in chunk_files:
            if os.path.exists(chunk_file):
                os.remove(chunk_file)
        for name_file in name_files:
            if os.path.exists(name_file):
                os.remove(name_file)
        if verbose:
            print("OK")
    if verbose:
            print(f"Finished processing BigWig file: {bw_path}")
    return {'regions': final_regions, 'matrix': final_matrix}

def join_strand_matrices(plus_dict, minus_dict):
    '''Combine score matrices for regions on + and - strands.

    :param plus_dict: Dictionary containing score matrices for regions on the + strand.
    :type plus_dict: dict
    :param minus_dict: Dictionary containing score matrices for regions on the - strand.
    :type minus_dict: dict
    :raises Exception: If keys in dictionaries do not match.
    :return: Dictionary containing combined score matrices for regions on both strands.
    :rtype: dict
    '''
    # Combine score matrices for regions on + and - strands.

    
    # /home/jm/repos/trxtools/trxtools/metaprofiles.py:143: FutureWarning:
    # The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. 
    # In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. 
    # To retain the old behavior, exclude the relevant entries before the concat operation.
    
    out_dict={}
    for key in list(plus_dict.keys()):
        key_min = key.replace("plus", "minus")
        key_min = key_min.replace("fwd", "rev")
        if key_min not in minus_dict.keys():
            raise Exception("Keys in dictionaries do not match: " + str(key_min))
        else:
            # This raises this FutureWarning:
            # /home/jm/repos/trxtools/trxtools/metaprofiles.py:143: FutureWarning:
            # The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. 
            # In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. 
            # To retain the old behavior, exclude the relevant entries before the concat operation.
            out_dict[key] = pd.concat([plus_dict[key], minus_dict[key_min]], ignore_index=True)
    return out_dict

def peak2matrice(bed_df=pd.DataFrame, peak_file_path='', 
                 g='references/GRCh38.primary_assembly.genome.cleaned.len', 
                 flank_5=0, flank_3=0, fillna=True):
    '''Get score matrices for positions in given regions (with optional flanks) from a single nroadPeak or narrowPeak file. 
    Matrix rows correspond to regions in BED. Columns correpond to nucleotide positions in regions + flanks. Flanks are strand-aware.
    
    :param bed_df: DataFrame with BED data (can use read_bed())
    :type bed_df: pandas.DataFrame
    :param peak_file_path: Path to peak file (broadPeak or narrowPeak format)
    :type peak_file_path: str
    :param g: Path to genome file (chromosome lengths), defaults to 'references/GRCh38.primary_assembly.genome.cleaned.len'
    :type g: str
    :param flank_5: length of 5'flank to extend BED regions by, strand-aware, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, strand-aware, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true replace NaN values with 0 (pybigwig returns positions with 0 coverage as NaN), defaults to True
    :type fill_na: bool, optional
    :return: DataFrame with the result score matrix
    :rtype: pandas.DataFrame
    '''
    
    peak_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pvalue', 'qValue','peak'] #broadPeak do not have 'peak' column
    peak_df = pd.read_csv(peak_file_path, sep='\t', header=None, names=peak_columns)
    a = BedTool.from_dataframe(peak_df)
    
    bed_df.index=bed_df.iloc[:, 3].rename('region')
    
    output_df = pd.DataFrame()
    for i,row in bed_df.iterrows():
        b = BedTool.from_dataframe(pd.DataFrame(row[:7]).T)
        b_extended = b.slop(s=True, l=flank_5, r=flank_3, g=g)
        
        peaks = a.intersect(b_extended).to_dataframe(names=peak_columns)
        
        start,end,strand = b_extended.to_dataframe()[['start','end','strand']].loc[0].tolist()
#         print(end-start)
        profile = pd.Series(index=np.arange(-flank_5,end-start-flank_5),data=0,name=i, dtype=float)

        ###code here how to convert selected beds into Series with the data
       
        for p,row in peaks.iterrows():
            p_start, p_end, p_strand = row['start'], row['end'], row['strand']
            # for "+" strand
            if strand == "+" and (p_strand == "." or p_strand == "+"):
                profile[p_start-start:p_end-start] = row['signalValue']
            
            # for "-" strand
            elif strand == "-" and (p_strand == "." or p_strand == "-"):
                profile[end-p_start:end-p_end] = row['signalValue'] ### !!! ### to be tested
            
        output_df = pd.concat([output_df, profile], axis=1)
        if fillna==True:
            output_df = output_df.fillna(0.0)
    output_df = output_df.T
    output_df.index.name = "region"
    return output_df.reset_index()

### level 0
def getMultipleMatrices(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, 
                        fill_na=True, pseudocounts=None, normalize_libsize=True,
                        align_3end=False, len_cutoff=None, skip_zeroes=False, chunk_baselimit=None,
                        verbose=False):
    
    '''Get score matrices for positions in given regions (with optional flanks) from multiple BigWig files.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.

    :param bw_paths_plus: list of paths to BigWig files (+ strand)
    :type bw_paths_plus: list
    :param bw_paths_minus: list of paths to BigWig files (- strand)
    :type bw_paths_minus: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
    :type flank_5: int, optional
    :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
    :type flank_3: int, optional
    :param fill_na: _description_, defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to datapoints, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
    :type normalize_libsize: bool, optional
    :param align_3end: if True, align profiles at the 3'end of features. Otherwise align at 5'end, defaults to False
    :type align_3end: bool, optional
    :return:  A dictionary containing score matrices for individual BigWig files. Dictionary keys are BigWig file names.
    :rtype: dict
    '''

    warn = False

    # Split strands once to avoid repeated operations
    bed_plus, bed_minus = bed_split_strands(bed_df)
    
    # Initialize output dictionary
    out_dict = {}
    
    # Process BigWig files one at a time to minimize memory usage
    for bw_plus, bw_minus in zip(bw_paths_plus, bw_paths_minus):
        try:
            # Process plus strand
            plus_result = matrixFromBigWig(
                bw_path=bw_plus,
                bed_df=bed_plus,
                flank_5=flank_5,
                flank_3=flank_3,
                fill_na=fill_na,
                pseudocounts=pseudocounts,
                align_3end=align_3end,
                len_cutoff=len_cutoff,
                skip_zeroes=skip_zeroes,
                chunk_baselimit=chunk_baselimit,
                verbose=verbose
                )
            
            # Process minus strand  
            minus_result = matrixFromBigWig(
                bw_path=bw_minus,
                bed_df=bed_minus,
                flank_5=flank_5,
                flank_3=flank_3,
                fill_na=fill_na,
                pseudocounts=pseudocounts,
                align_3end=align_3end,
                len_cutoff=len_cutoff,
                skip_zeroes=skip_zeroes,
                chunk_baselimit=chunk_baselimit,
                verbose=verbose
                )
            if verbose:
                print("Finished processing BigWigs")
            
            # Apply library size normalization if requested
            if normalize_libsize:
                if verbose:
                    print(f"Normalizing by library size...")
                # Calculate library sizes
                with pyBigWig.open(bw_plus) as bwp:
                    libsize_plus = int(bwp.header()['sumData'])
                with pyBigWig.open(bw_minus) as bwm:
                    libsize_minus = int(bwm.header()['sumData'])
                
                libsize_total = libsize_plus + libsize_minus
                
                # Normalize in-place to save memory
                if libsize_total > 0:
                    plus_result['matrix'] = plus_result['matrix'] / libsize_total
                    minus_result['matrix'] = minus_result['matrix'] / libsize_total
                if verbose:
                    print("OK")
            
            # Combine strands and store result
            if verbose:
                print(f"Combining strand matrices")
            if plus_result['regions'].size > 0 and minus_result['regions'].size > 0:
                all_regions = np.concatenate([plus_result['regions'], minus_result['regions']])
                # Pad shorter matrix to make them stack
                shape_plus = np.shape(plus_result['matrix'])
                shape_minus = np.shape(minus_result['matrix'])
                if shape_plus[1] > shape_minus[1]:
                    padded_minus = np.zeros((shape_minus[0], shape_plus[1]), dtype=np.float32)
                    padded_minus[:, :shape_minus[1]] = minus_result['matrix']
                    minus_result['matrix'] = padded_minus
                elif shape_plus[1] < shape_minus[1]:
                    padded_plus = np.zeros((shape_plus[0], shape_minus[1]), dtype=np.float32)
                    padded_plus[:,:shape_plus[1]] = plus_result['matrix']
                    plus_result['matrix'] = padded_plus
                all_matrices = np.vstack([plus_result['matrix'], minus_result['matrix']])
            elif plus_result['regions'].size == 0:
                all_regions = minus_result['regions']
                all_matrices = minus_result['matrix']
            elif minus_result['regions'].size == 0:
                all_regions = plus_result['regions']
                all_matrices = plus_result['matrix']
            
            # all_lengths = np.concatenate([plus_result['lengths'], minus_result['lengths']])
            # out_dict[bw_plus] = {'regions': all_regions, 'matrix': all_matrices, 'lengths': all_lengths}
            # Convert to df
            if verbose:
                print(f"Converting to DataFrame")
            out_dict[bw_plus] = pd.DataFrame(all_matrices, index=all_regions)
            # Clean up intermediate variables to free memory
            if verbose:
                print(f"Cleaning up intermediate variables")
            del plus_result, minus_result, all_regions, all_matrices
            # shift columns so position 0 aligns with to 5' or 3'ends
            if verbose:
                print("Aligning flanks...")
            if align_3end:
                # Because how matrixFromBigWig works, the 3'end is always at the end of the matrix
                # Then we just need to subtract max gene length and add 3'flank to set 0 at the correct position
                out_dict[bw_plus].columns = out_dict[bw_plus].columns - max(out_dict[bw_plus].columns) + flank_3
            else:
                out_dict[bw_plus].columns = out_dict[bw_plus].columns - flank_5
                       
        except Exception as e:
            print(f"Error processing {bw_plus} and {bw_minus}: {e}")
            warn = True
            continue
    if verbose:
        if warn:
            print(f"Execution finished with some errors, examine output for details")
        else:
            print(f"Finished!")
    return out_dict

def getMultipleMatricesSparse(bw_paths_plus, bw_paths_minus, bed_df, flank_5=0, flank_3=0, 
                        fill_na=True, pseudocounts=None, normalize_libsize=True,
                        align_3end=False, len_cutoff=None, skip_zeroes=False, chunk_baselimit=None,
                        verbose=False):
    
    '''
    getMultipleMatrices optimized for sparse data (features with mostly 0 coverage).
    Use when your bait protein binds in well defined sites (e.g TFs) as opposed to e.g. processive helicases or polymerases.
    Get score matrices for positions in given regions (with optional flanks) from multiple BigWig files.
    Matrix rows correspond to regions in BED.
    Columns correpond to nucleotide positions in regions + flanks.

    :param bw_paths_plus: list of paths to BigWig files (+ strand)
    :type bw_paths_plus: list
    :param bw_paths_minus: list of paths to BigWig files (- strand)
    :type bw_paths_minus: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
    :type flank_5: int, optional
    :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
    :type flank_3: int, optional
    :param fill_na: _description_, defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to datapoints, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
    :type normalize_libsize: bool, optional
    :param align_3end: if True, align profiles at the 3'end of features. Otherwise align at 5'end, defaults to False
    :type align_3end: bool, optional
    :return:  A dictionary containing score matrices for individual BigWig files. Dictionary keys are BigWig file names.
    :rtype: dict
    '''

    # Split strands once to avoid repeated operations
    bed_plus, bed_minus = bed_split_strands(bed_df)
    
    # Initialize output dictionary
    out_dict = {}
    
    # Process BigWig files one at a time to minimize memory usage
    for bw_plus, bw_minus in zip(bw_paths_plus, bw_paths_minus):
        try:
            # Process plus strand
            plus_result = matrixFromBigWigSparse(
                bw_path=bw_plus,
                bed_df=bed_plus,
                flank_5=flank_5,
                flank_3=flank_3,
                fill_na=fill_na,
                pseudocounts=pseudocounts,
                align_3end=align_3end,
                len_cutoff=len_cutoff,
                skip_zeroes=skip_zeroes,
                chunk_baselimit=chunk_baselimit,
                verbose=verbose
                )
            
            # Process minus strand  
            minus_result = matrixFromBigWigSparse(
                bw_path=bw_minus,
                bed_df=bed_minus,
                flank_5=flank_5,
                flank_3=flank_3,
                fill_na=fill_na,
                pseudocounts=pseudocounts,
                align_3end=align_3end,
                len_cutoff=len_cutoff,
                skip_zeroes=skip_zeroes,
                chunk_baselimit=chunk_baselimit,
                verbose=verbose
                )
            if verbose:
                print("Finished processing BigWigs")
            
            # Apply library size normalization if requested
            if normalize_libsize:
                if verbose:
                    print(f"Normalizing by library size...", end="")
                # Calculate library sizes
                with pyBigWig.open(bw_plus) as bwp:
                    libsize_plus = int(bwp.header()['sumData'])
                with pyBigWig.open(bw_minus) as bwm:
                    libsize_minus = int(bwm.header()['sumData'])
                
                libsize_total = libsize_plus + libsize_minus
                
                # Normalize in-place to save memory
                if libsize_total > 0:
                    plus_result['matrix'] = plus_result['matrix'] / libsize_total
                    minus_result['matrix'] = minus_result['matrix'] / libsize_total
                if verbose:
                    print("OK")
            
            # Combine strands and store result
            if verbose:
                print(f"Combining strand matrices")
            if plus_result['regions'].size > 0 and minus_result['regions'].size > 0:
                all_regions = np.concatenate([plus_result['regions'], minus_result['regions']])
                # Pad shorter matrix to make them stack
                shape_plus = np.shape(plus_result['matrix'])
                shape_minus = np.shape(minus_result['matrix'])
                if shape_plus[1] > shape_minus[1]:
                    padded_minus = lil_array(np.zeros((shape_minus[0], shape_plus[1]), dtype=np.float32))
                    padded_minus[:, :shape_minus[1]] = minus_result['matrix']
                    minus_result['matrix'] = coo_array(padded_minus)
                elif shape_plus[1] < shape_minus[1]:
                    padded_plus = lil_array(np.zeros((shape_plus[0], shape_minus[1]), dtype=np.float32))
                    padded_plus[:,:shape_plus[1]] = plus_result['matrix']
                    plus_result['matrix'] = coo_array(padded_plus)
                all_matrices = vstack([plus_result['matrix'], minus_result['matrix']])
            elif plus_result['regions'].size == 0:
                all_regions = minus_result['regions']
                all_matrices = minus_result['matrix']
            elif minus_result['regions'].size == 0:
                all_regions = plus_result['regions']
                all_matrices = plus_result['matrix']
            # all_lengths = np.concatenate([plus_result['lengths'], minus_result['lengths']])
            # out_dict[bw_plus] = {'regions': all_regions, 'matrix': all_matrices, 'lengths': all_lengths}
            # Convert to sparse df
            if verbose:
                print(f"Converting to DataFrame")
            out_dict[bw_plus] = pd.DataFrame.sparse.from_spmatrix(all_matrices, index=all_regions)
            # Clean up intermediate variables to free memory
            if verbose:
                print(f"Cleaning up intermediate variables")
            del plus_result, minus_result, all_regions, all_matrices
            # shift columns so position 0 aligns with to 5' or 3'ends
            if verbose:
                print("Aligning flanks...")
            if align_3end:
                out_dict[bw_plus].columns = out_dict[bw_plus].columns - max(out_dict[bw_plus].columns) + flank_3
            else:
                out_dict[bw_plus].columns = out_dict[bw_plus].columns - flank_5
                       
        except Exception as e:
            print(f"Error processing {bw_plus} and {bw_minus}: {e}")
            continue
    if verbose:
        print(f"Finished!")
    return out_dict

def normalizeToLibrary(metaprofile=pd.DataFrame(), bigwig_plus=[], bigwig_minus=[]):
    '''Normalize metaprofile to library size (sum of scores in a bigwig file).
    :param metaprofile: DataFrame with metaprofile data
    :type metaprofile: pandas.DataFrame
    :param bigwig_plus: list of paths to BigWig files for + strand
    :type bigwig_plus: list
    :param bigwig_minus: list of paths to BigWig files for - strand
    :type bigwig_minus: list
    :return: normalized metaprofile DataFrame
    :rtype: pandas.DataFrame
    '''
    out_df = metaprofile.copy()
    if len(bigwig_plus) != len(bigwig_minus):
        raise Exception("Number of + and - strand BigWig files must match")
    
    for bw_plus, bw_minus in zip(bigwig_plus, bigwig_minus):
        with pyBigWig.open(bw_plus) as bwp:
            lib_size_plus = int(bwp.header()['sumData'])
            bwp.close()
        with pyBigWig.open(bw_minus) as bwm:
            lib_size_minus = int(bwm.header()['sumData'])
            bwm.close()
        out_df[bw_plus] = out_df[bw_plus] / (lib_size_plus+lib_size_minus)
    return out_df
    


def getMultipleMatricesFromPeak(peak_paths=[], bed_df=pd.DataFrame,
                 g='references/GRCh38.primary_assembly.genome.cleaned.len', 
                 flank_5=0, flank_3=0):
    '''Get score matrices for positions in given regions (with optional flanks) from multiple peak files.
    Matrix rows correspond to regions in BED. Columns correpond to nucleotide positions in regions + flanks.
    :param peak_paths: list of paths to peak files (broadPeak or narrowPeak format)
    :type peak_paths: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param g: Path to genome file (chromosome lengths), defaults to 'references/GRCh38.primary_assembly.genome.cleaned.len'
    :type g: str
    :param flank_5: length of 5'flank to extend BED regions by, defaults to 0
    :type flank_5: int, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to 0
    :type flank_3: int, optional
    :return: dictionary of score matrices for each peak file
    :rtype: dict
    '''

    out_dict = {}
    for path in peak_paths:
        name = path.split("/")[-1].replace('_peaks.narrowPeak','').replace('_peaks.broadPeak','')
        out_dict[name] = peak2matrice(bed_df=bed_df, peak_file_path=path, 
                 g=g, flank_5=flank_5, flank_3=flank_3)
    
    return out_dict

def metaprofile(matrix_dict, agg_type='mean', normalize_internal=False, subset=None, name_col=None):
    '''Calculate metaprofiles from score matrices by aggregating each position in all regions.

    :param matrix_dict: Dict containing score matrices returned by get_multiple_matrices()
    :type matrix_dict: dict
    :param agg_type: Type of aggregation to use. Available: 'mean', 'median', 'sum', defaults to 'mean'
    :type agg_type: str, optional
    :param normalize: if true, normalize each profile internally (i.e. gene-wise) (x/sum(x)) before aggregating, defaults to False
    :type normalize: bool, optional
    :raises Exception: _description_
    :return: dataframe containing metaprofile values for each position in each of the input matrices (i.e. bigwig files)
    :rtype: pandas.DataFrame
    '''

    if agg_type not in ['mean', 'median', 'sum']:
        raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")

    if isinstance(subset, pd.DataFrame):
        if name_col is not None:
            matrix_dict = {key: value[value[name_col].isin(subset.index)] for key, value in matrix_dict.items()}
        else:
            matrix_dict = {key: value[value.index.isin(subset.index)] for key, value in matrix_dict.items()}
    elif isinstance(subset, list):
        if name_col is not None:
            matrix_dict = {key: value[value[name_col].isin(subset)] for key, value in matrix_dict.items()}
        else:
            matrix_dict = {key: value[value.index.isin(subset)] for key, value in matrix_dict.items()}
        # matrix_dict = {key: value[value['region'].isin(subset)] for key, value in matrix_dict.items()}


    if normalize_internal:
        if name_col is not None:
            dropped = {key: value.set_index(name_col).div(value.sum(axis=1,numeric_only=True),axis=0) for key, value in matrix_dict.items()}
        else:
            dropped = {key: value.div(value.sum(axis=1,numeric_only=True),axis=0) for key, value in matrix_dict.items()}
        return pd.DataFrame({key: value.agg(agg_type,numeric_only=True) for key, value in dropped.items()})

    else:
        return pd.DataFrame({key: value.agg(agg_type, numeric_only=True) for key, value in matrix_dict.items()})

def metaprofile_new(matrix_dict, agg_type='mean', normalize_internal=False, subset=None):
    '''Calculate metaprofiles from score matrices by aggregating each position in all regions.

    :param matrix_dict: Dict containing score matrices returned by get_multiple_matrices()
    :type matrix_dict: dict
    :param agg_type: Type of aggregation to use. Available: 'mean', 'median', 'sum', defaults to 'mean'
    :type agg_type: str, optional
    :param normalize_internal: if true, normalize each profile internally (i.e. gene-wise) (x/sum(x)) before aggregating, defaults to False
    :type normalize_internal: bool, optional
    :param subset: DataFrame or list containing the subset of regions to select, defaults to None
    :type subset: Union[pd.DataFrame, list], optional
    :raises Exception: _description_
    :return: dataframe containing metaprofile values for each position in each of the input matrices (i.e. bigwig files)
    :rtype: pandas.DataFrame
    '''

    if agg_type not in ['mean', 'median', 'sum']:
        raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")

    # Apply subset filtering if provided
    if isinstance(subset, pd.DataFrame):
        subset_regions = set(subset.index)
        matrix_dict = {key: {'regions': value['regions'][np.isin(value['regions'], list(subset_regions))],
                            'matrix': value['matrix'][np.isin(value['regions'], list(subset_regions))],
                            'lengths': value['lengths'][np.isin(value['regions'], list(subset_regions))]}
                      for key, value in matrix_dict.items()}
    elif isinstance(subset, list):
        subset_regions = set(subset)
        matrix_dict = {key: {'regions': value['regions'][np.isin(value['regions'], subset)],
                            'matrix': value['matrix'][np.isin(value['regions'], subset)],
                            'lengths': value['lengths'][np.isin(value['regions'], subset)]}
                      for key, value in matrix_dict.items()}

    result_dict = {}
    for key, value in matrix_dict.items():
        matrix = value['matrix']
        lengths = value['lengths']
        
        if normalize_internal:
            # Normalize each row by its sum, using actual lengths to avoid padding zeros
            normalized_matrix = np.zeros_like(matrix)
            for i, length in enumerate(lengths):
                row = matrix[i, :length]  # Only consider actual data, not padding
                row_sum = np.sum(row)
                if row_sum > 0:
                    normalized_matrix[i, :length] = row / row_sum
            matrix = normalized_matrix
        
        # Aggregate across all regions (rows) for each position (column)
        if agg_type == 'mean':
            # For mean, we need to consider only non-padded positions
            result = np.zeros(matrix.shape[1])
            counts = np.zeros(matrix.shape[1])
            for i, length in enumerate(lengths):
                result[:length] += matrix[i, :length]
                counts[:length] += 1
            # Avoid division by zero
            result = np.divide(result, counts, out=np.zeros_like(result), where=counts!=0)
        elif agg_type == 'median':
            # For median, create a masked array considering actual lengths
            masked_data = []
            max_len = matrix.shape[1]
            for pos in range(max_len):
                valid_values = []
                for i, length in enumerate(lengths):
                    if pos < length:
                        valid_values.append(matrix[i, pos])
                if valid_values:
                    masked_data.append(np.median(valid_values))
                else:
                    masked_data.append(0.0)
            result = np.array(masked_data)
        elif agg_type == 'sum':
            # For sum, we can sum all values but should consider only actual data
            result = np.zeros(matrix.shape[1])
            for i, length in enumerate(lengths):
                result[:length] += matrix[i, :length]
        
        result_dict[key] = result
    
    return pd.DataFrame(result_dict)

def regionScore(bw_paths_plus, bw_paths_minus, bed_df, agg_type='sum', flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True):
    '''Calculate coverage or statistic for multiple regions in multiple BigWig files.

    :param bw_paths_plus: list of paths to BigWig files (+ strand)
    :type bw_paths_plus: list
    :param bw_paths_minus: list of paths to BigWig files (- strand)
    :type bw_paths_minus: list
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions
    :type bed_df: pandas.DataFrame
    :param agg_type: operation to perform on region scores. Available options: 'sum', 'mean', 'median', defaults to 'sum'
    :type agg_type: str, optional
    :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
    :type flank_5: int, optional
    :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
    :type flank_3: int, optional
    :param fill_na: If true, NaNs will be replaced with zeroes (recommended, as pyBigWig reports positions with 0 score as NaN), defaults to True
    :type fill_na: bool, optional
    :param pseudocounts: pseudocounts to add to datapoints, defaults to None
    :type pseudocounts: float, optional
    :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
    :type normalize_libsize: bool, optional
    :return:  DataFrame with calculated scores. Rows are regions/genes, columns are BigWig files.
    :rtype: pandas.DataFrame
    '''

    if agg_type not in ['mean', 'median', 'sum']:
        raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")
    # Get score matrices
    outdict = getMultipleMatrices(
        bw_paths_plus=bw_paths_plus,
        bw_paths_minus=bw_paths_minus,
        bed_df=bed_df,
        flank_5=flank_5,
        flank_3=flank_3,
        fill_na=fill_na,
        pseudocounts=pseudocounts,
        normalize_libsize=normalize_libsize,
        align_3end=False
        )
    # Aggregate all positions per gene using chosen function
    outdict = {key:value.agg(func=agg_type, axis=1) for key, value in outdict.items()}
    out_df = pd.DataFrame(outdict)
    return out_df


## OLD VERSION
## kept for debugging and backwards compatibility
## remove in future
# def regionScore(bw_paths_plus, bw_paths_minus, bed_df, agg_type='sum', flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True):
#     '''Calculate coverage or statistic for multiple regions in multiple BigWig files.

#     :param bw_paths_plus: list of paths to BigWig files (+ strand)
#     :type bw_paths_plus: list
#     :param bw_paths_minus: list of paths to BigWig files (- strand)
#     :type bw_paths_minus: list
#     :param bed_df: dataframe in BED format containing genomic coordinates of target regions
#     :type bed_df: pandas.DataFrame
#     :param agg_type: operation to perform on region scores. Available options: 'sum', 'mean', 'median', defaults to 'sum'
#     :type agg_type: str, optional
#     :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
#     :type flank_5: int, optional
#     :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
#     :type flank_3: int, optional
#     :param fill_na: If true, NaNs will be replaced with zeroes (recommended, as pyBigWig reports positions with 0 score as NaN), defaults to True
#     :type fill_na: bool, optional
#     :param pseudocounts: pseudocounts to add to datapoints, defaults to None
#     :type pseudocounts: float, optional
#     :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
#     :type normalize_libsize: bool, optional
#     :return:  DataFrame with calculated scores. Rows are regions/genes, columns are BigWig files.
#     :rtype: pandas.DataFrame
#     '''

#     if agg_type not in ['mean', 'median', 'sum']:
#         raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")
#     # Get score matrices
#     outdict = getMultipleMatrices(
#         bw_paths_plus=bw_paths_plus,
#         bw_paths_minus=bw_paths_minus,
#         bed_df=bed_df,
#         flank_5=flank_5,
#         flank_3=flank_3,
#         fill_na=fill_na,
#         pseudocounts=pseudocounts,
#         normalize_libsize=normalize_libsize,
#         align_3end=False
#         )
#     # Aggregate all positions per gene using chosen function
#     result_dict = {}
#     for key, value in outdict.items():
#         regions = value['regions']
#         matrix = value['matrix']
#         lengths = value['lengths']
        
#         # Calculate aggregated score for each region
#         region_scores = []
#         for i, length in enumerate(lengths):
#             row_data = matrix[i, :length]  # Only consider actual data, not padding
#             if agg_type == 'sum':
#                 score = np.sum(row_data)
#             elif agg_type == 'mean':
#                 score = np.mean(row_data)
#             elif agg_type == 'median':
#                 score = np.median(row_data)
#             region_scores.append(score)
        
#         result_dict[key] = pd.Series(region_scores, index=regions)
    
#     out_df = pd.DataFrame(result_dict)
#     return out_df

### level 1
def binMultipleMatrices(mm={}, bins=[50, 10, 50], bed_df=pd.DataFrame(),
                        flank_5=None, flank_3=None, region_col=None):
    '''Bin multiple matrices of tRNA profiles into a single dataframe
      
    :param mm: dictionary of matrices of gene profiles 
    :type mm: dict
    :param bins: list of three integers, defaults to [50, 10, 50]
    :type bins: list, optional
    :param bed_df: dataframe in BED format containing genomic coordinates of target regions, defaults to pd.DataFrame()
    :type bed_df: pandas.DataFrame, optional
    :param flank_5: length of 5'flank to extend BED regions by, defaults to None
    :type flank_5: int, optional
    :type bed_df: pandas.DataFrame, optional
    :param flank_3: length of 3'flank to extend BED regions by, defaults to None
    :type flank_3: int, optional
    :type bed_df: pandas.DataFrame, optional
    :raises ValueError: if bins is not a list of three integers
    :raises ValueError: if input_value is not an integer or a DataFrame
    :return: dictionary of binned matrices of tRNA profiles
    :rtype: dict
    '''

    if len(bins) != 3:
        raise ValueError(f"{bins} takes three numbers")
        
    def process_input(input_value):
        if isinstance(input_value, pd.DataFrame):
            return tt.methods.bed2len(input_value)
        elif isinstance(input_value, int):
            return pd.Series(index=bed_df.iloc[:, 3].rename('region'), data=input_value) #uses region as index
        else:
            raise ValueError(f"{input_value} should be either an integer or a DataFrame")

    # Concatenate all the results into length_df
    length_df = pd.concat([process_input(flank_5),
                           tt.methods.bed2len(bed_df),
                           process_input(flank_3)], axis=1)

    results_mm = {} 
    # Iterate over the datasets (keys) of the mm dictionary
    for bw_key, value in mm.items():
        results_df = pd.DataFrame(columns=bed_df.iloc[:, 3].rename('region'))

        # bin individual profiles of tRNAs
        if region_col is not None:
            for region,row in value.set_index(region_col).iterrows():
                results_df[region] = tt.profiles.binCollect3(s1=row,
                                            lengths=length_df.loc[region].tolist(),
                                            bins=bins)
        elif 'region' in value.columns:
            for region,row in value.set_index('region').iterrows():
                results_df[region] = tt.profiles.binCollect3(s1=row,
                                            lengths=length_df.loc[region].tolist(),
                                            bins=bins)
        else:
            for region,row in value.iterrows():
                results_df[region] = tt.profiles.binCollect3(s1=row,
                                            lengths=length_df.loc[region].tolist(),
                                            bins=bins)
        
        results_df.index = np.arange(-bins[0],bins[1]+bins[2])
        results_mm[bw_key] = results_df.T.reset_index()
    
    return results_mm

## helper functions
def selectSubsetMM(matrix_dict, subset=None):
    '''Select a subset of regions from the matrix dictionary.

    :param matrix_dict: Dictionary containing score matrices.
    :type matrix_dict: dict
    :param subset: DataFrame or list containing the subset of regions to select, defaults to None
    :type subset: Union[pd.DataFrame, list], optional
    :return: Dictionary containing the selected subset of score matrices.
    :rtype: dict
    '''
    if isinstance(subset, pd.DataFrame):
        subset_regions = set(subset.index)
        return {key: {'regions': value['regions'][np.isin(value['regions'], list(subset_regions))],
                     'matrix': value['matrix'][np.isin(value['regions'], list(subset_regions))],
                     'lengths': value['lengths'][np.isin(value['regions'], list(subset_regions))]}
               for key, value in matrix_dict.items()}
    elif isinstance(subset, list):
        subset_regions = set(subset)
        return {key: {'regions': value['regions'][np.isin(value['regions'], subset)],
                     'matrix': value['matrix'][np.isin(value['regions'], subset)],
                     'lengths': value['lengths'][np.isin(value['regions'], subset)]}
               for key, value in matrix_dict.items()}
    else:
        return matrix_dict

# def regionScore2(bw_paths_plus, bw_paths_minus, bed_df, agg_type='sum', flank_5=0, flank_3=0, fill_na=True, pseudocounts=None, normalize_libsize=True):
#     """
#     Calculate coverage or statistic for multiple regions in multiple BigWig files.

#     :param bw_paths_plus: list of paths to BigWig files (+ strand)
#     :type bw_paths_plus: list
#     :param bw_paths_minus: list of paths to BigWig files (- strand)
#     :type bw_paths_minus: list
#     :param bed_df: dataframe in BED format containing genomic coordinates of target regions
#     :type bed_df: pandas.DataFrame
#     :param agg_type: operation to perform on region scores. Available options: 'sum', 'mean', 'median', defaults to 'sum'
#     :type agg_type: str, optional
#     :param flank_5: number of nt that input regions will be extended by on 5' side, defaults to 0
#     :type flank_5: int, optional
#     :param flank_3: number of nt that input regions will be extended by on 3' side, defaults to 0
#     :type flank_3: int, optional
#     :param fill_na: If true, NaNs will be replaced with zeroes (recommended, as pyBigWig reports positions with 0 score as NaN), defaults to True
#     :type fill_na: bool, optional
#     :param pseudocounts: pseudocounts to add to datapoints, defaults to None
#     :type pseudocounts: float, optional
#     :param normalize_libsize: normalization to library size (sum of scores in a bigwig file), defaults to True
#     :type normalize_libsize: bool, optional
#     :return:  DataFrame with calculated scores. Rows are regions/genes, columns are BigWig files.
#     :rtype: pandas.DataFrame
#     """    
#     if agg_type not in ['mean', 'median', 'sum']:
#         raise Exception("Wrong agg_type; available values: 'mean', 'median', 'sum'")
    
#     bed_plus, bed_minus = bed_split_strands(bed_df)
    
#     out_df_list = []
#     # Get scores from BigWig
#     for bigwig in zip(bw_paths_plus, bw_paths_minus):
#         out_plus = matrixFromBigWig(
#             bigwig[0],
#             bed_plus,
#             flank_5=flank_5,
#             flank_3=flank_3,
#             fill_na=fill_na,
#             pseudocounts=pseudocounts,
#             normalize_libsize=normalize_libsize
#             ).set_index('region')
#         out_minus = matrixFromBigWig(
#             bigwig[0],
#             bed_minus,
#             flank_5=flank_5,
#             flank_3=flank_3,
#             fill_na=fill_na,
#             pseudocounts=pseudocounts,
#             normalize_libsize=normalize_libsize
#             ).set_index('region')
#         # Merge strands and aggregate
#         out_merged = pd.concat([out_plus, out_minus]).agg(func=agg_type, axis=1)
#         out_merged.columns = [bigwig[0]]
#         print(out_merged)
#         out_df_list.append(out_merged)
#     # Join df's as columns
#     return pd.concat(out_df_list, axis=1)