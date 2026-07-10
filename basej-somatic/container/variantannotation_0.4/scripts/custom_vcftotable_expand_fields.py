import pandas as pd
import argparse
import vaex as vx
import os
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed
import dask.dataframe as dd

def process_and_export_chunk(df_chunk, chunk_id, temp_dir):
    ''' Process a chunk of a dataframe and export it to an HDF5 file. '''
    # Process the chunk
    new_columns = ['CSQ_Consequence', 'CSQ_Gene', 'CSQ_Feature_type', 'CSQ_BIOTYPE', 
                   'CSQ_HGVSc', 'CSQ_HGVSp', 'CSQ_cDNA_position', 'CSQ_CDS_position', 
                   'CSQ_Protein_position', 'CSQ_Existing_variation', 'CSQ_SIFT', 
                   'CSQ_PolyPhen', 'CSQ_gnomADe_AF', 'CSQ_gnomADg_AF', 'CSQ_MAX_AF', 
                   'CSQ_MAX_CLIN_SIG']
    indices = [1, 4, 5, 7, 10, 11, 12, 13, 14, 17, 36, 37, 46, 55, 66, 68]

    for index, column_name in zip(indices, new_columns):
        df_chunk[column_name] = df_chunk['CSQ'].apply(lambda x: x.split('|')[index] if len(x.split('|')) > index else None)

    df_chunk = df_chunk.drop(columns=['CSQ'])
    print(f'Chunk {chunk_id} processed')

    # Export the chunk to an HDF5 file
    temp_filename = os.path.join(temp_dir, f'chunk_{chunk_id}.hdf5')
    df_chunk.export_hdf5(temp_filename)
    print(f'Chunk {chunk_id} exported to {temp_filename}')
    return temp_filename

def process_chunks_in_parallel(dataframe, select_cols, chunk_size, temp_dir, num_cpus):
    file_list = []
    chunk_id = 0
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = []
        for df_chunk in vx.read_csv(dataframe, sep='\t', usecols=select_cols, chunk_size=chunk_size):
            futures.append(executor.submit(process_and_export_chunk, df_chunk, chunk_id, temp_dir))
            chunk_id += 1
        for future in concurrent.futures.as_completed(futures):
            file_list.append(future.result())
    return file_list

def get_sample_cols(headers, selected_cols):
    ''' Get sample columns to add to the selected cols '''
    for c in headers:
        if c.endswith('.GT'):
            selected_cols.append(c)
        if c.endswith('.QUAL'):
            selected_cols.append(c)
    return selected_cols

parser = argparse.ArgumentParser(description="Expand fields from vcftotable")
parser.add_argument("-d", "--dataframe", help="path to vcftotable")
parser.add_argument("-o", "--output", help="output filename for the processed vcftotable")
parser.add_argument("-p", "--num_cpus", help="num of processors to use", default=4, type=int)

if __name__ == '__main__':
    args = parser.parse_args()
    
    select_cols = [ "CHROM",
    "ID",
    "REF",
    "ALT",
    "TYPE",
    "CSQ",
    "POS",
    "QUAL" ]

    print('Adding sample columns to selected columns by reading 10 rows')
    df_small = pd.read_csv(args.dataframe, sep='\t', nrows=10)
    select_cols = get_sample_cols(df_small, select_cols)
    
    print('Reading the file in chunk and converting to hdf5')

    temp_dir = 'temp'
    chunk_size = 1000000
    os.makedirs(temp_dir, exist_ok=True)
    num_cpus = args.num_cpus

    file_list = process_chunks_in_parallel(args.dataframe, select_cols, chunk_size, temp_dir, num_cpus)
    print('Parallel processing the chunks to hdf5 completed.')
    
    master_df = vx.open_many(file_list)
    print('Merging the files to a single hdf5 file.')
    master_df.export_hdf5(args.output)


