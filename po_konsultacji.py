#!/usr/bin/env python3
import os
import gzip
from Bio import SeqIO
import pandas as pd
from orfs import rev_cmpl, translate
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--input',  required=True,
                    help='Folder with compressed sequences')
parser.add_argument('--output',  required=True,
                    help='Fasta nonredundant sequences')

args = parser.parse_args()


starts = set('ATG TTG CTG'.split())
stops = set('TAA TAG TGA'.split())

def find_orfs(seq, minlen):
    for strand in seq, rev_cmpl(seq):
        if strand in seq:
            comp="comp"
        else:
            comp="revcomp"    
        for frame in range(3):
            start = None
            for pos in range(frame, len(strand), 3):
                if start is None and strand[pos:pos+3] in starts:
                    start = pos
                elif start is not None and strand[pos:pos+3] in stops:
                    if pos - start + 1 >= minlen: 
                        yield strand[start:pos+3], start, pos+3, comp
                    start = None
def read_all_fna_gz_to_dataframe(folder_path):
    data = []
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fna.gz"):
            file_path = os.path.join(folder_path, file_name)
            
            with gzip.open(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequence_id = record.id
                    file_and_id = f"{file_name[:15]}_{sequence_id}"
                    data.append([str(record.seq), file_and_id])

    columns = ['Sequence', 'File_ID']
    df = pd.DataFrame(data, columns=columns)
    
    return df

def df_seq_to_df_orf(df_seq, min_orf_length=300):
    data = []

    for sequence, file_id in zip(df_seq['Sequence'], df_seq['File_ID'].astype(str)):
        for orf, start, stop, comp in find_orfs(sequence, min_orf_length):
            translated_orf = translate(orf)
            if 'X' not in translated_orf:
                data.append([translated_orf, file_id, start, stop, comp])

    columns = ['Translated_ORF', 'File_ID', 'Start', 'Stop', 'Comp']
    df_orfs = pd.DataFrame(data, columns=columns)

    return df_orfs

def nonredundant_sequences(df, sequence_column='Translated_ORF'):
    df_unique = df.drop_duplicates(subset=sequence_column)

    return df_unique
def write_df_to_fasta_file(df, file_name,output_folder = "pgps/output"):
    output_folder = "pgps/output"
    os.makedirs(output_folder, exist_ok=True)  # Create the output directory if it doesn't exist
    file_path = os.path.join(output_folder, file_name)

    with open(file_path, "w") as file:
        for index, row in df.iterrows():
            sequence_name = str(row['File_ID'])
            sequence = str(row['Translated_ORF'])
            start = str(row['Start'])
            stop = str(row['Stop'])
            comp = str(row['Comp'])
            file.write(">" + sequence_name +"_"+start+"_"+stop+"_"+comp+ "\n")
            file.write(sequence[:-1] + "\n")

def write_df_to_fasta_file(df, file_path):
    with open(file_path, "w") as file:
        for index, row in df.iterrows():
            sequence_name = str(row['File_ID'])
            sequence = str(row['Translated_ORF'])
            start = str(row['Start'])
            stop = str(row['Stop'])
            comp = str(row['Comp'])
            file.write(">" + sequence_name + "_" + start + "_" + stop + "_" + comp + "\n")
            file.write(sequence[:-1] + "\n")
    return df


write_df_to_fasta_file(nonredundant_sequences(df_seq_to_df_orf(read_all_fna_gz_to_dataframe(args.input))),args.output)

