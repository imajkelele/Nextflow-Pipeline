import os
import gzip
from Bio import SeqIO
from collections import namedtuple
import pandas as pd
import argparse

#Importuje plik pofiltrowany plik outputowy hmmer, outputem jest lista sekwencji i ścieżka do folderu do którego zapisywane są pliki gff3
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_nonredundant',  required=True,
                        help='Fasta file with nonredundant sequences')
    parser.add_argument('--input_filter_hmmer',  required=True,
                        help='Filtered Hmmer Output File')
    parser.add_argument('--output_fasta', required=True,
                        help='Fasta file with sequences that contain both Peptidase_M23 and SH3_5 domain')
    parser.add_argument('--output_gff3', required=True,
                        help='Output gff3 file')

    args = parser.parse_args()
    return args


def convert_df_to_gff3(df, output_gff3):
    with open(output_gff3, 'w') as output_handle:
        output_handle.write("##gff\n")
        for index, row in df.iterrows():
            seq_name = row['tname']
            source = 'HMMER'
            feature_type = row['qname']
            score = row['seqscore']

            # Modify this line according to the desired GFF3 format
            gff_line = f"{seq_name}\t{source}\t{feature_type}\t{score}\n"
            output_handle.write(gff_line)

import pandas as pd
from Bio import SeqIO

def find_matching_sequences_df(input_fasta, input_csv, output_fasta,output_gff3):
    # Wczytaj dane z pliku CSV do DataFrame
    df = pd.read_csv(input_csv)
    df = df.groupby('tname').filter(lambda group: any(group['qname'] == 'Peptidase_M23') and any(group['qname'] == 'SH3_5')) 

    # Znajdź sekwencje z pliku FASTA, których ID znajduje się w kolumnie 'tname' DataFrame
    matching_sequences = []
    with open(input_fasta, 'r') as fastafile:
        for record in SeqIO.parse(fastafile, 'fasta'):
            seq_id = record.id
            if seq_id in df['tname'].values:
                matching_sequences.append(record)

    # Zapisz znalezione sekwencje do nowego pliku FASTA
    with open(output_fasta, 'w') as outputfile:
        SeqIO.write(matching_sequences, outputfile, 'fasta')
    
    convert_df_to_gff3(df, output_gff3)



if __name__ == '__main__':
    args = parse_args()

    find_matching_sequences_df(args.input_nonredundant,args.input_filter_hmmer, args.output_fasta, args.output_gff3)
