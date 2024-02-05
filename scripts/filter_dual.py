import os
import gzip
from Bio import SeqIO
from collections import namedtuple
import pandas as pd
import argparse

#Importuje plik pofiltrowany plik outputowy hmmer, outputem jest lista sekwencji i ścieżka do folderu do którego zapisywane są pliki gff3
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='Filtered Hmmer Output File')
    parser.add_argument('--output_fasta', required=True,
                        help='Fasta file with sequences that contain both Peptidase_M23 and SH3_5 domain')
    parser.add_argument('--output_gff3_folder', required=True,
                        help='Output folder for gff3 files')

    args = parser.parse_args()
    return args

#skrypt bliźniaczy do readseqs z search_orfs.py
def read_seqs(file_path):
    NuclSeq = namedtuple('NuclSeq', ['header', 'sequence', 'file_id'])
    seq_data = []

    with gzip.open(file_path, "rt") if file_path.endswith(".fna.gz") else open(file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence_id = record.id
            file_and_id = f"{os.path.basename(file_path)[:15]}_{sequence_id}"
            seq_data.append(NuclSeq(record.description, str(record.seq), file_and_id))

    return seq_data

#zapisuje do pliku fasta sekwencje których sname matchuje te z filtrowanego wyniku hmmer
def create_final_fasta(df,output_file,genomes_folder="genomes_test"):
    matching_sequences = set()

    for index, row in df.iterrows():
        file = row['aname'] #aname- ex. GCA_001230465.1
        sname = row['sname'] #sname- ex. CTZB01000116.1
        file_path = f'{genomes_folder}/{file}_genomic.fna.gz' #w pliku aname
        sequences = read_seqs(file_path)

        for seq in sequences:
            if sname in seq.header: #szuka sekwencji odpowiadającym sname
                matching_sequences.add(f'>{seq.header}\n{seq.sequence}\n')

    with open(output_file, 'w') as output_handle:
        output_handle.write(''.join(matching_sequences)) #zapisuje je do jednego finalnego pliku fasta


def convert_df_to_gff3(df, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    # Grupowanie DataFrame według unikalnych wartości 'sname'
    grouped_df = df.groupby('sname')

    for sname, group in grouped_df:
        # Tworzenie nazwy pliku GFF3 dla danego 'sname'
        output_file_path = os.path.join(output_folder, f'{sname}.gff3')

        with open(output_file_path, 'w') as output_handle: #to trzeba dostoswać do tego co chcemy mieć w plikach gff3
            output_handle.write("##gff\n")
            for index, row in group.iterrows():
                seq_name = row['sname']
                source = 'HMMER'
                feature_type = row['qname']
                start_pos = row['start']
                end_pos = row['stop']
                comp = row['comp']
                score = row['seqscore']
                strand = '+' if start_pos <= end_pos else '-'

                # Tworzenie linii GFF3 i zapis do pliku
                gff_line = f"{seq_name}\t{source}\t{feature_type}\t{start_pos}\t{end_pos}\t{score}\t{strand}\n"
                output_handle.write(gff_line)

#modyfikacje pofiltrowanego pliku outputowego hmmer
def extract_dual(hmmer_out,output_fasta,output_gff3):
    df = pd.read_csv(hmmer_out)

    df[['GCA', 'aname', 'sname', 'start', 'stop', 'comp']] = df['tname'].str.split('_', expand=True) #podzielenie pierwszwej kolumny żeby wyizolować aname,sname, start,stop i comp

    df['aname'] = df['GCA'] + '_' + df['aname']  # po podzieleniu przez "_" łączenie aname z powrotem do ex. GCA_001230465.1

    df = df.drop(columns=['GCA', 'tname']) #usunięcie nadmiaru kolumn

    column_order = ['aname', 'sname', 'start', 'stop', 'comp'] + [col for col in df.columns if col not in ['aname', 'sname', 'start', 'stop', 'comp']] #ustawienie w dobrej kolejności

    df = df[column_order]

    #wyizolowanie sekwencji które się duplikują i przynajmniej dla jednego rekordu występuje "Peptidase_M23" i przynajmniej dla jednego "SH3_5"
    df_duplicates = df.groupby('sname').filter(lambda group: any(group['qname'] == 'Peptidase_M23') and any(group['qname'] == 'SH3_5')) 

    create_final_fasta(df_duplicates,output_fasta) #tworzenie finalnego fasta
    convert_df_to_gff3(df_duplicates,output_gff3)  #tworzenie plików gff3


if __name__ == '__main__':
    args = parse_args()
    extract_dual(args.input, args.output_fasta, args.output_gff3_folder)