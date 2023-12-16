import pandas as pd
from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--input',  required=True,
                    help='ORF fasta files')
parser.add_argument('--output',  required=True,
                    help='Nonredundant fasta file')

args = parser.parse_args()
def fasta_to_dataframe(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))

    data = {'Header': [record.id for record in records],
            'Sequence': [str(record.seq) for record in records]}

    df = pd.DataFrame(data)
    return df
def nonredundant_sequences(df, sequence_column='Sequence'):
    df_unique = df.drop_duplicates(subset=sequence_column)

    return df_unique
def write_df_to_fasta_file(df, output_file):
    with open(output_file, "w") as file:
        for index, row in df.iterrows():
            sequence_name = str(row['Header'])
            sequence = str(row['Sequence'])
            file.write(">" + sequence_name + "\n")
            file.write(sequence[:-1] + "\n")
    return df

write_df_to_fasta_file(nonredundant_sequences(fasta_to_dataframe(args.input)),args.output)