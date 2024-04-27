import argparse
from Bio import SeqIO
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='CD-Hit results')
    parser.add_argument('--output',  required=True,
                        help='CD-Hit result with organism names')
    args = parser.parse_args()
    return args

def read_fasta_to_dataframe(file_path):
    sequences = {'SequenceID': [], 'Sequence': []}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences['SequenceID'].append(record.id)
        sequences['Sequence'].append(str(record.seq))
    return pd.DataFrame(sequences)


def change_name(input, output, assembly_summary="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"):
    df = pd.read_csv(assembly_summary, sep='\t', skiprows=1)

    df.rename(columns={'#assembly_accession': "assembly_accession"}, inplace=True)

    fasta_df = read_fasta_to_dataframe(input)

    fasta_df[['Prefix', 'Accession', 'ID1', 'Start', 'End', 'Comp']] = fasta_df['SequenceID'].str.split('_', expand=True)

    fasta_df['Accession'] = fasta_df['Prefix'] + '_' + fasta_df['Accession']

    fasta_df.drop(columns=['SequenceID','Prefix'], inplace=True)

    merged_df = pd.merge(fasta_df, df, left_on='Accession', right_on='assembly_accession', how='left')


    with open(output, 'w') as file:
        for index, row in merged_df.iterrows():
            row['organism_name'] = row['organism_name'].replace(' ', '_') ############## _
            file.write(f">{row['organism_name']}_{row['Accession']}_{row['Start']}\n") 
            file.write(f"{row['Sequence']}\n")

if __name__ == '__main__':
    args = parse_args()

    change_name(args.input, args.output)