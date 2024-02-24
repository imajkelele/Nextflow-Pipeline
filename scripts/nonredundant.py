from Bio import SeqIO
import argparse
from collections import namedtuple


#parsing input and output as a path
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='ORF fasta files')
    parser.add_argument('--output',  required=True,
                        help='Nonredundant fasta file')

    args = parser.parse_args()
    return args

#searching for nonredundant sequences
def read_nonr_seqs(file_path, output_path):
    seq_hashes = set() #set stores unique values
    
    records = list(SeqIO.parse(file_path, "fasta")) #reading the fasta file
    
    f = open(output_path, 'w') #writing to the file

    for record in records:
        seq_hash = hash(str(record.seq)) #hashing converts sequences into a numerical sequence, optimises sequence comparison
        if seq_hash not in seq_hashes: 
            seq_hashes.add(seq_hash) #adding a new value to the set if it is not already there
            f.write(f'>{record.id}\n{record.seq}\n') #saving to file

    f.close()


if __name__ == '__main__':
    args = parse_args()
    nonr_seqs = read_nonr_seqs(args.input, args.output)
    
