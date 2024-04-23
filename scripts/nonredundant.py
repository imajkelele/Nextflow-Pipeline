from Bio import SeqIO
import argparse
from collections import namedtuple

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='ORF fasta files')
    parser.add_argument('--output',  required=True,
                        help='Nonredundant fasta file')

    args = parser.parse_args()
    return args


def read_nonr_seqs(file_path, output_path):
    seq_hashes = set()
    
    records = list(SeqIO.parse(file_path, "fasta"))
    
    f = open(output_path, 'w')

    for record in records:
        seq_hash = hash(str(record.seq))
        if seq_hash not in seq_hashes:
            seq_hashes.add(seq_hash)
            f.write(f'>{record.id}\n{record.seq}\n')

    f.close()


if __name__ == '__main__':
    args = parse_args()
    nonr_seqs = read_nonr_seqs(args.input, args.output)
    
