from Bio import SeqIO
import argparse
import hashlib

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='ORF fasta file')
    parser.add_argument('--output', required=True, help='Nonredundant fasta file')
    return parser.parse_args()

def calculate_sequence_hash(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()

def read_nonr_seqs(file_path, output_path):
    seq_hashes = set()
    with open(output_path, 'w') as output_file:
        with open(file_path, 'r') as input_file:
            for record in SeqIO.parse(input_file, "fasta"):
                seq_hash = calculate_sequence_hash(str(record.seq))
                if seq_hash not in seq_hashes:
                    seq_hashes.add(seq_hash)
                    output_file.write(f'>{record.id}\n{record.seq}\n')

if __name__ == '__main__':
    args = parse_args()
    read_nonr_seqs(args.input, args.output)
