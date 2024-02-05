from Bio import SeqIO
import argparse
from collections import namedtuple


#parsowanie inputu i outputu jako ścieżki
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='ORF fasta files')
    parser.add_argument('--output',  required=True,
                        help='Nonredundant fasta file')

    args = parser.parse_args()
    return args

#szukanie sekwencji nieredundantnych 
def read_nonr_seqs(file_path, output_path):
    seq_hashes = set() #set przechowuje unikalne wartości 
    
    records = list(SeqIO.parse(file_path, "fasta")) #odczytywanie pliku fasta
    
    f = open(output_path, 'w') #zapisywanie do pliku

    for record in records:
        seq_hash = hash(str(record.seq)) #hashowanie zamienia sekwencje na ciąg numeryczny, optymalizuje porównywanie sekwencji
        if seq_hash not in seq_hashes: 
            seq_hashes.add(seq_hash) #dodanie nowej wartości do setu, jeśli jeszcze tam nie występuje
            f.write(f'>{record.id}\n{record.seq}\n') #zapisanie do pliku

    f.close()


if __name__ == '__main__':
    args = parse_args()
    nonr_seqs = read_nonr_seqs(args.input, args.output)
    
