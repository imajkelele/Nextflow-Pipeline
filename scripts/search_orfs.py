#!/usr/bin/env python3
import os
import gzip
from Bio import SeqIO
from lib.orfs import rev_cmpl, translate
import argparse
from collections import namedtuple

#defining codons
starts = set('ATG TTG CTG'.split())
stops = set('TAA TAG TGA'.split())

#parsing input and output as a path
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='Fna.gz File')
    parser.add_argument('--output',  required=True,
                        help='ORFs')

    args = parser.parse_args()
    return args


def read_seqs(file_path):
    NuclSeq = namedtuple('NuclSeq', ['sequence', 'file_id']) #creating a tuple with the possibility to access elements by name
    seq_data = []
    
    handle = gzip.open(file_path, "rt") if file_path.endswith(".fna.gz") else open(file_path) #file extraction
    
    #importing to tuple sequence - sekwencja, ex. file_id - GCA_000276045.1_AKHH01000074.1. 
    for record in SeqIO.parse(handle, "fasta"):
        sequence_id = record.id
        file_and_id = f"{os.path.basename(file_path)[:15]}_{sequence_id}" #adds to file_id, the filename as the last 15 characters in the filename
        seq_data.append( NuclSeq(str(record.seq), file_and_id) )
        
    handle.close()

    return seq_data
    
#seeks ORF in 6 scenarios based on start and stop codons, defined earlier
def find_orfs(seq, minlen):
    for comp, strand in ('comp', seq), ('revcomp', rev_cmpl(seq)):
        for frame in range(3):
            start = None
            for pos in range(frame, len(strand), 3):
                if start is None and strand[pos:pos+3] in starts:
                    start = pos
                elif start is not None and strand[pos:pos+3] in stops:
                    if pos - start + 1 >= minlen: 
                        yield strand[start:pos+3], start, pos+3, comp
                    start = None

#searches for an orf of a given length and stores it in a list consisting of namedtuples
def search_seqs(seq_data, min_orf_length=300):
    ORF = namedtuple('ORF', ['translation', 'file_id', 'start', 'stop', 'comp'])
    orf_data = []

    for sequence, file_id in seq_data:
        for orf, start, stop, comp in find_orfs(sequence, min_orf_length):
            translated_orf = translate(orf)
            if 'X' not in translated_orf: #Checking for a stop codon inside the orf
                orf_data.append( ORF(translated_orf, file_id, start, stop, comp) )

    return orf_data

#Saving the translated orf to a fasta file
def save_orfs(orf_data, file_path):
    with open(file_path, "w") as file:
        for orf in orf_data:
            file.write(f'>{orf.file_id}_{orf.start}_{orf.stop}_{orf.comp}\n') #ex. >GCA_001045625.1_JSAJ01000001.1_7185_8160_comp
            file.write(orf.translation[:-1] + "\n")


if __name__ == '__main__':
    args = parse_args()
    seqs = read_seqs(args.input)
    orfs = search_seqs(seqs)
    save_orfs(orfs, args.output)

