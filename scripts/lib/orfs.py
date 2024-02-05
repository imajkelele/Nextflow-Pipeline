from collections import defaultdict

# IUPAC DNA ambiguous alphabet. A dictionary.
dna_ambiguous = {
    'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A',
    'Y' : 'R', 'R' : 'Y', 'W' : 'W', 'S' : 'S',
    'K' : 'M', 'M' : 'K', 'D' : 'H', 'V' : 'B',
    'H' : 'D', 'B' : 'V', 'N' : 'N', 'X' : 'X',
    '-': '-'
}

# Bacterial codon table (codon table No. 11).
# A default dictionary, in case a requested key (codon)
# does not exist, returns 'X' (unknown amino acid).
codon_tab = defaultdict(lambda: 'X', {
    'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
    'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
    'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
    'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

    'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
    'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
    'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
    'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

    'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT'     : 'S',
    'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
    'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
    'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

    'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
    'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
    'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
    'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'
})

# Bacterial start and stop codons. Sets.
starts = set( 'ATG TTG CTG'.split() )
stops  = set( 'TAA TAG TGA'.split() )

def rev_cmpl(seq):
    '''Given a sequence returns its reverse complement.'''
    revcmpl = ''.join( dna_ambiguous[nt] for nt in seq[::-1] )
    return revcmpl

def find_orfs(seq, minlen):
    '''Given a seqeunce and an ORF minimal lenght
       searches for ORFs in both strands. Yields
       every ORF as a string (nucleotide sequnce).'''
    for strand in seq, rev_cmpl(seq):
        for frame in range(3):
            start = None
            for pos in range(frame, len(strand), 3):
                if start is None and strand[pos:pos+3] in starts:
                    start = pos
                elif start is not None and strand[pos:pos+3] in stops:
                    if pos - start + 1 >= minlen:
                        yield strand[start:pos]
                    start = None
                    
def translate(seq):
    '''Given a sequence returns its translation
       as a string (protein sequence) using
       the bacterial codon table (11).'''
    trans =  ''.join(
        codon_tab[ seq[pos:pos+3] ]
        for pos in range(0, len(seq), 3)
    )
    return trans