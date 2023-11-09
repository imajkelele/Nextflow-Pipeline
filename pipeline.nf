#!/usr/bin/env nextflow
params.sequences_folder = '/scratch/home-users/idzik/pgps/genomes_test'
params.fasta = '/scratch/home-users/idzik/pgps/output/output_nonredundant.faa'
params.hmm_domains = file('/scratch/home-users/idzik/domains.hmm')
params.hmm = '/scratch/home-users/idzik/pgps/output/hmmer_results.out'
params.pandas = '/scratch/home-users/idzik/pgps/output/output_pandas.csv'
sequences_folder = params.sequences_folder
fasta = params.fasta
hmm_domains = file(params.hmm_domains)
hmm = params.hmm
pandas = params.pandas


process Nonredundant {
    input:
        path sequences_folder  
        path fasta
    output:
        path fasta
    script:
    """
    python3 /scratch/home-users/idzik/po_konsultacji.py --input ${sequences_folder} --output ${fasta}
    """
}
process Search {
    input:
        path fasta 
        path hmm_domains
        path hmm
    output:
        path hmm
    script:
    """
    hmmsearch -E 0.01 --domE  0.01 --incE  0.001  --incdomE   0.001 -o /dev/null  --domtblout ${hmm}  ${hmm_domains} ${fasta}
    """

}

process Hmm_filter {
    input:
        path hmm 
        path pandas
    output:
        path pandas
    script:
    """
    python3 /scratch/home-users/idzik/hmm_to_pandas.py --input ${hmm} --output ${pandas}
    """
}

workflow {
    Hmm_filter(Search(Nonredundant(sequences_folder, fasta),hmm_domains,hmm), pandas)
  //  Nonredundant(sequences_folder)
    //Search(output_fasta, hmm_domains)
    //Hmm_filter(hmm)
}
