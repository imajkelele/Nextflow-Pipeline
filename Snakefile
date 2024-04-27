sequences, = glob_wildcards('genomes_all/{sequence}.fna.gz')
rule all:
    input:
        'output/tree.pdf'
        
rule searchorfs:
    conda:
        'envs/biopandas.yml'
    input:
        'genomes_all/{sequence}.fna.gz'
    output:
        'output/orfs/{sequence}_orf.faa'

    log:
        'log/orfs/{sequence}_log_orf.log'
    shell:
        '''python3 scripts/search_orfs.py \
		--input {input} \
		--output {output} \
		> {log} 2>&1
        '''
        
rule merge:
    params:
        mask = rules.searchorfs.output[0].replace('{sequence}', '*')
    input:
        expand(rules.searchorfs.output, sequence=sequences)
    output:
        'output/merged_orfs.faa'
    log:
        'log/merge.log'
    shell:
        '''
        for file in {input}; do
        dd if="$file" of="{output}" conv=notrunc oflag=append 2>> "{log}"
        done
        '''

rule nonredundant:
    conda:
        'envs/biopandas.yml'
    input:
        rules.merge.output
    output:
        'output/nonredundant.faa'
    log:
        'log/nonredundant.log'
    shell:
        '''python3 scripts/nonredundant.py \
		--input {input} \
		--output {output} \
		> {log} 2>&1
        '''

rule hmmersearch:
    conda:
        'envs/hmmer.yml'
    params:
        E = 0.01,
        domE = 0.01,
        incE = 0.001,
        incdomE = 0.001,
        o = '/dev/null'
    threads:
        16
    input:
        hmm_domains   = 'input/domains.hmm',
        fasta = rules.nonredundant.output
    output:
        'output/hmm_result.out'
    log:
        'log/hmm_result.log'
    shell:
        '''hmmsearch --cpu {threads}                     \
                   -E    {params.E}                      \
                  --domE  {params.domE}                  \
                  --incE      {params.incE}              \
                  --incdomE   {params.incdomE}           \
                  --noali                                \
                  --notextw                              \
                  --acc                                  \
                  -o {params.o}                         \
                  --domtblout {output}                  \
                  {input.hmm_domains}                   \
                  {input.fasta}                         \
                   > {log} 2>&1
        '''

rule filter:
    conda:
        'envs/biopandas.yml'
    input:
        rules.hmmersearch.output
    output:
        'output/filter_hmmer.csv'   
    log:
        'log/filter.log'
    shell:
        """
        python3 scripts/filter_hmmer.py \
        --input {input} \
        --output {output}
        > {log} 2>&1\
        """

rule results:
    conda:
        'envs/biopandas.yml'
    input:
        nonredundant = rules.nonredundant.output,
        filter_output = rules.filter.output
        
    output:
        fasta = 'output/results.faa',
        gff3 = 'output/results.gff3'
    log:
        'log/results.log'
    shell:
        """
        python3 scripts/filter_dual.py \
        --input_nonredundant {input.nonredundant} \
        --input_filter_hmmer {input.filter_output} \
        --output_fasta {output.fasta} \
        --output_gff3 {output.gff3} \
        > {log} 2>&1\
        """
rule changenames:
    conda:
        'envs/pandas.yml'
    input:
        rules.results.output.fasta
    output:
        'output/result_organism_names.faa'
    log:
        'log/changenames.log'
    shell:
        """
        python3 scripts/organism_names.py \
		--input {input} \
		--output {output} \
		> {log} 2>&1
        """

rule cdhit:
    conda:
        'envs/cdhit.yml'
    input:
        rules.changenames.output
    output:
        'output/cdhit_result.faa'
    log:
        'log/cdhit.log'
    shell:
        """
        cd-hit -i {input} -o {output} -s 0.9 -aS 0.9
        """

#    TODO: cluster protein seqs from rules.results.output.fasta using cd-hit,
#          minimal mutual coverage 0.9, minimal sequence identity 0.9


rule mafft:
    conda:
        'envs/mafft.yml'
    input:
        rules.cdhit.output
    output:
        'output/cdhit_result_aligned.faa'
    log:
        'log/mafft.log'
    shell:
        """
        mafft {input} > {output}
        """

rule fasttree:
    conda:
        'envs/fasttree.yml'
    input:
        rules.mafft.output
    output:
        'output/fasttree_out.tree'
    log:
        'log/fasttree.log'
    shell:
        """
        FastTree {input} > {output}
        """
        
rule figtree:
    conda:
        'envs/figtree.yml'
    input:
        rules.fasttree.output
    output:
        'output/tree.pdf'
    log:
        'log/figtree.log'
    shell:
        """
        figtree -graphic PDF {input} {output}
        """
#rule raxml:
#    conda:
#        'envs/raxml.yml'
#    input:
#        rules.mafft.output
#    output:
#        'output/raxml'
#    params:
#        path = "/scratch/home-users/idzik/Snakemake_MB"
#    log:
#        'log/raxml.log'
#    shell:
#        """
#        raxmlHPC-PTHREADS-SSE3 -T 6 -p 12345 -s {input} -w {params.path} -m PROTCATJTT -n tree
#        """

#    TODO: build a tree for representative sequences given by cd-hit using
#          raxmlHPC-PTHREADS-SSE3

#rule dist:

#    TODO: calculate % species representation across clusteres given by cd-hit,
#          e.g. Cluster 0: 90% S. aureus, 10% S. epidermidis, etc.

#rule tree:
#    TODO: visualise the tree computed by raxmlHPC-PTHREADS-SSE3, include the
#          information on the number of sequences and the top species
#          in a given cluster, e.g. Cluster 0 (13), S. aureus (90%)

