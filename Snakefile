sequences, = glob_wildcards('genomes_test/{sequence}.fna.gz')
rule all:
    input:
        'output/results.faa'
        
rule searchorfs:
    conda:
        'envs/biopy.yml'
    input:
        'genomes_test/{sequence}.fna.gz'
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
        cat {input} > {output} 2> {log}
        '''

rule nonredundant:
    conda:
        'envs/biopy.yml'
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
        'envs/pandas.yml'
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
        'envs/pandas.yml'
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
