sequences, = glob_wildcards('pgps/genomes_test/{sequence}.fna.gz')
rule all:
    input:
        'pgps/output/filter_pandas.csv'
        

rule orfs:
    input:
        'pgps/genomes_test/{sequence}.fna.gz'
    output:
        'pgps/output/orfs/{sequence}_orf.faa'

    log:
        'pgps/log/{sequence}_log_orf.log'
    shell:
        '''python3 /scratch/home-users/idzik/pgps/scripts/create_orf.py \
		--input {input} \
		--output {output} \
		> {log} 2>&1
        '''        
rule merge:
    params:
        mask = rules.orfs.output[0].replace('{sequence}', '*')
    input:
        expand(rules.orfs.output, sequence=sequences)
    output:
        'pgps/output/merged_orfs.faa'
    log:
        'pgps/log/merge.log'
    shell:
        '''
        cat {input} > {output} 2> {log}
        '''

rule nonredundant:
    input:
        rules.merge.output
    output:
        'pgps/output/output_nonredundant.faa'
    log:
        'pgps/log/log_nonredundant.log'
    shell:
        '''python3 /scratch/home-users/idzik/pgps/scripts/nonredundant_orfs.py \
		--input {input} \
		--output {output} \
		> {log} 2>&1
        '''

rule search:
    params:
        E = 0.01,
        domE = 0.01,
        incE = 0.001,
        incdomE = 0.001,
        o = '/dev/null'
    input:
        hmm_domains   = 'pgps/domains.hmm',
        fasta = rules.nonredundant.output
    output:
        'pgps/output/hmm_result.out'
    log:
        'pgps/log/blastn/hmm_result.log'
    shell:
        '''hmmsearch -E    {params.E}                   \
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
    input:
        rules.search.output
    output:
        'pgps/output/filter_pandas.csv'   
    log:
        'pgps/log/blastn/filter.log'
    shell:
        """
        python3 /scratch/home-users/idzik/pgps/scripts/hmm_to_pandas.py \
        --input {input} \
        --output {output}\
        """
