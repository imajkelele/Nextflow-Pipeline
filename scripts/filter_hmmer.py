import argparse
import pandas as pd
import re
#imports the output file from hmmer and returns the results as csv
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input',  required=True,
                        help='Hmm file')
    parser.add_argument('--output', required=True,
                        help='Pandas result')

    args = parser.parse_args()
    return args


def hmm_to_pandas(hmm_file,output):
    with open(hmm_file, 'r') as file:
        lines = file.readlines()

    lines = lines[3:-10] #import of lines without description lines

    data = [re.split(r'\s+', re.sub(r'\s+', ' ', line.strip())) for line in lines] #removal of whitespace characters and data splitting

    df = pd.DataFrame(data)

    df['description of target'] = df.iloc[:, 22:30].apply(lambda x: ''.join(map(str, x)), axis=1) #excess columns by dividing

    df = df.drop(df.columns[22:30], axis=1) #removal of excess columns

    column_names = ["tname", "tacc", "tlen", "qname", "qacc", "qlen", "E-value", "seqscore", "seqbias",
                "#", "of", "c-Evalue", "i-Evalue", "domscore", "dombias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"] #naming the columns

    df.columns = column_names
    numerical_columns = ["tlen", "qlen", "E-value", "seqscore", "seqbias",
                "#", "of", "c-Evalue", "i-Evalue", "domscore", "dombias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"] #designation of columns with numerical data

    df[numerical_columns] = df[numerical_columns].apply(pd.to_numeric, errors='coerce')
    df['qcov'] = ((df['hmm_to'] - df['hmm_from'] + 1) / df['qlen']) #qcov calculation
    df['tcovq'] = ((df['env_to'] - df['env_from'] + 1) / df['tlen']) #tcovq calculation

    condition1 = df['i-Evalue'] <= 0.01 #filtering by i-evalue <= 0.01
    condition2 = df['qcov'] >= 0.80 #filtering by qcov > 0.80
    df = df[condition1 & condition2] #applying filters
    
    df.to_csv(output, index=False) #writing to the file


if __name__ == '__main__':
    args = parse_args()
    hmm_to_pandas(args.input, args.output)
    
