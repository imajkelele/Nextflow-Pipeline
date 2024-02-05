import argparse
import pandas as pd
import re
#importuje plik outputowy z hmmer i zwraca wyniki jako csv
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

    lines = lines[3:-10] #importowanie linii z pominięciem linii opisowych

    data = [re.split(r'\s+', re.sub(r'\s+', ' ', line.strip())) for line in lines] #usunięcie białych znaków i podział

    df = pd.DataFrame(data)

    df['description of target'] = df.iloc[:, 22:30].apply(lambda x: ''.join(map(str, x)), axis=1) #nadmiar kolumn przez podzielenie 

    df = df.drop(df.columns[22:30], axis=1) #usunięcie nadmiaru kolumn

    column_names = ["tname", "tacc", "tlen", "qname", "qacc", "qlen", "E-value", "seqscore", "seqbias",
                "#", "of", "c-Evalue", "i-Evalue", "domscore", "dombias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"] #nadanie nazw kolumn

    df.columns = column_names
    numerical_columns = ["tlen", "qlen", "E-value", "seqscore", "seqbias",
                "#", "of", "c-Evalue", "i-Evalue", "domscore", "dombias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"] #oznaczenie kolumn z danymi numerycznymi

    df[numerical_columns] = df[numerical_columns].apply(pd.to_numeric, errors='coerce')
    df['qcov'] = ((df['hmm_to'] - df['hmm_from'] + 1) / df['qlen']) #wyliczenie qcov
    df['tcovq'] = ((df['env_to'] - df['env_from'] + 1) / df['tlen']) #wyliczenie tcovq

    condition1 = df['i-Evalue'] <= 0.01 #filtrowanie po i-evalue <= 0.01
    condition2 = df['qcov'] >= 0.80 #filtrowanie po qcov > 0.80
    df = df[condition1 & condition2] #zastosowanie filtrowania
    
    df.to_csv(output, index=False) #zapis do pliku


if __name__ == '__main__':
    args = parse_args()
    hmm_to_pandas(args.input, args.output)
    
