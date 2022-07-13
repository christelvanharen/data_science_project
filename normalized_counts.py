"""
Author: Moshtach Ismail
This script converts the raw counts to normalized counts.
"""
import pandas as pd
from sklearn import preprocessing
import sys

# de input bestand is matrix.txt
# de output bestand is genormaliseerde_countsAllsamples.txt

def een(raw_counts, normalised_counts):
    """
    The raw counts file is read in using pandas, so that the
    counts can be normalized with the sklearn.preprocessing.
    """
    # read file into pandas dataframe
    df = pd.read_table(f"{raw_counts}",
                       sep="\t", usecols=range(1, 7),
                       skiprows=1)
    
    x = df.values  # gives numpy array back
    column = ''
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x) # normalized values
    df = pd.DataFrame(x_scaled) # normalized values into pandas dataframe
    # makes new file with new variables 
    df.to_csv(f"{normalised_counts}", sep='\t', mode='w')


def main(raw_counts, normalised_counts):
    een(raw_counts, normalised_counts)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
