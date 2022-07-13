"""
Author: Moshtach Ismail
This script converts the raw counts to normalized counts.
"""

import pandas as pd
from sklearn import preprocessing
import sys

# de input bestand is matrix.txt
# de output bestand is

def test(raw_counts):
    """
    the first value of each row will be put in a list.
    This is so it can be added to the file after
    normalizing the values.

    :param raw_counts: this is the path to the file with the raw counts
    :return: c: list with the first value of each row (Gene names)
    """
    c = [] # gene names
    # reads the file with the raw counts
    df = pd.read_table(raw_counts, sep='\t',
                       error_bad_lines=False)

    q = df.iloc[:,0] # only gets first value of each row
    for x in q:
        c.append(x)

    return c

def een(raw_counts,c):
    """
    The raw counts file is read in using pandas, so that the
    counts can be normalized with the sklearn.preprocessing.
    """
    df = pd.read_table(raw_counts,
                       sep="\t", usecols=range(1, 7), # only gets
                       # columns with values
                       skiprows=1) # skips header

    x = df.values  # gives numpy array back
    column = ''
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)  # normalized values
    df = pd.DataFrame(x_scaled)  # normalized values into pandas
    # dataframe
    # add the the column names to the new dataframe
    df.columns = ['A549_0_1', 'A549_0_2', 'A549_0_3','A549_25_1',
                   'A549_25_2','A549_25_3']
    df['gene_name'] = c # a new column will be add to the dataframe,
    # the gene names

    # shift column gene_name to first position
    first_column = df.pop('gene_name')

    # insert column using insert(position,column_name,
    # first_column) function
    df.insert(0, 'gene_name', first_column)
    # make new file with normalized data and good column names in order
    df.to_csv("genormaliseerde_waarden_nieuw.txt", sep='\t', mode='w',
              index=False)


def main():
    raw_counts = "/Users/mushtaaqismail/PycharmProjects/data_science_project/matrix.txt"
    c = test(raw_counts)
    een(raw_counts,c)

main()



# def een(raw_counts, normalised_counts):
#     """
#     The raw counts file is read in using pandas, so that the
#     counts can be normalized with the sklearn.preprocessing.
#     """
#     # read file into pandas dataframe
#     df = pd.read_table(f"{raw_counts}",
#                        sep="\t", usecols=range(1, 7),
#                        skiprows=1)
#
#     x = df.values  # gives numpy array back
#     column = ''
#     min_max_scaler = preprocessing.MinMaxScaler()
#     x_scaled = min_max_scaler.fit_transform(x)  # normalized values
#     df = pd.DataFrame(
#         x_scaled)  # normalized values into pandas dataframe
#     # makes new file with new variables
#     df.to_csv(f"{normalised_counts}", sep='\t', mode='w')
#
#
# def main(raw_counts, normalised_counts):
#     een(raw_counts, normalised_counts)
#
#
# if __name__ == '__main__':
#     main(sys.argv[1], sys.argv[2])
