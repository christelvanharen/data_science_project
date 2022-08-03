import pandas as pd
from sklearn import preprocessing


def een():
    df = pd.read_table("matrix.txt",
                       sep="\t", usecols=range(1, 7),
                       skiprows=1)
    # print(df)
    x = df.values  # returns a numpy array
    column = ''
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled)
    # print(df)
    df.to_csv("genormaliseerde_rawcounts_output.txt", sep='\t',
              mode='w')


def main():
    een()


main()
