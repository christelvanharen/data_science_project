"""
Author: Michelle Memelink, Christel van Haren and Gijsbert Keja
This file uses the normalised rawcounts and a test matrix for an eQTL analysis.
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def counts(filename_rawcounts):
    """
    This file uses the normalised rawcounts file from the
    normalised_counts.py. It makes a 2D list from every row in the file.
    """
    counter = 0
    with open(filename_rawcounts, 'r') as c:
        first_lines = c.readline() #skips the first row
        count = []
        for lines in c:
            inf = lines.strip("\n") #removes the enter at the end
            inf2 = inf.split("\t")
            inf2 = list(map(float, inf2[1:])) #skips the first column and makes every element a float
            count.append(inf2)
            print(inf2)
        return count


def regressieanalyse(filename_genotypes):
    """
    The file used is a test matrix.
    """
    with open(filename_genotypes, 'r') as g:
        first_line = g.readline() #skips the first row
        genotype = []
        for line in g:
            info2 = line.strip("\n") #removes the enter at the end
            info = info2.split("\t")
            info = list(map(int, info[1:])) #skips the first colums
            print(info)
            genotype.append(info)

        return genotype


def analysis(genotype, count, directory):
    """
    This function does the eQTL analysis. It calculates the
    coefficient and if there are certain combinations of the
    genotypes than that plot will be skipped.
    """
    counts = 0
    for i in genotype:
        for k in count:
            x = np.array(i) #genotype
            y = np.array(k) #phenotype

            if x.__contains__(0) and x.__contains__(
                    1) or x.__contains__(
                0) and x.__contains__(2):

                slope, intercept, r, p, std_err = stats.linregress(x, y)

                def myfunc(x):
                    """
                    This calculates the coefficient and returns it.
                    """
                    return slope * x + intercept #y=ax+b

                mymodel = list(map(myfunc, x))

                # if the counts(phenotype) have all 0's in a row,
                # the plot will not be made.
                if mymodel != [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]:
                    plt.scatter(x, y)
                    plt.plot(x, mymodel, color='black')
                    plt.title("Regression Analysis")
                    plt.xlabel("Genotype")
                    plt.ylabel("Phenotype")
                    plt.savefig(f"{directory}/eQTL{counts}.jpg")
                    plt.clf()
                    plt.close()
                    counts += 1
                else:
                    continue
            else:
                continue


def main():
    filename_rawcounts = "final_normalised_counts.txt"
    filename_genotypes = "resultaten_snpsel_definitief_testen.txt"
    directory = "Plaats hier je directory waar je je plotjes wilt hebben"
    count = counts(filename_rawcounts)
    genotype = regressieanalyse(filename_genotypes)
    analysis(genotype, count, directory)

main()
