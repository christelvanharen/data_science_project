"""
Author: Michelle Memelink, Christel van Haren,  Gijsbert Keja and
Moshtach Ismail
This file uses the normalised rawcounts and a test matrix for an eQTL analysis.
"""
import sys
import matplotlib.pyplot as plt
import numpy
from scipy import stats

def filter(filename_rawcounts, filename_genotypes):
    """
    This function looks at both files and gets all the counts from the
    gene names that are in snpsel_definitief and puts these in a 2D
    list.
    :return:
    """
    count = [] # lijst met alle counts van de gennamen waarvan
    # genotype bekend is
    key_words = []
    next = []
    with open(filename_genotypes, 'r') as b:
        for lines in b:
            keywords = lines.split("\t")[0]
            key_words.append(keywords)

    with open(filename_rawcounts, 'r') as f:
        all_lines = f.read().splitlines()
        for line in all_lines:
            if any(kw in line for kw in key_words): # if any of the word in key words is in line print it
                inf = line.strip("\n")  # removes the enter at the end
                inf2 = inf.split("\t")
                inf2 = list(map(float, inf2[1:]))  # skips the first column and makes every element a float
                print(inf2)
                count.append(inf2)
              #  next.append(kw)

    return count, key_words

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
            genotype.append(info)

        return genotype

def analysis(genotype, count, directory, key_words):
    """
    This function does the eQTL analysis. It calculates the
    coefficient and if there are certain combinations of the
    genotypes than that plot will be skipped.
    """
    counts = 0

    arr = numpy.array(genotype)
    arr2 = numpy.array(count)
    for x, y, a in zip(arr, arr2, key_words[1:]):
        if x.__contains__(0) and x.__contains__(1) or x.__contains__(0) and x.__contains__(2):
            slope, intercept, r, p, std_err = stats.linregress(x,y)

            def myfunc(x):
                """
                This calculates the coefficient and returns it.
                """
                return slope * x + intercept  # y=ax+b

            mymodel = list(map(myfunc, x))

                # if the counts(phenotype) have all 0's in a row,
                # the plot will not be made.
            if mymodel != [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]:
                plt.scatter(x, y, color="lime", label="cancer patients ")
                plt.plot(x, mymodel, color='black')
                plt.plot(x[0], y[0], 'ob',label="not-cancer patients")
                # first element are cancerpatients
                plt.plot(x[-1], y[-1], 'ob') #
                # last element are cancerpatients
                plt.title(f"Regression Analysis {a}") #gene name in 
                # title
                plt.xlabel("Genotype")
                plt.legend(loc="upper left") # legenda toevoegen
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
    filename_rawcounts = "/Users/mushtaaqismail/PycharmProjects/" \
                         "data_science_project/final_normalised_counts.txt"
    filename_genotypes = "/Users/mushtaaqismail/PycharmProjects/" \
                         "data_science_project/resultaten_snpsel_definitief.txt"
    directory = "/Users/mushtaaqismail/PycharmProjects" \
                "/data_science_project/plt"
    count, key_words = filter(filename_rawcounts, filename_genotypes)
    genotype = regressieanalyse(filename_genotypes)
    analysis(genotype, count, directory, key_words)

main()
