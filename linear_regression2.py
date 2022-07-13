"""
Author: Christel van Haren and Michelle Memelink
This file uses the normalised rawcounts and a test matrix for an eQTL analysis.
"""

import img as img
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from PIL import Image


def counts():
    """
    This file uses the normalised rawcounts file from the
    normalised_counts.py. It makes a 2D list from every row in the file.
    """
    with open("genormaliserde_rawcounts.txt", 'r') as c:
        first_lines = c.readline() #skips the first row
        count = []
        for lines in c:
            inf = lines.strip("\n") #removes the enter at the end
            inf2 = inf.split("\t")
            inf2 = list(map(float, inf2))
            count.append(inf2[1:]) #skips the first column

        return count


def regressieanalyse():
    """
    The file used is a test matrix.
    """
    with open("regressieanalyse test.txt", 'r') as g:
        first_line = g.readline() #skips the first row
        genotype = []
        for line in g:
            info2 = line.strip("\n") #removes the enter at the end
            info = info2.split("\t")
            info = list(map(int, info))
            genotype.append(info[1:]) #skips the first colums

        return genotype


def analysis(genotype, count):
    """
    This function does the eQTL analysis. It calculates the
    coefficient and if there are certain combinations of the
    genotypes than that plot will be skipped.
    """
    for i in genotype:
        for k in count:
            x = np.array(i) #genotype
            y = np.array(k) #phenotype

            if x.__contains__(0) and x.__contains__(
                    1) or x.__contains__(
                0) and x.__contains__(2):

                slope, intercept, r, p, std_err = stats.linregress(x, y)
                if stats.ttest_ind(x, y) <= 0.05:
                    def myfunc(x):
                        """
                        This calculates the coefficient and returns it.
                        """
                        return slope * x + intercept  # y=ax+b

                    mymodel = list(map(myfunc, x))

                    # if the counts(phenotype) have all 0's in a row,
                    # the plot will not be made.
                    if mymodel != [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]:
                        plt.scatter(x, y)
                        plt.plot(x, mymodel, color='black')
                        plt.title("Regression Analysis")
                        plt.xlabel("Genotype")
                        plt.ylabel("Phenotype")
                        # for n in range(0, len(mymodel)):
                        #     img.save(f"RegressionPlots/eQTL{n}.jpg")
                        #     cv2_plt_imshow.imshow(f"RegressionPlots/eQTL{n}.jpg", mymodel)
                        #     cv2.imwrite(f"RegressionPlots/eQTL{n}.jpg", img)
                        plt.show()
                    else:
                        continue
                else:
                    continue
            else:
                continue


def main():
    count = counts()
    genotype = regressieanalyse()
    analysis(genotype, count)


main()
