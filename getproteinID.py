"""
Author: Christel van Haren
Description: This file extracts the proteins IDs from each chromosome
file.
Date: 16-07-2022
"""
import sys

# input file = chr10_protein.txt
# output file = proteinIDs.txt

def protein(input_file, output_file):
    """
    Extracts the accession ID from the FASTA file of the chromosome.
    """
    proteinIDs = []
    with open(input_file, "r") as c:
        for line in c:
            if line.startswith(">"):
                ids = line.split("|")
                proteinIDs.append(ids[1])
                write_file = open(output_file, "w")
                write_file.write(str(proteinIDs))
                # correct adjustments made for in file
                with open(output_file, "r") as p:
                    for lines in p:
                        id = lines.replace("', '", "\n")
                        ids = id.removeprefix("['").removesuffix("']")
                        print(ids)
                        write_file = open(output_file, "w")
                        write_file.write(ids)
                write_file.close()


def main(input_file, output_file):
    protein(input_file, output_file)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
