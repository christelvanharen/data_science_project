"""
Author: Christel van Haren
Date: 16 july 2022
"""

def protein():
    proteinIDs = []
    with open("chr10_protein.txt", "r") as c:
        for line in c:
            if line.startswith(">"):
                ids = line.split("|")
                # print(ids[1])
                proteinIDs.append(ids[1])
                write_file = open("proteinIDs.txt", "w")
                write_file.write(str(proteinIDs))
                write_file.close()

def file():
    with open("proteinIDs.txt", "r") as p:
        for line in p:
            id = line.replace("', '", "\n")
            ids = id.removeprefix("['").removesuffix("']")
            print(ids)
            write_file = open("proteinIDs.txt", "w")
            write_file.write(ids)
            write_file.close()


def main():
    protein()
    file()


main()