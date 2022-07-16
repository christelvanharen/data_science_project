def protein():
    """
    """
    proteinIDs = []
    with open("chr10_protein.txt", 'r') as c:
        for line in c:
            if line.startswith(">"):
                ids = line.split("|")
                # print(ids[1])
                proteinIDs.append(ids[1])
                write_file = open("proteinIDs.txt", "w")
                write_file.write(str(proteinIDs))
                write_file.close()


def main():
    protein()


main()