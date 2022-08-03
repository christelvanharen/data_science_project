"""
 Author: Christel van Haren
 Date: 16 july 2022
 """

def protein():
    try:
        proteinIDs = []
        with open("chr10_protein.txt", 'r') as c:
            for line in c:
                if line.startswith(">"):
                    ids = line.split("|")
                    gn = ids[2].split("GN=")
                    gen_id = gn[1].split(" ")
                    # print(gen_id[0])
                    proteinIDs.append(gen_id[0])
                    duplicates = list(set(proteinIDs))
                    write_file = open("proteinIDs.txt", "w")
                    write_file.write(str(duplicates[:600]))
                write_file.close()
    except IndexError:
        pass


def file():
    with open("proteinIDs.txt", "r") as p:
        for line in p:
            gen = line.replace("', '", "\n")
            ids = gen.removeprefix("['").removesuffix("']")
            print(ids)
            write_file = open("proteinIDs.txt", "w")
            write_file.write(ids)
        write_file.close()


def main():
    protein()
    file()


main()