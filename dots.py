def dots(file, new_file):
    write_file = open(new_file, "w")
    with open(file, "r") as c:
        for line in c:
            gene_name = line.split("\t")
            # print(gene_name[0])
            no_dots = gene_name[0].split(".")
            # print(no_dots[0])
            new_column = line.replace(gene_name[0], no_dots[0])
            print(new_column)
            write_file.write(new_column)
    write_file.close()

def main():
    file = "genormaliseerde_waarden_nieuw.txt"
    new_file = "final_normalised_counts.txt"
    dots(file, new_file)


main()