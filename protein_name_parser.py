"""
Author: Gijsbert Keja
This file parses the gennames given in a .txt file and writes
the gennames to a new file
"""

def main(chr_bestand, chr_output):
    lijst_zonder_duplicates = file_reader_splitter(chr_bestand)
    file_writer(chr_output, lijst_zonder_duplicates)

def file_reader_splitter(chr_bestand):
    """
    This file parses the file with the gene names and also removes the
    duplicates
    :param chr_bestand: This is the name of the file which needs to be read
    :type chr_bestand: str
    :return: lijst_zonder_duplicates
    :rtype: list
    """
    chr_lezen = open(chr_bestand, 'r')
    lijst_met_duplicates = []
    lijst_zonder_duplicates = []
    teller = 0
    for x in chr_lezen:
        if x.startswith('>'):
            teller = teller + 1
            line = x.split("|") 
            lijn = line[2].split()
            positie_gennaam = len(lijn) - 3 # de vaste positie van het gennaam is van de laatste positie in de zin 3 naar links
            gennaam = lijn[positie_gennaam]
            gen = gennaam.split("=") # in de gennaam zat er de = voor dus deze moet weg
            lijst_met_duplicates.append(gen[1])
    lijst_zonder_duplicates = list(set(lijst_met_duplicates)) # door het gebruik van een set worden alle duplicaten verwijdert
    return lijst_zonder_duplicates

def file_writer(chr_output, lijst_zonder_duplicates):
    """
    :param chr_output: This is the name of the file which needs to be created
    :type chr_output: str
    :param lijst_zonder_duplicates: this is a list with all gene names
    :type lijst_zonder_duplicates: list
    :return: N/A
    :rtype: N/A
    """
    new_file = open(chr_output, "w")
    for x in lijst_zonder_duplicates:
        new_file.write(x)
        new_file.write("\n") # voegt enters toe aan het nieuwe bestand want anders worden alle string achter elkaar geplakt

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])




