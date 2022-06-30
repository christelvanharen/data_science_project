"""
Author: Gabe van den Hoeven
Date: 10-06-2022
Description: This script is used to search for proteins that have
interactions with any of a given list of proteins and to find Single
Nucleotide Polymorphisms (SNPs) in these interaction proteins using the dbSNP.
If the snp is considered to be pathogenic in the dbSNP it is put in a list
that is later made into a TSV file.
"""

import datetime
import http.client
import re
from Bio import Entrez
import stringdb
import sys

import id2ensemblID


class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def __getattr__(self, attr):
        return getattr(self.terminal, attr)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)


def main(filename_in, filename_out):
    sys.stdout = Logger("console_log.txt")
    start_time = datetime.datetime.now()
    print(f"Started script at {start_time}.")
    genes = []
    print("Retrieving ribosomal protein names...")
    with open(filename_in, "r") as file:
        for line in file:
            line = line.strip()
            genes.append(line)
    interaction_proteins = get_proteins(genes)
    all_snps = get_data(interaction_proteins)
    write_to_file(all_snps, filename_out)
    end_time = datetime.datetime.now()
    print(f"Done at {end_time}. Results can be found in"
          f" {filename_out}.")
    print(f"Time elapsed while script was running: {end_time - start_time}.")


def get_proteins(genes):
    """Takes a list of gene symbols, converts them to ensembl IDs,
    then finds interaction proteins by querying the stringdb.

    :parameter: genes - list of gene symbols for ribosomal proteins.
    :return: all_interaction_proteins - list of all unique
    ribosome-interaction proteins.
    """
    ensembl_ids = []
    result = id2ensemblID.id2ensembl(genes, genus_name="homo",
                                     species_name="sapiens")
    for ensembl_id in result:
        ensembl_ids.append(ensembl_id[2])
    print("Retrieving interaction proteins...")
    result = stringdb.get_interaction_partners(ensembl_ids, limit=1000,
                                               required_score=900)
    protein_list = result["preferredName_B"]
    all_interaction_proteins = list(dict.fromkeys(protein_list))
    print(f"{len(all_interaction_proteins)} proteins were found.")
    return all_interaction_proteins


def get_data(proteins):
    """For each interaction protein, retrieves all SNPs from the SNPdb.

    :parameter: proteins - list of interaction proteins.
    :return: all_snps - nested list with information of each SNP that was
    found.
    """
    print("Searching for SNPs...")
    all_snps = []
    for protein in proteins:
        print(f"Current protein: {protein}. Index: {proteins.index(protein)}")
        try:
            Entrez.email = ""
            handle = Entrez.esearch(db="SNP", term=protein)
            record = Entrez.read(handle)
            count = int(record["Count"])
            retmax = 10000
            retstart = 0
            while retstart < count:
                handle = Entrez.esearch(db="SNP", term=protein, retmax=retmax,
                                        retstart=retstart)
                record = Entrez.read(handle)
                id_list = record["IdList"]
                handle = Entrez.efetch(db="SNP", id=id_list, rettype="xml",
                                       retmax=retmax)
                try:
                    snp_data = handle.read()
                except http.client.IncompleteRead as ir:
                    snp_data = ir.partial
                snp_data = snp_data.decode("utf-8")
                snp_data = snp_data.split("\n")
                del snp_data[0]
                del snp_data[-1]
                print(f"{len(snp_data)} SNPs were found.")
                all_snps = filter_snp_data(snp_data, all_snps)
                retstart += retmax
        except ValueError as ve:
            print(f"ValueError: {ve.__traceback__}")
        except IndexError as ie:
            print(f"IndexError: {ie.__traceback__}")
    return all_snps


def filter_snp_data(snp_data, all_snps):
    """Filters the SNP data in a list from XML format.

    :parameter: snp_data - list unfiltered SNP data.
    :return: all_snps - nested list with information for each SNP that was
    found.
    """
    for line in snp_data:
        tmp_list = []
        line = line.split("><")
        for element in line:
            if element.startswith("SNP_ID>"):
                tmp_list.append(element.replace("SNP_ID>", "").replace(
                    "</SNP_ID", ""))
            elif element.startswith("CLINICAL_SIGNIFICANCE"):
                tmp_list.append(element.replace("CLINICAL_SIGNIFICANCE>", "").
                                replace("</CLINICAL_SIGNIFICANCE", ""))
            elif element.startswith("NAME>"):
                tmp_list.append(element.replace("NAME>", "").replace(
                    "</NAME", ""))
            elif element.startswith("CHRPOS>"):
                tmp_list.append(element.replace("CHRPOS>", "").replace(
                    "</CHRPOS", ""))
            elif element.startswith("SNP_CLASS>"):
                tmp_list.append(element.replace("SNP_CLASS>", "").replace(
                    "</SNP_CLASS", ""))
            elif element.startswith("DOCSUM>"):
                match = re.search("SEQ=\[(.+)]", element)
                if match:
                    tmp_list.append(match.group(1))
        if not tmp_list == []:
            all_snps.append(tmp_list)
    return all_snps


def write_to_file(file_lines, filename_out):
    """Writes filtered SNP data from a list to a new file.

    :parameter: file_lines - list with lines to be written to the file.
    """
    print("Writing results to file...")
    with open(filename_out, "w") as file:
        file.write("SNP_ID\tClinical_significance\tGene\tMutation"
                   "\tVariation\tChromosome_Position\n")
        for line in file_lines:
            try:
                if line[1] == "pathogenic" and line[4] == "snv":
                    new_line = f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t" \
                               f"{line[4]}\t{line[5]}\n"
                    file.write(new_line)
                    print(f"Writing SNP: {line[0]}. Index: {file_lines.index(line)}")
            except IndexError as e:
                print(f"IndexError. Line {file_lines.index(line)}\n{line}")
                pass


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
