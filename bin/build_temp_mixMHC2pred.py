#!/usr/bin/python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

import argparse, sys
import os.path
import time
import csv


def final_out(infileAssoc, outfile):
    fGene = []
    files = []
    gene1 = []
    gene2 = []
    bkpoint1 = []
    bkpoint2 = []
    split_read1 = []
    split_read2 = []
    disco_reads = []
    closest_bkp1 = []
    closest_bkp2 = []
    ftype = []
    confidence = []
    stopCod = []
    alleles = []
    peptide = []
    rank = []
    gene1 = []
    gene2 = []
    i = 0

    with open(infileAssoc) as in_file:
        for line in in_file:
            fGene.append(line.split("#")[0])
            gene1.append(line.split("#")[0].split("_")[0])
            gene2.append(line.split("#")[0].split("_")[1])
            files.append(line.split("#")[1])
            ftype.append(line.split("#")[2])
            confidence.append(line.split("#")[3])
            stopCod.append(line.split("#")[4])
            bkpoint1.append(line.split("#")[5])
            bkpoint2.append(line.split("#")[6])
            split_read1.append(line.split("#")[7])
            split_read2.append(line.split("#")[8])
            disco_reads.append(line.split("#")[9])
            closest_bkp1.append(line.split("#")[10])
            closest_bkp2.append(line.split("#")[11].replace("\n", ""))

    with open(outfile, "+w") as out_file:
        out_file.write(
            "Fusion\tGene1\tGene2\tBreakpoint1\tBreakpoint2\tSplit_Reads1\tSplit_Reads2\tDiscordant_Reads\tClosest_Breakpoint1\tClosest_Breakpoint2\tHLA_Type\tFusion_Peptide\tRank\tEvent_Type\tStop_Codon\tConfidence\n"
        )
        for assoc_file in files:
            tries = 0
            while not os.path.exists(assoc_file):
                time.sleep(1)
                if tries == 600:
                    break
                tries = tries + 1

            with open(assoc_file) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter="\t")
                for row in csv_reader:
                    if row[0].startswith("#") or row[0].startswith("Peptide"):
                        pass
                    else:
                        if "__" in row[1]:
                            if len(row[1].split("__")) == 2:
                                allele1 = (
                                    row[1]
                                    .split("__")[0]
                                    .replace("_", "*", 1)
                                    .replace("_", ":")
                                )
                                allele2 = (
                                    row[1]
                                    .split("__")[1]
                                    .replace("_", "*", 1)
                                    .replace("_", ":")
                                )
                                allele = allele1 + "/" + allele2
                                alleles.append(allele)
                            elif len(row[1].split("__")) == 3:
                                allele1 = (
                                    row[1]
                                    .split("__")[0]
                                    .replace("_", "*", 1)
                                    .replace("_", ":")
                                )
                                allele2 = (
                                    row[1]
                                    .split("__")[1]
                                    .replace("_", "*", 1)
                                    .replace("_", ":")
                                )
                                allele3 = (
                                    row[1]
                                    .split("__")[2]
                                    .replace("_", "*", 1)
                                    .replace("_", ":")
                                )
                                allele = allele1 + "/" + allele2 + "/" + allele3
                                alleles.append(allele)
                        elif "__" not in row[1]:
                            allele = row[1].replace("_", "*", 1).replace("_", ":")
                            alleles.append(allele)
                        peptide.append(row[0])
                        rank.append(row[2])
                        out_file.write(
                            fGene[i].replace("_", "-")
                            + "\t"
                            + gene1[i]
                            + "\t"
                            + gene2[i]
                            + "\t"
                            + bkpoint1[i]
                            + "\t"
                            + bkpoint2[i]
                            + "\t"
                            + split_read1[i]
                            + "\t"
                            + split_read2[i]
                            + "\t"
                            + disco_reads[i]
                            + "\t"
                            + closest_bkp1[i]
                            + "\t"
                            + closest_bkp2[i]
                            + "\t"
                            + allele
                            + "\t"
                            + row[0]
                            + "\t"
                            + row[2]
                            + "\t"
                            + ftype[i]
                            + "\t"
                            + stopCod[i]
                            + "\t"
                            + confidence[i]
                            + "\n"
                        )
            i += 1


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        usage="build_temp.py [-h] -a {path/to/Associations.txt/input} -o {/path/to/output/}"
    )
    parser.add_argument(
        "-a", "--AssociationsFile", help="Assosciations tmp file", required=True
    )
    parser.add_argument("-o", "--outDir", help="Output dir", required=True)
    args = parser.parse_args()

    inFile = args.AssociationsFile
    outFile = args.outDir + "_final.tsv"

    final_out(inFile, outFile)
