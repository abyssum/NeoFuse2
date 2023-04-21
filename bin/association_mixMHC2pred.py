#!/usr/bin/python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

import argparse, sys, os
import csv

# function to create files with peptides in mixMHC2pred input format
def get_peps(
    xenoInFile,
    fusionGene,
    bkpoint1,
    bkpoint2,
    split_read1,
    split_read2,
    disco_reads,
    closest_bkp1,
    closest_bkp2,
    peps,
    ftype,
    confidence,
    stopCod,
    out_files,
):

    with open(xenoInFile) as in_file:
        for line in in_file:
            gene1 = line.split("#", 1)[0].split("\t", 1)[0].split(" - ")[0]
            gene2 = line.split("#", 1)[0].split("\t", 1)[0].split(" - ")[1]
            fgene = f"{gene1}_{gene2}"
            fusionGene.append(fgene)
            peps.append(
                line.split("#", 1)[0].split("\t", 1)[1].replace("\t", "\n").upper()
            )
            ftype.append(line.split("#")[1])
            confidence.append(line.split("#")[2])
            stopCod.append(line.split("#")[3])
            bkpoint1.append(line.split("#")[4])
            bkpoint2.append(line.split("#")[5])
            split_read1.append(line.split("#")[6])
            split_read2.append(line.split("#")[7])
            disco_reads.append(line.split("#")[8])
            closest_bkp1.append(line.split("#")[9])
            closest_bkp2.append(line.split("#")[10].replace("\n", ""))

    fileID = 1
    for i in range(0, len(fusionGene)):
        tag = fusionGene[i].split(",")[0].split("(")[0] + "_" + str(fileID)
        out_file = xenoInFile.replace("xeno", "%s" % tag).replace("tsv", "txt")
        with open(out_file, "w") as mmp_out:
            mmp_out.write(peps[i])
            out_files.append(out_file)
            fileID += 1

    return (
        fusionGene,
        bkpoint1,
        bkpoint2,
        split_read1,
        split_read2,
        disco_reads,
        closest_bkp1,
        closest_bkp2,
        peps,
        ftype,
        confidence,
        stopCod,
        out_files,
    )
    in_file.close()

# Function to retrive HLAs predicted by HLAHD
def seek_hla(hlahdInFile):
    with open(hlahdInFile) as in_file:
        hl = []
        for line in in_file:
            h = line.replace("\n", " ")
            hl.append(h)
        hla = " ".join(hl)
    return hla


# Function to build the associations and mixMHC2pred run intermidiate files
def tmp_out_pep(
    out_files,
    tmpOutFile1,
    tmpOutFile2,
    outDir,
    fGenes,
    bkpoint1,
    bkpoint2,
    split_read1,
    split_read2,
    disco_reads,
    closest_bkp1,
    closest_bkp2,
    peptides,
    hla,
    ftype,
    confidence,
    stopCod,
    cores,
):
    gene_file = []
    counter = 2  # used for parallelizing mixMHC2pred jobs
    fileID = 1  # used in printing fusion gene names and corresponding peptide sequences at the MHCFlurry run file
    with open(tmpOutFile1, "+w") as out_file:
        for i in range(0, len(fGenes)):
            out = (
                outDir
                + "_"
                + fGenes[i].split(",")[0].split("(")[0]
                + "_"
                + str(fileID)
                + ".tsv"
            )
            if not hla:  # check if there are HLAs supported by mixMHC2pred
                pass
            else:
                gene_file.append(
                    fGenes[i]
                    + "#"
                    + out
                    + "#"
                    + ftype[i]
                    + "#"
                    + confidence[i]
                    + "#"
                    + stopCod[i]
                    + "#"
                    + bkpoint1[i]
                    + "#"
                    + bkpoint2[i]
                    + "#"
                    + split_read1[i]
                    + "#"
                    + split_read2[i]
                    + "#"
                    + disco_reads[i]
                    + "#"
                    + closest_bkp1[i]
                    + "#"
                    + closest_bkp2[i]
                )  # get the lines for the associations file
                if counter <= cores:
                    out_file.write(
                        """MixMHC2pred -i %s -o %s -a %s &\n"""
                        % (out_files[i], out, hla)
                    )
                else:
                    out_file.write(
                        """MixMHC2pred -i %s -o %s -a %s &\n"""
                        % (out_files[i], out, hla)
                    )
                    out_file.write("wait\n")
                    counter = 1
            counter += 1
            fileID += 1

    # write the associations file - to be used later for multiplexing the output files into a single final output
    with open(tmpOutFile2, "+w") as out_file2:
        for j in gene_file:
            out_file2.write("%s\n" % j)


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        usage="association.py [-h] -x {/path/to/xeno/input/} -l {path/to/OptiType/input} -o {/path/to/output/} -c (NCores)"
    )
    parser.add_argument("-x", "--XenoFile", help="Xeno file", required=True)
    parser.add_argument("-l", "--hlahdFile", help="Optitype file", required=True)
    parser.add_argument("-o", "--outDir", help="Output dir", required=True)
    parser.add_argument(
        "-c", "--cores", help="Number of cores", type=int, required=False
    )
    args = parser.parse_args()
    in_file_xeno = args.XenoFile
    hlahd_file = args.hlahdFile
    outDir = args.outDir
    cores = args.cores
    fusionGene = []
    peps = []
    ftype = []
    confidence = []
    stopCod = []
    bkpoint1 = []
    bkpoint2 = []
    split_read1 = []
    split_read2 = []
    disco_reads = []
    closest_bkp1 = []
    closest_bkp2 = []
    out_files = []
    # hla=[]
    outFile1 = outDir + "_TMP_OUT.sh"
    outFile2 = outDir + "_ASSOCIATIONS_OUT.txt"

    get_peps(
        in_file_xeno,
        fusionGene,
        bkpoint1,
        bkpoint2,
        split_read1,
        split_read2,
        disco_reads,
        closest_bkp1,
        closest_bkp2,
        peps,
        ftype,
        confidence,
        stopCod,
        out_files,
    )
    hla = seek_hla(hlahd_file)
    tmp_out_pep(
        out_files,
        outFile1,
        outFile2,
        outDir,
        fusionGene,
        bkpoint1,
        bkpoint2,
        split_read1,
        split_read2,
        disco_reads,
        closest_bkp1,
        closest_bkp2,
        peps,
        hla,
        ftype,
        confidence,
        stopCod,
        cores,
    )
