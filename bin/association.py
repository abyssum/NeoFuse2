#!/usr/bin/python
"""
Requirements:
    * Python >= 3.6.2

Copyright (c) 2020 Georgios Fotakis <georgios.fotakis@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

import argparse, sys
import csv

# Function to retrieve fusion gene names, peptide sequences and confidence level
def seek_fusePep(
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
):
    with open(xenoInFile) as in_file:
        for line in in_file:
            gene1 = line.split("#", 1)[0].split("\t", 1)[0].split(" - ")[0]
            gene2 = line.split("#", 1)[0].split("\t", 1)[0].split(" - ")[1]
            fgene = f"{gene1}_{gene2}"
            fusionGene.append(fgene)
            peps.append(
                line.split("#", 1)[0].split("\t", 1)[1].replace("\t", " ").upper()
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
    in_file.close()

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
    )


# Function to retrive HLAs predicted by OptiType
def seek_hla(optiInFile, outDir, hla):
    with open(optiInFile) as in_file:
        hl = []
        uns = []
        for line in in_file:
            h = "HLA-%s" % line
            h = h.replace("\n", " ")
            hl.append(h)
        hla = " ".join(hl)
        unsupported = "\n".join(uns)
    in_file.close()
    out = outDir + "_unsupported.txt"
    with open(out, "+w") as out:
        out.write(
            "The following HLA types were predicted by OptiType, but are not supported by MHCFlurry:\n"
        )
        out.write(unsupported)
    out.close()
    return hla


# Function to build the associations and Mhcflurry run intermidiate files
def tmp_out_pep(
    xenoInFile,
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
    counter = 2  # used for parallelizing mhcflurry jobs
    fileID = 1  # used in printing fusion gene names and corresponding peptide sequences at the MHCFlurry run file
    if "xeno_8.tsv" in xenoInFile:
        postfix = "_8.tsv"
    elif "xeno_9.tsv" in xenoInFile:
        postfix = "_9.tsv"
    elif "xeno_10.tsv" in xenoInFile:
        postfix = "_10.tsv"
    elif "xeno_11.tsv" in xenoInFile:
        postfix = "_11.tsv"
    with open(tmpOutFile1, "+w") as out_file:
        for i in range(0, len(fGenes)):
            out = (
                outDir
                + "_"
                + fGenes[i].split(",")[0].split("(")[0]
                + "_"
                + str(fileID)
                + postfix
            )
            if not hla:  # check if there are supported by MHCFlurry HLAs
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
                        """mhcflurry-predict --affinity-only --alleles %s --peptides %s --out %s --models /home/neofuse/.local/share/mhcflurry/4/2.0.0/models_class1_pan/models.combined &\n"""
                        % (hla, peptides[i], out)
                    )
                else:
                    out_file.write(
                        """mhcflurry-predict --affinity-only --alleles %s --peptides %s --out %s --models /home/neofuse/.local/share/mhcflurry/4/2.0.0/models_class1_pan/models.combined &\n"""
                        % (hla, peptides[i], out)
                    )
                    out_file.write("wait\n")
                    counter = 1
            counter += 1
            fileID += 1
    out_file.close()
    # write the associations file - to be used later for multiplexing the output files into a single final output
    with open(tmpOutFile2, "+w") as out_file2:
        for j in gene_file:
            out_file2.write("%s\n" % j)
    out_file2.close()


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        usage="association.py [-h] -x {/path/to/xeno/input/} -l {path/to/OptiType/input} -o {/path/to/output/} -c (NCores)"
    )
    parser.add_argument("-x", "--XenoFile", help="Xeno file", required=True)
    parser.add_argument("-l", "--OptiFile", help="Optitype file", required=True)
    parser.add_argument("-o", "--outDir", help="Output dir", required=True)
    parser.add_argument(
        "-c", "--cores", help="Number of cores", type=int, required=False
    )
    args = parser.parse_args()
    fGene = []
    bkpoint1 = []
    bkpoint2 = []
    split_read1 = []
    split_read2 = []
    disco_reads = []
    closest_bkp1 = []
    closest_bkp2 = []
    peps = []
    ftype = []
    confidence = []
    stopCod = []
    hla = ""
    xenoFile = args.XenoFile
    optiInFile = args.OptiFile
    outDir = args.outDir
    cores = args.cores
    outFile1 = outDir + "_TEST_OUT.sh"
    outFile2 = outDir + "_ASSOCIATIONS_OUT.txt"

    seek_fusePep(
        xenoFile,
        fGene,
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
    )
    hla = seek_hla(optiInFile, outDir, hla)

    tmp_out_pep(
        xenoFile,
        outFile1,
        outFile2,
        outDir,
        fGene,
        bkpoint1,
        bkpoint2,
        split_read1,
        split_read2,
        disco_reads,
        closest_bkp1,
        closest_bkp1,
        peps,
        hla,
        ftype,
        confidence,
        stopCod,
        cores,
    )
