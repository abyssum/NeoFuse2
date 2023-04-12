# NeoFuse2
NeoFuse is a user-friendly pipeline for the prediction of fusion neoantigens from tumor RNA-seq data.

NeoFuse takes single-sample FASTQ files of RNA-seq reads (single- or paired-end) as input and predicts putative fusion neoantigens through five main analytical modules based on state-of-the-art computational tools:

* Genotyping of class-I and II Human Leukocyte Antigen (HLA) genes at 4-digit resolution using OptiType (Szolek et al., 2014) and HLA-HD (cite).
* Prediction of fusion peptides using Arriba (cite), together with confidence scores reflecting the likelihood that a fusion is caused by a tumor-specific genomic rearrangement and is not due to technical artifacts.
* Prediction of the binding affinity of fusion peptides to HLA types, quantified as half maximal inhibitory concentration (IC50) and percentile rank, using MHCflurry2 (O’Donnell et al., 2018), netMHCpan (Jurtz et al., 2017), and netMHCpanII (cite).
* Quantification of gene expression levels, as transcripts per million (TPM), using STAR (Dobin et al., 2013) and featureCounts (Liao et al., 2014).
* Neoantigen prioritization based on IC50 binding affinity and confidence score, and annotation of each neoantigen with: IC50, percentile rank, confidence score, binding HLA type, expression of the fusion and HLA genes in TPM, and information about the presence of a premature stop codon that might cause nonsense mediated decay of the fusion transcript.
