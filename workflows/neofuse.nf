/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScrnaseq.initialise(params, log)

def checkPathParamList = [
    params.input, params.multiqc_config, params.fasta, params.gtf,
    params.transcript_fasta, params.salmon_index, params.kallisto_index,
    params.star_index, params.txp2gene, params.barcode_whitelist, params.cellranger_index,
    params.universc_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { FASTQC_CHECK      } from '../subworkflows/local/fastqc'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
// TODO: Are this channels still necessary?
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
(protocol, chemistry, other_parameters) = WorkflowScrnaseq.formatProtocol(params.protocol, params.aligner)

// general input and params
ch_input = file(params.input)
ch_genome_fasta = params.fasta ? file(params.fasta) : []
ch_gtf = params.gtf ? file(params.gtf) : []
ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta): []
ch_txp2gene = params.txp2gene ? file(params.txp2gene) : []
ch_multiqc_alevin = Channel.empty()
ch_multiqc_star = Channel.empty()
if (params.barcode_whitelist) {
    ch_barcode_whitelist = file(params.barcode_whitelist)
} else if (params.protocol.contains("10X")) {
    ch_barcode_whitelist = file("$baseDir/assets/whitelist/10x_${chemistry}_barcode_whitelist.txt.gz", checkIfExists: true)
} else {
    ch_barcode_whitelist = []
}


//kallisto params
ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
kb_workflow = params.kb_workflow

workflow NEOFUSE {

    ch_versions     = Channel.empty()
    ch_mtx_matrices = Channel.empty()

    // Check input files and stage input data
    ch_fastq = INPUT_CHECK( ch_input ).reads

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Run FastQC
    ch_multiqc_fastqc = Channel.empty()
    if (!params.skip_fastqc) {
        FASTQC_CHECK ( ch_fastq )
        ch_versions       = ch_versions.mix(FASTQC_CHECK.out.fastqc_version)
        ch_multiqc_fastqc = FASTQC_CHECK.out.fastqc_zip
    } else {
        ch_multiqc_fastqc = Channel.empty()
    }

    ch_filter_gtf = GTF_GENE_FILTER ( ch_genome_fasta, ch_gtf ).gtf