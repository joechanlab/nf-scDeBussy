#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
include { PREPROCESSING } from './modules/preprocessing'
include { LIMMA } from './modules/limma'
include { CELLALIGNDTW } from './modules/CellAlignDTW'
include { REPORT } from './modules/report'

workflow {
    // Access the samplesheet
    sample_sheet = file(params.samplesheet)

    // Read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // Create a channel from the paths
    out = Channel.from(sample_sheet_data).map { row ->
        def (cluster_ordering, path) = row
        return tuple(cluster_ordering, path)
    }

    PREPROCESSING(out)

    LIMMA(PREPROCESSING.out.cluster_ordering, PREPROCESSING.out.output_path)
    
    CELLALIGNDTW(LIMMA.out.cluster_ordering, LIMMA.out.output_path)

    //REPORT(CELLALIGNDTW.out.cluster_ordering, CELLALIGNDTW.out.output_path)

}
