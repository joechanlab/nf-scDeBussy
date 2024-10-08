#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CELLALIGNDTW } from './modules/CellAlignDTW'

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

    CELLALIGNDTW(out)

}
