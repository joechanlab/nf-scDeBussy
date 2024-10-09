#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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

    CELLALIGNDTW(out)

    REPORT(CELLALIGNDTW.out.cluster_ordering, CELLALIGNDTW.out.summary_df_path, CELLALIGNDTW.out.aggregated_curves_path, CELLALIGNDTW.out.scores_df_path)

}
