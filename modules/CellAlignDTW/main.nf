process CELLALIGNDTW {
    conda '/usersoftware/chanj3/tslearn'
    publishDir "${params.outdir}/CellAlignDTW/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(cluster_ordering), path(path)

    output:
    val cluster_ordering, emit: cluster_ordering
    path "${cluster_ordering}_summary_df.csv", emit: summary_df_path
    path "${cluster_ordering}_aggregated_curves.csv", emit: aggregated_curves_path
    path "${cluster_ordering}_scores_df.csv", emit: scores_df_path

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/run_CellAlignDTW.py \
        --path ${path} \
        --clusters ${cluster_ordering} \
        --gene_list ${params.gene_list} \
    """

}