process CELLALIGNDTW {
    conda '/usersoftware/chanj3/tslearn'
    publishDir "${params.outdir}/CellAlignDTW/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(cluster_ordering), path(path), val(gene_list)

    output:
    val cluster_ordering, emit: cluster_ordering
    path cluster_ordering, emit: output_path

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    mkdir ${cluster_ordering}
    /usersoftware/chanj3/tslearn/bin/python ${baseDir}/bin/run_CellAlignDTW.py \
        --path ${path} \
        --outpath ${cluster_ordering} \
        --clusters ${cluster_ordering} \
        --gene_list ${gene_list} \
        --n_splines ${params.n_splines} \
        --lam ${params.lam}
    """

}