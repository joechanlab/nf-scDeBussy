process PREPROCESSING {
    conda '/usersoftware/chanj3/tslearn'
    publishDir "${params.outdir}/preprocessing/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(cluster_ordering), path(path)

    output:
    val cluster_ordering, emit: cluster_ordering
    path cluster_ordering, emit: output_path

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    export SCIPY_ARRAY_API=1
    mkdir -p ${cluster_ordering}
    /usersoftware/chanj3/tslearn/bin/python ${baseDir}/bin/preprocessing.py \
        --path ${path} \
        --outpath ${cluster_ordering} \
        --clusters ${cluster_ordering} \
    """
}