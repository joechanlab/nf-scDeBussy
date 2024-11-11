process LIMMA {
    conda '/usersoftware/chanj3/ArchR'
    publishDir "${params.outdir}/limma/", mode: 'copy'
    cache 'lenient'

    input:
    val cluster_ordering
    path path

    output:
    val cluster_ordering, emit: cluster_ordering
    path cluster_ordering, emit: output_path

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    mkdir -p ${cluster_ordering}
    /usersoftware/chanj3/ArchR/bin/Rscript ${baseDir}/bin/limma.R \
        ${path} \
        ${cluster_ordering} \
        ${cluster_ordering}
    """
}
