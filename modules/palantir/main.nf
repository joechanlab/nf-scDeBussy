process PALANTIR {
    label 'process_medium'
    conda "/usersoftware/chanj3/tslearn"
    publishDir "${params.outdir}/palantir/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(patient), path(path)

    output:
    path "palantir.${patient}.ipynb", emit: report_ipynb
    path "palantir.${patient}.ipynb", emit: report_html

    script:
    """
    export HOME=\$PWD
    /usersoftware/chanj3/tslearn/bin/python -m ipykernel install --user --name tslearn --display-name "[conda] tslearn"
    /usersoftware/chanj3/tslearn/bin/papermill ${baseDir}/bin/palantir.ipynb palantir.${patient}.ipynb \
    -p patient ${patient}
    /usersoftware/chanj3/tslearn/bin/jupyter nbconvert --to html palantir.${patient}.ipynb
    """
}
