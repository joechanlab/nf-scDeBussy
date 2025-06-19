process REPORT {
    label 'process_medium'
    conda "/usersoftware/chanj3/tslearn"
    publishDir "${params.outdir}/report/", mode: 'copy'
    cache 'lenient'

    input:
    val cluster_ordering
    path output_path

    output:
    path "${cluster_ordering}_report.ipynb", emit: report_ipynb
    path "${cluster_ordering}_report.html", emit: report_html
    path "${cluster_ordering}_enrichR.csv", emit: enrichR_csv

    script:
    """
    export HOME=\$PWD
    export SCIPY_ARRAY_API=1
    /usersoftware/chanj3/tslearn/bin/python3 -m ipykernel install --user --name tslearn --display-name "[conda] tslearn"
    /usersoftware/chanj3/tslearn/bin/papermill ${baseDir}/bin/report.ipynb ${cluster_ordering}_report.ipynb \
    -p summary_df_path ${output_path}/${cluster_ordering}_summary_df.csv \
    -p aggregated_curves_path ${output_path}/${cluster_ordering}_aggregated_curves.csv \
    -p scores_df_path ${output_path}/${cluster_ordering}_scores_df.csv \
    -p hvg_df_path ${output_path}/hvg_genes_${cluster_ordering}.txt \
    -p cluster_ordering ${cluster_ordering}
    /usersoftware/chanj3/tslearn/bin/jupyter nbconvert --to html ${cluster_ordering}_report.ipynb
    """
}
