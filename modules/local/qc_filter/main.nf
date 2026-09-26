process QC_FILTER {
    label 'process_medium'

    input:
    path adata

    output:
    path 'filtered.zarr',   emit: adata
    path 'qc_report.json',  emit: report
    path 'versions.yml',    emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    lincs-pipeline qc --adata ${adata} ${args} --out filtered.zarr --report qc_report.json

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
