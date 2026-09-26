process COLLECT_METRICS {
    label 'process_single'

    input:
    path metrics

    output:
    path 'metrics.tsv',   emit: table
    path 'summary.tsv',   emit: summary
    path 'versions.yml',  emit: versions

    script:
    """
    lincs-pipeline collect --metrics ${metrics} --out metrics.tsv --summary summary.tsv

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
