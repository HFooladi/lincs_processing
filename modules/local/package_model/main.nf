process PACKAGE_MODEL {
    tag "${split_type}"
    label 'process_single'

    input:
    tuple val(split_type), path(split_json), path(models), path(metrics)

    output:
    path 'packages/*',   emit: bundle
    path 'versions.yml', emit: versions

    script:
    """
    lincs-pipeline package \\
        --split_json ${split_json} \\
        --models ${models} \\
        --metrics ${metrics} \\
        --outdir packages

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
