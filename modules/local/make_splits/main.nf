process MAKE_SPLITS {
    tag "${split_type}"
    label 'process_single'

    input:
    path adata
    each split_type

    output:
    tuple val(split_type), path("${split_type}.parquet"), path("${split_type}.json"), emit: split
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    lincs-pipeline split \\
        --adata ${adata} \\
        --split_type ${split_type} \\
        ${args} \\
        --out_parquet ${split_type}.parquet \\
        --out_json ${split_type}.json

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
