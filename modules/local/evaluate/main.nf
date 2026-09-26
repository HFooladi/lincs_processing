process EVALUATE {
    tag "${split_type}:seed${seed}"
    label 'process_single'

    input:
    tuple val(split_type), val(seed), path(predictions), path(train_summary), path(split)
    path adata

    output:
    tuple val(split_type), val(seed), path("metrics_${split_type}_seed${seed}.json"), emit: metrics
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    lincs-pipeline evaluate \\
        --adata ${adata} \\
        --split ${split} \\
        --predictions ${predictions} \\
        --train_summary ${train_summary} \\
        ${args} \\
        --out metrics_${split_type}_seed${seed}.json

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
