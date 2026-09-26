process TRAIN {
    tag "${split_type}:seed${seed}"
    label 'process_gpu'

    input:
    tuple val(split_type), path(split), path(split_json), val(seed)
    path adata
    path cell_features

    output:
    tuple val(split_type), val(seed), path("model_seed${seed}"), path("predictions_seed${seed}.npz"), path("train_summary_seed${seed}.json"), emit: run
    path 'versions.yml', emit: versions

    script:
    def cells = cell_features ? "--cell_features ${cell_features}" : ''
    def args = task.ext.args ?: ''
    """
    lincs-pipeline train \\
        --adata ${adata} \\
        --split ${split} \\
        --split_json ${split_json} \\
        --seed ${seed} \\
        ${cells} ${args} \\
        --outdir .

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
