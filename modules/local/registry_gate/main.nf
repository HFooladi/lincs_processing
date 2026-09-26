process REGISTRY_GATE {
    label 'process_single'

    input:
    path summary
    path packages, stageAs: 'packages/*'

    output:
    path 'gate_decision.json', emit: decision
    path 'versions.yml',       emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    lincs-pipeline gate --summary ${summary} --packages packages ${args} --out gate_decision.json

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
