process INGEST {
    tag "${release}"
    label 'process_high_memory'

    input:
    tuple val(release), path(gctx), path(sig_info), path(gene_info), path(pert_info), path(sig_metrics), path(drug_info)

    output:
    path 'dataset.zarr',          emit: adata
    path 'ingest_manifest.json',  emit: manifest
    path 'versions.yml',          emit: versions

    script:
    def metrics = sig_metrics ? "--sig_metrics ${sig_metrics}" : ''
    def drugs = drug_info ? "--drug_info ${drug_info}" : ''
    def args = task.ext.args ?: ''
    """
    lincs-pipeline ingest \\
        --release ${release} \\
        --gctx ${gctx} \\
        --sig_info ${sig_info} \\
        --gene_info ${gene_info} \\
        --pert_info ${pert_info} \\
        ${metrics} ${drugs} ${args} \\
        --out dataset.zarr \\
        --manifest ingest_manifest.json

    lincs-pipeline versions --process ${task.process} > versions.yml
    """
}
