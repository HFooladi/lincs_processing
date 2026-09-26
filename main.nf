#!/usr/bin/env nextflow
/*
 * LINCS L1000 perturbation-response prediction.
 * See README.md ("Perturbation-prediction workflow") for usage.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { LINCS_PERTURB } from './workflows/lincs_perturb'

def optionalFile(value) {
    return value ? file(value, checkIfExists: true) : []
}

def csvList(value) {
    return value.toString().tokenize(',')*.trim()
}

workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)

    if (!params.gctx || !params.sig_info || !params.gene_info || !params.pert_info) {
        error "Please provide --gctx, --sig_info, --gene_info and --pert_info (or use -profile test)"
    }
    if (params.release != 'beta' && !params.sig_metrics) {
        error "--release ${params.release} needs --sig_metrics"
    }
    if (params.tracker == 'mlflow') {
        // SQLite needs its directory to exist before parallel tasks open it.
        file(params.mlflow_dir).mkdirs()
    }

    ch_raw = channel.of([
        params.release,
        file(params.gctx, checkIfExists: true),
        file(params.sig_info, checkIfExists: true),
        file(params.gene_info, checkIfExists: true),
        file(params.pert_info, checkIfExists: true),
        optionalFile(params.sig_metrics),
        optionalFile(params.drug_info),
    ])

    LINCS_PERTURB(
        ch_raw,
        csvList(params.splits),
        csvList(params.seeds).collect { seed -> seed.toInteger() },
        optionalFile(params.cell_features),
    )
}
