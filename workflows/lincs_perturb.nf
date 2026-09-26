/*
 * Ingest -> QC & filter -> Split -> Fit & train -> Evaluate -> Package -> Gate
 *
 * Training fans out over a channel of (split, seed) tuples, one task each,
 * and the results are collected into one table across splits.
 */

include { INGEST          } from '../modules/local/ingest'
include { QC_FILTER       } from '../modules/local/qc_filter'
include { MAKE_SPLITS     } from '../modules/local/make_splits'
include { TRAIN           } from '../modules/local/train'
include { EVALUATE        } from '../modules/local/evaluate'
include { COLLECT_METRICS } from '../modules/local/collect_metrics'
include { PACKAGE_MODEL   } from '../modules/local/package_model'
include { REGISTRY_GATE   } from '../modules/local/registry_gate'

workflow LINCS_PERTURB {
    take:
    ch_raw          // [release, gctx, sig_info, gene_info, pert_info, sig_metrics|[], drug_info|[]]
    split_types     // list of split types
    seeds           // list of training seeds
    cell_features   // path or []

    main:
    INGEST(ch_raw)
    QC_FILTER(INGEST.out.adata)
    // One dataset version feeds every downstream task: make it a value channel.
    ch_adata = QC_FILTER.out.adata.first()
    MAKE_SPLITS(ch_adata, split_types)

    // Fan-out: every split file x every seed is one training task.
    ch_jobs = MAKE_SPLITS.out.split.combine(channel.fromList(seeds))
    TRAIN(ch_jobs, ch_adata, cell_features)

    ch_eval = TRAIN.out.run
        .map { split_type, seed, _model, preds, summary -> [split_type, seed, preds, summary] }
        .combine(MAKE_SPLITS.out.split.map { split_type, split, _json -> [split_type, split] }, by: 0)
    EVALUATE(ch_eval, ch_adata)

    COLLECT_METRICS(EVALUATE.out.metrics.map { _split_type, _seed, metrics -> metrics }.collect())

    // Package: every seed's model and metrics of one split, grouped together.
    ch_models = TRAIN.out.run
        .map { split_type, seed, model, _preds, _summary -> [split_type, model] }
        .groupTuple()
    ch_metrics = EVALUATE.out.metrics
        .map { split_type, _seed, metrics -> [split_type, metrics] }
        .groupTuple()
    ch_package = MAKE_SPLITS.out.split
        .map { split_type, _split, json -> [split_type, json] }
        .join(ch_models)
        .join(ch_metrics)
    PACKAGE_MODEL(ch_package)

    REGISTRY_GATE(COLLECT_METRICS.out.summary, PACKAGE_MODEL.out.bundle.collect())

    ch_versions = INGEST.out.versions
        .mix(QC_FILTER.out.versions, MAKE_SPLITS.out.versions, TRAIN.out.versions)
        .mix(EVALUATE.out.versions, COLLECT_METRICS.out.versions)
        .mix(PACKAGE_MODEL.out.versions, REGISTRY_GATE.out.versions)
        .map { yml -> yml.text }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true, newLine: true)

    emit:
    metrics  = COLLECT_METRICS.out.table
    summary  = COLLECT_METRICS.out.summary
    packages = PACKAGE_MODEL.out.bundle
    decision = REGISTRY_GATE.out.decision
    versions = ch_versions
}
