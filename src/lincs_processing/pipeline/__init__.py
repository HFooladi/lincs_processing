"""Perturbation-prediction pipeline steps driven by the Nextflow workflow.

Each step is a module with a plain Python API plus an entry in
``lincs_processing.pipeline.cli`` so ``lincs-pipeline <step>`` can call it from a
Nextflow process. Heavy dependencies (anndata, torch, rdkit, mlflow) come from
the ``pipeline`` dependency group and are imported lazily by the CLI.
"""
