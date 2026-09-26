# Raw LINCS files (DVC-tracked, not in git)

Put one release per folder and track it with DVC so its hash is versioned:

```bash
# CLUE 2020 beta release (https://clue.io/data/CMap2020#LINCS2020)
data/raw/beta/level5_beta_trt_cp_n720216x12328.gctx
data/raw/beta/siginfo_beta.txt
data/raw/beta/geneinfo_beta.txt
data/raw/beta/compoundinfo_beta.txt

uv run dvc add data/raw/beta
git add data/raw/beta.dvc data/raw/.gitignore
```

For GEO releases use `data/raw/gse92742/` (or `gse70138/`) with the level-5
`.gctx`, `sig_info`, `sig_metrics`, `pert_info` and `gene_info` tables, and run
the workflow with `--release gse92742 --sig_metrics ...`.
