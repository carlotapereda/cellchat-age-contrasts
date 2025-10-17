# CellChat Age-Contrast Analyses
Author: Carlota Pereda Serras
Last Edited: 10/17/2025

Scripts and precomputed CellChat objects to compute pathway- and LR-level age contrasts
across sex (F/M) and APOE genotype (E33/E44).

![Workflow](./workflow.png)

## Contents
- `createCCObjs.R` — build CellChat objects (DB v2) by sex/age/genotype.
- `netDiffInteraction_generateGraphs.R` — differential network visuals (Panel A).
- `heatmap_pathways.R` — pathway-level log2FC matrices (Panel B).
- `heatmap_LR.R` — GABA-A/B LR log2FC matrices for inhibitory neurons (Panel C).
- `netDiffinteraction_adapted.R` — adapted plotting function used by Panel A.

## Running
Use `Rscript <script>.R` from this directory. Scripts log to corresponding `panel*/...*.log`.

## Data
`.rds` objects are tracked via Git LFS.



