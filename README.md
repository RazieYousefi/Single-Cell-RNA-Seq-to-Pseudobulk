# GSE193531 single-cohort pipeline

A minimal, executable, single-cohort extract of one cohort's analysis pipeline from the
multiple myeloma disease-continuum signature project. The cohort was selected at random
from the project's three cohorts (GSE193531, GSE124310, GSE271107) for this release.

This repository covers one cohort's **within-cohort** pipeline only: raw count matrix to
cell-level QC, normalization/clustering, cell-type annotation, plasma-cell identity
confirmation, malignant/normal classification, and the patient-level Healthy-vs-MM
pseudobulk differential-expression test. It stops there. It does **not** include the
project's downstream, cross-cohort steps (concordance analysis across cohorts, the fixed
51-gene candidate signature, continuum/stage characterization, or pathway attribution) --
those combine this cohort with the other two and are out of scope for a single-cohort
release.

## Layout

The six notebooks live directly in this folder (`gitHub/`) and expect `data/` and
`results/` as **siblings of this folder** -- i.e. one level up, matching how the main
project is laid out (`.../MM_Progression/{gitHub, data, results, ...}`). Every notebook's
second code cell sets `PROJECT = Path("../")`, so it resolves to that project root and
reads/writes `PROJECT / "data" / ...` and `PROJECT / "results" / ...` there, not inside
`gitHub/` itself. If you move this folder out on its own (e.g. publishing it as a
standalone repository), create `data/` and `results/` as its siblings, matching this same
layout.

## Data

GSE193531 is distributed on GEO as a single dense UMI count matrix (genes x cells) plus a
companion cell-level metadata CSV, rather than per-sample CellRanger archives. Download
the cohort's files from GEO accession **GSE193531** and place them at (one level up from
this folder):

```
data/raw/GSE193531/GSE193531_umi-count-matrix.csv
data/raw/GSE193531/GSE193531_cell-level-metadata.csv
```

`GSE193531_umi-count-matrix.csv` has genes as rows and cells as columns, with count
values written as floats (`"0.0"`, not `"0"`). `GSE193531_cell-level-metadata.csv` must
contain, at minimum, an `index` column (cell ID matching the count matrix's column
names), `sample_ID`, `disease_stage` (values `NBM`/`MGUS`/`SMM`/`MM`), and
`normal_or_neoplastic` (the original study's own per-cell malignant/normal call, used by
notebook 5).

Processed `.h5ad` objects produced by each notebook land under `data/processed/GSE193531/`
(one level up) and are not checked into this repository -- they range from roughly 300 MB
to 700 MB for this cohort and regenerate in minutes to tens of minutes when the notebooks
are run.

## Pipeline

Run in order (each notebook reads the previous one's output):

| Step | Notebook | Produces |
|---|---|---|
| 1 | `01_qc.ipynb` | `data/processed/GSE193531/GSE193531_filtered.h5ad` -- hybrid (per-sample adaptive + fixed-floor) QC-filtered cells, doublets removed (Scrublet) |
| 2 | `02_normalization_clustering.ipynb` | `..._normalized.h5ad` -- library-size-normalized + log1p full gene set, HVGs flagged, PCA/UMAP/Leiden clustering |
| 3 | `03_cell_type_annotation.ipynb` | `..._annotated.h5ad` -- per-cluster cell-type calls (8 canonical marker sets, DE-cross-checked) |
| 4 | `04_plasma_cell_identity.ipynb` | `..._plasma_cells.h5ad` -- confirmed plasma-cell subset |
| 5 | `05_malignant_normal_classification.ipynb` | `..._plasma_cells_classified.h5ad` -- per-cell malignant/normal label (re-attached from the original study, kappa/lambda-checked) |
| 6 | `06_pseudobulk_deseq2.ipynb` | `results/GSE193531/pseudobulk_deg/GSE193531_healthy_vs_mm_deseq2_results.csv` -- patient-level Healthy-vs-MM pseudobulk DEG table |

Open each in Jupyter (or run with `jupyter nbconvert --to notebook --execute <name>.ipynb`)
and run top to bottom, in order, from a kernel whose working directory is this `gitHub/`
folder. Each notebook creates its own output directories under `data/processed/GSE193531/`
and `results/GSE193531/` and raises a clear `FileNotFoundError` if the previous notebook's
output is missing.

### Two steps involve a documented manual-review decision, not a mechanical computation

Leiden cluster numbering is not reproducible across runs (it depends on the exact
software versions and random seed in your environment). Two notebooks in this pipeline
make a judgment call about specific cluster numbers, and both ship with a conservative,
clearly-marked default that you should review and edit against your own run's printed
evidence rather than trust blindly:

- **`02_normalization_clustering.ipynb`** -- flagging non-plasma contaminant / QC-caution
  clusters (`MYELOID_CLUSTERS`, `B_CELL_CLUSTERS`, `QC_CAUTION_CLUSTERS`, in the "Cluster
  Classification" cell). Default: nothing flagged. This label (`cell_type_candidate`) is a
  supplementary sanity check used by notebook 3, not a hard input to any later computation.
- **`04_plasma_cell_identity.ipynb`** -- assigning each Leiden cluster as Plasma-like /
  Non-plasma / Ambiguous (`PLASMA_LIKE_CLUSTERS` / `NON_PLASMA_CLUSTERS` /
  `AMBIGUOUS_CLUSTERS`). Default: clusters already labeled `Plasma_cell` by notebook 3 are
  plasma-like; everything else is Ambiguous (excluded). This one **does** determine which
  cells reach every later step, so review the printed per-cluster marker evidence (tables
  and dotplot, in the notebook) before trusting it.

## R / DESeq2 setup (notebook 6 only)

`06_pseudobulk_deseq2.ipynb` calls Bioconductor's DESeq2 from Python via `rpy2`. This
needs, separately from `requirements.txt`:

1. A working R installation.
2. Bioconductor's `DESeq2` package:
   ```r
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("DESeq2")
   ```
3. Optionally, `apeglm` (used for log-fold-change shrinkage if present; the notebook falls
   back to the unshrunken DESeq2 estimate if it isn't installed):
   ```r
   BiocManager::install("apeglm")
   ```
4. `pip install rpy2` (already in `requirements.txt`), with `rpy2` able to find your R
   installation (usually automatic if `R` is on your `PATH`).

**Note on notebook 6's testing status:** this notebook is a direct port of the project's
own DESeq2 notebook. Notebooks 1-5 of this pipeline were executed end-to-end, top to
bottom, against a synthetic dataset while preparing this repository (confirming the
`data`/`results` path layout and every step's logic) and ran without error. Notebook 6's
pseudobulk construction and gene-filtering cells were verified the same way, chained
directly off notebooks 1-5's real output. The DESeq2/`rpy2` cell itself could not be
executed in the environment this repository was prepared in (no R installation available
there) -- the original project notebook carries the same caveat for the same reason. Run
it once carefully in your own environment and check the output against DESeq2's documented
interface before trusting it for anything downstream.

## Method summary

Cell-level QC combines per-sample adaptive thresholds (median minus 3 median absolute
deviations, in log1p space, for UMI counts and detected genes) with fixed guardrails
(minimum 500 UMIs, minimum 200 detected genes, maximum 15% mitochondrial transcript
percentage) and per-sample Scrublet doublet removal. Normalization is library-size
(target sum 1e4) plus log1p. Cell types are called by scoring 8 canonical marker sets per
cell (`sc.tl.score_genes`), aggregating to the Leiden-cluster level, and cross-checking
against each cluster's own unbiased differentially expressed genes. Plasma-cell identity
is confirmed at the cluster level against plasma and non-plasma lineage markers (this
cohort is CD138+-sorted, so cluster review -- not an automatic marker-threshold formula --
resolves genuinely ambiguous ambient-RNA-vs-mislabeling cases). Malignant/normal status
uses the original study's own per-cell label (re-attached from the raw GEO metadata, since
it was not carried into any processed object upstream), with kappa/lambda immunoglobulin
light-chain restriction computed as an independent, sample-level QC check -- never as a
second vote that overrides the primary label. Differential expression is patient-level
pseudobulk (raw UMI counts summed per patient, avoiding cell-level pseudoreplication) via
DESeq2, comparing each Healthy patient's `Normal`-labeled cells against each MM patient's
`Neoplastic`-labeled cells.

## License / provenance

Code extracted and adapted, for minimal release, from the analysis pipeline behind
"A multicompartment single-cell molecular signature associated with the multiple myeloma
disease continuum." Raw sequencing data is third-party (GEO accession GSE193531) and is
not redistributed here.
# Single-Cell-RNA-Seq-to-Pseudobulk
