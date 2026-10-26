# GSE193531 single-cohort pipeline

An executable, single-cohort extract of analysis pipeline from the
multiple myeloma disease-continuum signature project (GSE193531).

This repository covers: raw count matrix to
cell-level QC, normalization/clustering, cell-type annotation, plasma-cell identity
confirmation, malignant/normal classification, and the patient-level Healthy-vs-MM
pseudobulk differential-expression test.

## Layout

The six notebooks live directly here.

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
values written as floats (`"0.0"`). `GSE193531_cell-level-metadata.csv` must
contain, at minimum, an `index` column (cell ID matching the count matrix's column
names), `sample_ID`, `disease_stage` (values `NBM`/`MGUS`/`SMM`/`MM`), and
`normal_or_neoplastic` (the original study's own per-cell malignant/normal call, used by
notebook 5).

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

Open each in Jupyter and run top to bottom, in order, from a kernel whose working directory is this `gitHub/`
folder. Each notebook creates its own output directories under `data/processed/GSE193531/`
and `results/GSE193531/` and raises a clear `FileNotFoundError` if the previous notebook's
output is missing.


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
disease continuum." Raw sequencing data is GEO accession GSE193531.

