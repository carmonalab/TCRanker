# Ranking TCRs to predict tumor-reactivity by paired scRNA and scTCR data

`TCRanker` is a tool to prioritize the TCRs that are more likely to be tumor-reactive.

It takes as input scRNA-seq data paired with scTCR-seq data, filters cell type wanted (CD8+ T cells by default), and ranks TCRs in a tumor sample according to their likelihood to recognize tumor antigens. To rank TCRs, TCRanker evaluates T cell clonal expansion and T cell transcriptional features associated with tumor reactivity, including exhaustion and proliferation level, in each T cell clonotype.

<br/>

<img src="./docs/TCRanker_workflow.svg">
</p>
<p align = "center" style="color:grey">
TCRanker workflow
</p>

<br/>

## Installation
`TCRanker`  has dependency on:

- R (>= 4.2.0)
- [UCell](https://github.com/carmonalab/UCell) (>= 2.0.0)
- [scGate](https://github.com/carmonalab/scGate) (>= 1.0.2)

To install `TCRanker`, run:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
remotes::install_github("carmonalab/scGate")
remotes::install_github("carmonalab/TCRanker")
```

<br/>

## Data preparation

Currently, `TCRanker` supports both `Seurat` and `SingleCellExperiment` object as the input query. Optimally, the TCR data should be stored in the query as a column of metadata.

If you want to know how to prepare the data set as input query, we provide you an demonstration of how to prepare a basic Seurat object containing expression matrix together with TCR meta data [here](https://carmonalab.github.io/TCRanker.demo/preparation.html)

<br/>

## TCRanker
The function has following parameters:
```
TCRanker(
  query,
  tcr,
  signature = "default",
  assay = NULL,
  group = "none",
  exhaustion = TRUE,
  proliferation = TRUE,
  species = "auto",
  FUN = "mean",
  minClonSize = 5,
  filterCell = "CD8T",
  strictFilter = TRUE,
  keepObject = FALSE,
  ...
)
```
Only `query` and `tcr` are indispensable.

Your could refer [TCRanker.demo](https://carmonalab.github.io/TCRanker.demo/demo.html) for more detailed illustration of parameters.

<br/>

## Output
The default output would be a data frame with following structure:

| clonotype | size | freq | exhaustion.score | exhaustion.ranking | proliferation.score | proliferation.ranking |
|:---------:|:----:|:----:|:----------------:|:------------------:|:-------------------:|:---------------------:|
|    ...    | ...  | ...  |        ...       |          ...       |       ...           |           ...         |

If `group` was provided, it'd be the column right after `clonotype`

If custom gene signature(s) were provided, there would be 2 extra columns for each signature. They would be named "SignatureName.score" and "SignatureName.ranking" and positioned prior to the exhaustion and proliferation columns.

<br/>

## Sequential analysis
In order to conduct further in vitro validation or upon the ranking result, you might need to assemble the full nucleotide sequence based on the output. You can check [here](https://carmonalab.github.io/TCRanker.demo/stitchr.html) to see how to use [stitchr](https://github.com/JamieHeather/stitchr) ([Heather et al., (2022)](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac190/6553689)) as an example.