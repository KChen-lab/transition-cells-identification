# Identifying and characterizing transition cells in developmental processes from scRNA-seq data 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14027784.svg)](https://doi.org/10.5281/zenodo.14027784)

Decoding cellular state transitions is crucial for understanding complex biological processes in development and disease. While recent advancements in single-cell RNA sequencing offer insights into cellular trajectories, existing tools primarily study expressional rather than regulatory state shifts. We present a novel framework utilizing gene-pair correlations to detect transition cells from scRNA-seq data. Applying our approach to various contexts, including tissue regeneration, preinvasive lesions and humoral responses post-vaccination, reveals transition cells and their distinct gene expression profiles. Our study sheds light onto the underlying molecular mechanisms driving cellular state transitions enhancing our ability to identify therapeutic targets for disease interventions.

Our method models gene expressions using stochastic differential equations (SDEs), and identifies transition cells through transition index calculated based on gene pair-wise correlation coefficients: <br />

![workflow](https://github.com/KChen-lab/transition-cells-identification/blob/main/images/workflow.png)

Through characterizing cellular state transitions, we can identify genes potentially useful for diagnosis, prognosis and therapeutics.

For more informations and mathematical details please refer to our manuscript

# Installation
To install the developmental version from GitHub:

```
if(!require(devtools)) install.packages("devtools");
devtools::install_github("KChen-lab/transition-cells-identification",force=T)
```
To load the installed package in R:
```
library(CellTran)
```
# Usage
[Here](https://github.com/KChen-lab/transition-cells-identification/blob/main/example/identify_transition_cells_using_simulation_data.ipynb) is an example about how to apply our method in scRNA-seq data. The data used in the example is simulated by SERGIO with true stable and transition labels and can be found in the [data/](https://github.com/KChen-lab/transition-cells-identification/tree/main/data) folder.
## Inputs
```data```: A Seurat object containing both normalized count matrix and metadata such as neighboring results for a single-cell dataset.<br/>
```highly_variable_gene```: A gene list containing top most variable genes. Calculated by ```var``` function by default. Can be set manually if particular genes need to be included or excluded. E.g. cell cycle-related genes can be excluded if focusing on non-cell cycle transitions.<br/>
```group```: The cell group highly_variable_genes selected based on. To capture different expression patterns in different potential stable states, we calculate GPPCCs using group specific top most variable genes. The group should be a column name of `data@meta.data`, such as 'seurat_clusters','cell_type','time_point'. Using 'seurat_clusters' by default. <br/>
```n_neighbor```: The number of neighboring cells used to calculate GPPCCs. This number choice is empirically determined from the data, to achieve a good tradeoff between temporal resolution and estimation accuracy. Usually 200~300 neighbors achieve a good performance. Using 300 by default.<br/>
```n_gene```: The number of top most variable genes used in calculating GPPCCs. Using 50 by default.<br/>
```return_pearson```: Whether report GPPCCs or not. Using False by default.<br/>
## Multiple datasets comparison
The transition index calculated using this method reflects the relative likelihood of a cell to be a transition cell. We assume there are both transition cells and stable cells in the dataset. We first find the archetype of a transition cell and a stable cell, and calculate the transition indices for all other cells based on these two archetypal cells. When comparing cell transition indices from multiple datasets, users can either first combine cells together and then calculate transition indices (recommended) or include the same archetypal stable and transition cells for all the datasets.      
