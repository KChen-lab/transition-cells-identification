% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transition_index.R
\name{gene_pearson}
\alias{gene_pearson}
\title{Calculate GPPCCs}
\usage{
gene_pearson(
  data,
  highly_variable_gene = NULL,
  n_neighbor = 300,
  n_gene = 50,
  ...
)
}
\arguments{
\item{data}{Seurat object containing both count matrix and metadata such as neighboring
results for a single-cell dataset.}

\item{highly_variable_gene}{A gene list containing top most variable genes. Calculated
by \code{var} function by default.}

\item{n_neighbor}{The number of neighboring cells used to calculate GPPCCs. 300 by default.}
}
\value{
A matrix of GPPCCs. Rows are cells and columns are gene pairs.
}
\description{
This is the function calculating gene pair-wise Pearson's correlation coefficient (GPPCCs).
The GPPCCs calculated using our method are assumed from the cells that exhibit similar
regulatory profiles and around equilibrium points. Thus, for each cell, we select the cell
together with its nearest neighbors in Principal Component Analysis (PCA) dimensions based
on cosine distance to calculate GPPCCs. We also use scLink in the calculation to deal with
the gene expression values affected by dropout.
}
