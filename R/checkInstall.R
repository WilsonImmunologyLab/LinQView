if (!require("Seurat")) install.packages("Seurat")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("cowplot")) install.packages("cowplot")
if (!require("FastKNN")) install.packages("FastKNN")
if (!require("reshape2")) install.packages("reshape2")
if (!require("umap")) install.packages("umap")
if (!require("emstreeR")) install.packages("emstreeR")
if (!require("scales")) install.packages("scales")
if (!require("Rtsne")) install.packages("Rtsne")
if (!require("grDevices")) install.packages("grDevices")
if (!require("Rcpp")) install.packages("Rcpp")
if (!require("igraph")) install.packages("igraph")
if (!require("ElPiGraph.R")) {
  if (!require("devtools")) {install.packages("devtools")}
  devtools::install_github("Albluca/distutils")
  devtools::install_github("Albluca/ElPiGraph.R")
}
