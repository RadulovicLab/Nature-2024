# Header
This is the R code for scRNAseq data analysis for paper "Formation of Memory Assemblies through the DNA Sensing Tlr9 Pathway". Inputs are five 10x scRNAseq samples, raw fastq and matrix data were uploaded to GEO: [GSE254780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254780). 
## System Requirements
### OS Requirements
This code has been tested on the following system:
* Item Ventura (13.3.1)
### R Dependencies
```
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(SeuratWrappers)
library(doBy)
library(fgsea)
library(data.table)
library(pheatmap)
library(gprofiler2)
library(scDblFinder)
library(SoupX)
library(lsa)
library(msigdbr)
```
The code was run with R version 4.3.2 (2023-10-31)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.3.1
