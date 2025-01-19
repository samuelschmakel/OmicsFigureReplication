# Omics Figure Replication
Analysis of RNA-seq dataset from pre-print. Paper can be found [here](https://pubmed.ncbi.nlm.nih.gov/39677760/). 

## Motivation
This project was undertaken to demonstrate my skills in bioinformatics and data visualization. The goal was to replicate graphics used in a recently submitted paper, including PCA plots, heatmaps for pairwise DEG comparisons, and Gene Ontology (GO) visualizations. These visualizaitons are critical tools for interpreting high-throughput omics data and drawing biologically meaningful conclusions.

## Installation
### Install R
Ensure that R is installed on your system. You can download the latest version of R from the [CRAN website](https://www.r-project.org/).

### Clone the repository
Clone the repository to your local machine with the following command:
```bash
   git clone https://github.com/samuelschmakel/OmicsFigureReplication.git
   cd OmicsFigureReplication
```

### Install Required R Packages
Install the required R packages by running the following commands in R:
```R
  install.packages(c("ggplot2", "pheatmap", "GEOquery", "DESeq2", "clusterProfiler", "org.Hs.eg.db"))
```

## Usage
Open `analysis.R` in RStudio (or your preferred R environment) and execute the code.
