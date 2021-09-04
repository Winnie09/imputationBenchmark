# Benchmark for single-cell RNA-seq imputation methods

### For the community

**NEWS** 

Since we have received many emails requesting processed data or imputed data that was used in this benchmark work (thanks for the interest!), we decided to submit the data as an ExperimentData package to [Bioconductor](https://bioconductor.org/packages/3.12/BiocViews.html#___ExperimentData). We appreciate your patience and will update you as soon as it is released. For now, please feel free to access the following sharing Dropbox folder to get the data used in the codes.

**[sharing Dropbox folder for data](https://www.dropbox.com/sh/w3yg2nucnng5v1u/AAAM8Ym_KU9XF4z51RT81eNEa?dl=0). It consists of three subfolders: 

* bulkrna: bulk RNA-seq files.
* **processed** (most readers are interested in this folder): the single-cell RNA-seq data we processed and passed into imputation methods. Each subfolder contains multiple files. These are essentially the gene expression matrix for the same data, except that they are in different formats to meet the input requirements of different imputation methods (as mentioned in our manuscript). For example, "genebycell" means the each row represents a gene and each column represents a cell; "norm_" means the gene expression values have been normalized by library size. 
* raw: raw data of pbmc and hca bone marrow data. 

**[sharing Dropbox folder for realDE/result] (https://www.dropbox.com/sh/zztcqfnx9nhth4u/AAAx14KLezX6dBEzR9M8aKXZa?dl=0). It consists of three subfolders:

* cellbench
* GSE81861
* hca

We also received many emails from researchers who want to apply
this evaluation pipeline in some projects, or introduce these analyses in courses. If you have questions about the codes location of any analyses, or want to 
have a copy of the processed data in this work, please email me at whou10@jhu.edu.  We are more than happy to provide these. We appreicate your interest and patience.

### Publication

**Hou, W.**, Ji, Z., Ji, H.\* and Hicks, S.C.\*, (2020). *A Systematic Evaluation of Single-cell RNA-sequencing Imputation Methods*. 
Genome Biology 21, 218 (2020), doi: 10.1186/s13059-020-02132-x.  Links to: [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x), 
[Code](https://github.com/Winnie09/imputationBenchmark), [Twitter](https://twitter.com/GenomeBiology/status/1298976169484681219).
