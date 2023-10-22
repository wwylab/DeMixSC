# DeMixSC for bulk RNA-seq deconvolution

DeMixSC is a generalizable framework designed to leverage single-cell sequencing and a small benchmark dataset for bulk RNA-seq deconvolution.
Please refer to our [preprint](https://www.biorxiv.org/content/10.1101/2023.10.10.561733) for details of DeMixSC.

# Framework overview
DeMixSC offers accurate cell-type deconvolution for large bulk RNA-seq datasets through a two-tier procedure.  

First, DeMixSC utilizes a specifically designed benchmark dataset to identify and adjust genes with high inter-platform discrepancies (Fig.1A). 

Second, to deconvolve a large unmatched bulk RNA-seq dataset, DeMixSC aligns it with the benchmark dataset (Fig.1B) and then employs a refined weighted non-negative least square (wNNLS) framework for deconvolution (Fig.1C).

Note: DeMixSC requires a matched tissue type between the small benchmark dataset and the large targeted cohort.

<figure>
  <img src="./figures/framework.png" width="800px"/>
  <figcaption>Rescource: Figure 1, Workflow of DeMixSC (Guo et al. 2023, doi: https://doi.org/10.1101/2023.10.10.561733).</figcaption>
</figure>

# Usage

DeMixSC is implemented in R. To use DeMixSC, please install the required packages as detailed below: 

```r
# Install BiocManager if necessary
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")     
      
# Install required pacakges
  BiocManager::install("sva")
  BiocManager::install("preprocessCore")
  install.packages("doParallel")
  install.packages("nnls")
  install.packages('Seurat')
```

Next, please install the latest version of DeMixSC from our GitHub repository:

```r
# install devtools if necessary
  if (!"devtools" %in% rownames(installed.packages())) {
    install.packages('devtools')
  }

# Install our DeMixSC package with the following command
  devtools::install_github('wwylab/DeMixSC')
```

You can then use the provided functions to conduct deconvolution analyses. We provide a [tutorial](https://sphingosine.github.io/tutorial/DeMixSC.html) on using DeMixSC for deconvolving the benchmarking dataset. 


# License
This package is licensed under GNU Affero General Public License v3.0. See the LICENSE file for details. 

# Contact
For questions or comments about DeMixSC, please contact the package maintainers at Shuai Guo <SGuo3@mdanderson.org> and Xiaoqian Liu <XLiu31@mdanderson.org>. If you find a bug or have a feature request, please submit an issue on the GitHub repository.



