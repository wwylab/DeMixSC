% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeMixSC.deconvolution.R
\name{DeMixSC.core}
\alias{DeMixSC.core}
\title{The main functions of DeMixSC}
\usage{
DeMixSC.core(
  Y,
  X,
  iter.max = iter.max,
  eps = eps,
  print.sample = print.sample,
  output.file = output.file
)
}
\arguments{
\item{Y}{A normalized expression vector for the sample.}

\item{X}{A reference expression matrix, with genes on the rows and cell types on the columns.}

\item{iter.max}{Maximum number of iterations for the algorithm. Default is 2000.}

\item{eps}{Convergence threshold. The iterations stop when the change between iterations is below this threshold.
Default is 0.001.}

\item{print.sample}{Sample index for which iterative proportions should 
be recorded and written to the output file. This can be useful for 
monitoring the performance on a specific sample. Default is sample 1.}

\item{output.file}{Name of the output file where iterative data for 
the selected sample will be written.}
}
\value{
\item{proportion}{A cell-type proportion vector for all cell types in the sample.}
 \item{weights}{A weight vector for all genes in the sample.}
 \item{converge}{Number of iterations to reach convergence or a message indicating maximum iterations reached.}
}
\description{
This function estimates the cell type proportions for each bulk sample.
The function performs the DeMixSC algorithm by iteratively estimating cell-type proportions and updating weights.
}
\author{
Shuai Guo, Xiaoqian Liu, Wenyi Wang
}
