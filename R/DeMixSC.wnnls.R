#' @title  The weighted non-negative least square (wNNLS) framework of DeMixSC
#'
#' @author Shuai Guo, Xiaoqian Liu, Wenyi Wang

############################################################

#' @title Estimate cell type proportions using wNNLS
#'
#' @description
#' This function estimates the cell type proportions for each bulk sample.
#' The function performs the DeMixSC algorithm by iteratively estimating cell-type proportions and updating weights.
#'
#'
#' @param j Index used for parallel processing.
#'
#' @param MAT A normalized expression matrix.
#'
#' @param R A cell-type-specific gene expression matrix, where rows are genes and columns are cell types.
#'
#' @param c A numeric value representing the tuning parameter for the weight function in the weighted non-negative least squares algorithm. Default is 2.
#'
#' @param iter.max An integer specifying the maximum number of iterations for the deconvolution algorithm.
#' The algorithm will stop if this number is reached, even if convergence (steady results over iterations) has not been achieved. Default is 2000.
#'
#' @param eps A numeric value representing the convergence threshold for the deconvolution algorithm.
#' The algorithm will stop if the difference between iterations falls below this threshold. Default is 0.001.
#'
#' @param print.sample An integer specifying the index of the sample whose iterative results should be printed during the deconvolution process.
#' This is useful for monitoring the algorithm’s progress on a specific sample. Default is 1.
#' Users can set `print.sample = NULL` to disable saving iteration details.
#'
#' @param output.file A character string specifying the name of the file where the algorithm’s iterative results will be saved.
#' This file stores intermediate deconvolution results for the specified sample. Default is `"./DeMixSC.save.iter.txt"`.
#'
#'
#' @return
#'  \item{proportion}{A cell-type proportion vector for all cell types in the sample.}
#'  \item{converge}{Number of iterations to reach convergence or a message indicating maximum iterations reached.}
#'
#'
#' @export

DeMixSC.wnnls = function(j, MAT, R, c, iter.max = iter.max, eps = eps, print.sample = print.sample, output.file = output.file){

  # Initialize the model using an unweighted matrix.
  Y = MAT[,j]
  k = ncol(R)
  est = nnls(R, Y)

  # Compute the initial weights for all genes; these weights are updated during each iteration.
  sqrt.weight = sqrt(1/((est$residuals)^2 + (est$fitted)^2 + c))

  # Apply weights to both the sample vector (Y) and the reference matrix (X).
  weighted.Y = Y * (sqrt.weight)
  weighted.R = R * as.matrix((sqrt.weight))[,rep(1,k)]
  weighted.est = nnls(weighted.R, weighted.Y)

  # Calculate the initial estimated cell type proportion.
  proportion = weighted.est$x/sum(weighted.est$x)

  # Record the initial estimates for the specified sample.
  if (!is.null(print.sample)) {
    if (print.sample == j) {
      message <- paste0(Sys.time()," - Processed sample: ",j," - the initial estimates: \n",
                        paste(colnames(R), collapse = ", "), "\n",
                        paste(round(proportion,4), collapse = ", "), "\n")
      write(message, file = output.file, append = TRUE)
    }
  }

  # Iterate to refine the estimates.
  for(iter in 1:iter.max){
    sqrt.weight = sqrt(1/( (weighted.est$residuals)^2 + (weighted.est$fitted)^2 + c ) )
    weighted.Y = Y * sqrt.weight
    weighted.R = R * as.matrix(sqrt.weight)[,rep(1,k)]
    weighted.est = nnls(weighted.R, weighted.Y)
    proportion.iter = weighted.est$x/sum(weighted.est$x)

    # Record the estimates for the specified sample during each iteration.
    if (!is.null(print.sample)) {
      if (print.sample == j) {
        message <- paste0(Sys.time()," - Processed sample ",j," - the estimates in iteration ",iter,":\n",
                          paste(colnames(R), collapse = ", "), "\n",
                          paste(round(proportion.iter,4), collapse = ", "), "\n")
        write(message, file = output.file, append = TRUE)
      }
    }

    # Check convergence; if reached, return the results.
    if(sum(abs(proportion.iter - proportion)) < k*eps){
      proportion = proportion.iter
      return(list(proportion = proportion,
                  converge = paste0('converge at ', iter)))
    }
    proportion = proportion.iter

  }
  return(list(proportion = proportion,
              converge = paste0('reach iter.max ', iter.max)))
}

#' @title Parallel execution of DeMixSC.wnnls for estimating cell type proportions
#'
#' @description
#' Executes the DeMixSC.wnnls function in parallel for estimating cell-type proportions.
#'
#'
#' @param bulk.mat A gene expression matrix to be deconvolved, where rows are genes and columns are samples.
#'
#' @param ref A cell-type-specific gene expression matrix, where rows are genes and columns are cell types.
#'
#' @param c A numeric value representing the tuning parameter for the weight function in the weighted non-negative least squares algorithm. Default is 2.
#'
#' @param iter.max An integer specifying the maximum number of iterations for the deconvolution algorithm.
#' The algorithm will stop if this number is reached, even if convergence (steady results over iterations) has not been achieved. Default is 2000.
#'
#' @param eps A numeric value representing the convergence threshold for the deconvolution algorithm.
#' The algorithm will stop if the difference between iterations falls below this threshold. Default is 0.001.
#'
#' @param nthread An integer specifying the number of threads to use for parallel processing (executing multiple operations simultaneously).
#' By default, it uses one less than the total number of available CPU cores.
#'
#' @param print.sample An integer specifying the index of the sample whose iterative results should be printed during the deconvolution process.
#' This is useful for monitoring the algorithm’s progress on a specific sample. Default is 1.
#' Users can set `print.sample = NULL` to disable saving iteration details.
#'
#' @param output.file A character string specifying the name of the file where the algorithm’s iterative results will be saved.
#' This file stores intermediate deconvolution results for the specified sample. Default is `"./DeMixSC.save.iter.txt"`.
#'
#'
#' @return
#'  \item{proportion}{A matrix of cell-type proportions, with cell types on the rows and samples on the columns.}
#'  \item{converge}{A vector indicating the number of iterations until convergence for each sample.}
#'
#'
#' @note
#' By default, this function utilizes the maximum number of available cores minus one for parallel computation.
#' This is designed to maximize computational efficiency.
#' However, users should be cautious when running other tasks concurrently to avoid potential system slowdowns or overloading.
#' If resource concerns arise or if you're operating on a shared system, users might consider specifying a lower number of threads through the `nthread` parameter.
#'
#' @export

DeMixSC.wnnls.parallel = function(bulk.mat,
                                  ref,
                                  c = 2,
                                  iter.max = 2000,
                                  eps = 0.001,
                                  nthread = (detectCores()[1]-1),
                                  print.sample = 1,
                                  output.file = "./DeMixSC.save.iter.txt"){

  # Re-ensure consistency in genes between reference and target bulk.
  common.gene = intersect(rownames(bulk.mat), rownames(ref))
  bulk.mat = bulk.mat[match(common.gene, rownames(bulk.mat)), ]
  ref = ref[match(common.gene, rownames(ref)), ]

  # Log information for iterative records (if enabled).
  if (!is.null(print.sample)) {
    cat("Note: DeMixSC typically converges in <600 iterations under the default convergence threshold (esp=0.001).\n")
    cat("For demonstration, proportions from each iteration of the first sample are recorded by default.\n")
    cat(paste0("Iterative data will be stored in the file: \"",output.file,"\".\n"))
    write("", file = output.file)
  }

  cat("Note: DeMixSC defaults to using (max.cores-1). Be cautious to avoid system overload. Consider adjusting `nthread` if needed.\n")

  # Set up the parallel back-end.
  cl = makeCluster(nthread)

  # Start parallel deconvolution.
  DeMixSC.res = NULL
  cat("Starting parallel computation... \n")
  registerDoParallel(cl)
  DeMixSC.res = pbapply::pblapply(X = 1:ncol(bulk.mat), FUN = DeMixSC.wnnls, cl = nthread,
                                  MAT = bulk.mat, R = ref, c = c, iter.max = iter.max, eps = eps, print.sample = print.sample, output.file = output.file)
  stopCluster(cl)
  cat("Parallel computation completed. \n")

  # Compile the estimated cell type proportions.
  combined.proportion = do.call(cbind, lapply(DeMixSC.res, function(x) x$proportion))
  rownames(combined.proportion) = colnames(ref)
  colnames(combined.proportion) = colnames(bulk.mat)

  # Compile the convergence information.
  combined.converge = unlist(lapply(DeMixSC.res, function(x) x$converge))
  names(combined.converge) = colnames(bulk.mat)

  # Deconvolution is completed and return results.
  cat("Deconvolution tasks completed! \n")
  return(list(cell.type.proportions = combined.proportion,
              converge = combined.converge))

}
