#' @title  The main functions of DeMixSC
#'
#' @author Shuai Guo, Xiaoqian Liu, Wenyi Wang

############################################################

#' @title Estimate cell type proportions using DeMixSC
#' 
#' @description
#' This function estimates the cell type proportions for each bulk sample.
#' The function performs the DeMixSC algorithm by iteratively estimating cell-type proportions and updating weights.
#' 
#' @param Y A normalized expression vector for the sample.
#' @param X A reference expression matrix, with genes on the rows and cell types on the columns.
#' @param iter.max Maximum number of iterations for the algorithm. Default is 2000.
#' @param eps Convergence threshold. The iterations stop when the change between iterations is below this threshold.
#' Default is 0.001.
#' @param print.sample Sample index for which iterative proportions should 
#' be recorded and written to the output file. This can be useful for 
#' monitoring the performance on a specific sample. Default is sample 1.
#' @param output.file Name of the output file where iterative data for 
#' the selected sample will be written.
#' 
#' @return
#'  \item{proportion}{A cell-type proportion vector for all cell types in the sample.}
#'  \item{weights}{A weight vector for all genes in the sample.}
#'  \item{converge}{Number of iterations to reach convergence or a message indicating maximum iterations reached.}
#'
#' @export

DeMixSC.core = function(Y, X, iter.max = iter.max, eps = eps, print.sample = print.sample, output.file = output.file){ 
  
  # Initialize the model using an unweighted matrix.
  k = ncol(X)
  est = nnls(X, Y)
  
  # Compute the initial weights for all genes; these weights are updated during each iteration.
  sqrt.weight = sqrt(1/((est$residuals)^2 + (est$fitted)^2 + 2))
  
  # Apply weights to both the sample vector (Y) and the reference matrix (X).
  weighted.Y = Y * (sqrt.weight)
  weighted.X = X * as.matrix((sqrt.weight))[,rep(1,k)]
  weighted.est = nnls(weighted.X, weighted.Y)
  
  # Calculate the initial estimated cell type proportion.
  proportion = weighted.est$x/sum(weighted.est$x)
  
  # Record the initial estimates for the specified sample.
  if (!is.null(print.sample)) {
    if (print.sample == j) {
      message <- paste0(Sys.time()," - Processed sample: ",j," - the initial estimates: \n", 
                        paste(colnames(X), collapse = ", "), "\n",
                        paste(round(proportion,4), collapse = ", "), "\n")
      write(message, file = output.file, append = TRUE)
    }
  }
  
  # Iterate to refine the estimates.
  for(iter in 1:iter.max){
    sqrt.weight = sqrt(1/( (weighted.est$residuals)^2 + (weighted.est$fitted)^2 + 2 ) )
    weighted.Y = Y * sqrt.weight
    weighted.X = X * as.matrix(sqrt.weight)[,rep(1,k)]
    weighted.est = nnls(weighted.X, weighted.Y)
    proportion.iter = weighted.est$x/sum(weighted.est$x)
    
    # Record the estimates for the specified sample during each iteration
    if (!is.null(print.sample)) {
      if (print.sample == j) {
        message <- paste0(Sys.time()," - Processed sample ",j," - the estimates in iteration ",iter,":\n", 
                          paste(colnames(X), collapse = ", "), "\n",
                          paste(round(proportion.iter,4), collapse = ", "), "\n")
        write(message, file = output.file, append = TRUE)
      }
    }
    
    # Check convergence; if reached, return the results.
    if(sum(abs(proportion.iter - proportion)) < k*eps){
      proportion = proportion.iter
      sqrt.weight.final = round(sqrt.weight, 4)
      return(list(proportion = proportion,
                  weights = sqrt.weight.final,
                  converge = paste0('converge at ', iter)))
    }
    proportion = proportion.iter
    sqrt.weight.final = round(sqrt.weight, 4)
  }
  return(list(proportion = proportion,
              weights = sqrt.weight.final,
              converge = paste0('reach iter.max ', iter.max)))
}

#' @title Parallel execution of DeMixSC.core for deconvolution
#' 
#' @description
#' Executes the DeMixSC.core function in parallel for estimating cell-type proportions.
#' 
#' @param bulk.mat A expression matrix, which genes on the rows and samples on the columns.
#' @param ref Reference matrix of expression signatures, with genes on the rows and cell types on the columns.
#' @param iter.max Maximum number of iterations for the algorithm. Default is 2000.
#' @param eps Convergence threshold. The iterations stop when the change between iterations is below this threshold.
#' Default is 0.001.
#' @param nthread Number of threads allocated for deconvolution. Default is using one less than total available threads.
#' @param print.sample Sample index for which iterative proportions should 
#' be recorded and written to the output file. This can be useful for 
#' monitoring the performance on a specific sample. Default is sample 1.
#' @param output.file Name of the output file where iterative data for 
#' the selected sample will be written.
#'
#' @return 
#'  \item{proportion}{A matrix of cell-type proportions, with cell types on the rows and samples on the columns.}
#'  \item{converge}{A vector indicating the number of iterations until convergence for each sample.}
#'  \item{weights}{A matrix of gene weights with genes, weith genes on the rows and samples on the columns.}
#'  
#' @note 
#' By default, this function utilizes the maximum number of available cores minus one for parallel computation. 
#' This is designed to maximize computational efficiency. 
#' However, users should be cautious when running other tasks concurrently to avoid potential system slowdowns or overloading.
#' If resource concerns arise or if you're operating on a shared system, users might consider specifying a lower number of threads through the `nthread` parameter.
#' 
#' @export

DeMixSC.deconvolution = function(bulk.mat, ref, iter.max = 2000, eps = 0.001, nthread = (detectCores()[1]-1), 
                                 print.sample = 1, output.file = "./DeMixSC.save.iter.txt"){
  
  common.gene = intersect(rownames(bulk.mat), rownames(ref))
  bulk.mat = bulk.mat[match(common.gene, rownames(bulk.mat)), ]
  ref = ref[match(common.gene, rownames(ref)), ]
  DeMixSC.res = NULL
  
  message("Note: DeMixSC defaults to using (max.cores-1). Be cautious to avoid system overload. Consider adjusting `nthread` if needed.")
  
  # Set up the parallel backend.
  cl = makeCluster(nthread)
  
  # Log information for iterative records (if enabled).
  if (!is.null(print.sample)) {
    message("Note: DeMixSC typically converges in <600 iterations under the default convergence threshold (esp=0.001).\n")
    message("For demonstration, proportions from each iteration of the first sample are recorded by default.\n")
    message(paste0("Iterative data will be stored in the file: \"",output.file,"\".\n"))
    write("", file = output.file)
  }
  
  # Start parallel deconvolutions.
  message("Starting parallel computation... \n") 
  registerDoParallel(cl)
  DeMixSC.res = foreach(j = 1:ncol(bulk.mat), .export = "DeMixSC.core", .packages = c("nnls", "base")) %dopar% {
    res = DeMixSC.core(Y = bulk.mat[,j], X = ref, iter.max = iter.max, eps = eps, print.sample = print.sample, output.file = output.file)
    res
  }
  stopCluster(cl)
  message("Parallel computation completed. \n") 
  
  # Compile results.
  combined.proportion = do.call(cbind, lapply(DeMixSC.res, function(x) x$proportion))
  rownames(combined.proportion) = colnames(ref)
  colnames(combined.proportion) = colnames(bulk.mat)
  
  combined.converge = unlist(lapply(DeMixSC.res, function(x) x$converge))
  names(combined.converge) = colnames(bulk.mat)
  
  combined.weights = do.call(cbind, lapply(DeMixSC.res, function(x) x$weights))
  colnames(combined.weights) = colnames(bulk.mat)
  
  message("Deconvolution tasks completed! \n") 
  return(list(cell.type.proportions = combined.proportion,
              converge = combined.converge,
              weights = combined.weights))
  
} 















