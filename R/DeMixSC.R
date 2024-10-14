#' @title The main function of DeMixSC package
#'
#' @author Shuai Guo, Xiaoqian Liu, Wenyi Wang
#'
#' @description
#' Implements an integrated workflow for preprocessing and deconvolving bulk RNA-seq datasets using single-cell reference data.
#'
#'
#' @param option A character string specifying the mode of operation.
#' Options include:
#' - `"retina"`: Runs the deconvolution using the provided retina benchmark data and consensus reference.
#' - `"hgsc"`: Runs the deconvolution using the provided high-grade serous ovarian carcinoma (HGSC) benchmark data and consensus reference.
#' - `"user.defined"`: Allows the user to provide their own reference, benchmark, and bulk cohort data for deconvolution.
#'
#' @param benchmark.mode Indicates whether DeMixSC should run in benchmark mode.
#' When set to TRUE (the default), DeMixSC will use predefined benchmark datasets for deconvolution.
#' If set to FALSE, the user must provide their own target bulk matrix (`mat.target`) for deconvolution in real data mode.
#'
#' @param mat.target A gene expression matrix of target bulk cohort data to be deconvolved, where rows are genes and columns are samples.
#' If `NULL`, the function runs in benchmark mode (i.e., deconvolves the benchmark data `mat.a` and `mat.b`).
#'
#' @param mat.a A gene expression matrix 'a' from the benchmark dataset (e.g., bulk matrix), where rows are genes and columns are samples.
#' Used when `option = "user.defined"`.
#'
#' @param mat.b A gene expression matrix 'b' from the benchmark dataset (e.g., pseudo-bulk matrix), where rows are genes and columns are samples.
#' Used when `option = "user.defined"`.
#'
#' @param reference A cell-type-specific gene expression matrix, where rows are genes and columns are cell types.
#'
#' @param min.expression A numeric value that defines the minimum expression threshold for filtering genes.
#' Genes with expression values below this threshold across all samples will be removed. Default is 1.
#'
#' @param scale.factor A numeric value used to scale the gene expression data.
#' This helps normalize the data by making the total expression in each sample comparable. Default is 1e5.
#'
#' @param adjp.cutoff A numeric value representing the adjusted p-value cutoff (a threshold for statistical significance after accounting for multiple tests).
#' Genes with adjusted p-values greater than this threshold are not considered significantly differentially expressed. Default is 0.05.
#'
#' @param log2fc.cutoff A numeric value defining the log2 fold change cutoff for selecting differentially expressed genes.
#' Genes with a log2 fold change smaller than this threshold are not considered. Default is 0.75.
#'
#' @param top.ranked.genes An integer value indicating the number of top-ranked differentially expressed genes to be selected based on their mean expression.
#' These genes are used for technological discrepancy adjustments. Default is 5000.
#'
#' @param adj.factor A numeric value used as a scaling factor to adjust the expression levels of discrepancy genes. Default is 1000.
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
#' **Benchmark Mode:**
#' When `mat.target` is not provided (i.e., benchmark mode), the function returns:
#' \item{bulk.prop}{A matrix of estimated cell-type proportions for the bulk benchmark data (`mat.a`), where rows are cell types and columns are samples.}
#' \item{pseudobulk.prop}{A matrix of estimated cell-type proportions for the pseudo-bulk benchmark data (`mat.b`), where rows are cell types and columns are samples.}
#'
#' **Real Data Mode:**
#' When `mat.target` is provided, the function returns:
#'  \item{proportion}{A matrix of estimated cell-type proportions for the provided bulk cohort (`mat.target`), with rows are cell types and columns are samples.}
#'  \item{converge}{A vector indicating the number of iterations for the deconvolution algorithm to converge for each sample in the bulk cohort.}
#'
#'
#' @export

DeMixSC = function(option = c("retina", "hgsc", "user.defined"),
                   benchmark.mode = TRUE,
                   mat.target = NULL,
                   mat.a = NULL,
                   mat.b = NULL,
                   reference = NULL,
                   min.expression = 1,
                   scale.factor = 1e5,
                   adjp.cutoff = 0.05,
                   log2fc.cutoff = 0.75,
                   top.ranked.genes = 5000,
                   adj.factor = 1000,
                   c=2, iter.max = 2000,
                   eps = 0.001,
                   nthread = (detectCores()[1]-1),
                   print.sample = 1,
                   output.file = "./DeMixSC.save.iter.txt") {

  ###############################################
  # Check if data is provided appropriately
  ###############################################

  if (benchmark.mode && is.null(mat.target)) {

    # Benchmark mode is enabled and no target matrix (mat.target) is provided.
    # DeMixSC will proceed in benchmark mode as expected.
    cat("Data input check: \n", "Benchmark mode is enabled, and no `mat.target` is provided.\n", "DeMixSC will proceed in benchmark mode.\n")

  } else if (!benchmark.mode && !is.null(mat.target)) {

    # Benchmark mode is disabled, and mat.target is provided.
    # DeMixSC will proceed in real data mode as expected.
    cat("Data input check: \n", "Benchmark mode is disabled, and `mat.target` is provided.\n", "DeMixSC will proceed in real data mode.\n")

  } else if (benchmark.mode && !is.null(mat.target)) {

    # Benchmark mode is enabled, but mat.target is also provided.
    # This may indicate that the user forgot to disable benchmark mode.
    # DeMixSC will automatically detect the target matrix and proceed in real data mode.
    cat("Data input check: \n", "Benchmark mode is enabled, but `mat.target` is provided.\n", "DeMixSC will switch to real data mode.\n")

  } else {

    # Benchmark mode is disabled, but no mat.target is provided.
    # DeMixSC cannot proceed without the target matrix, so it stops and prompts the user.
    stop("Invalid data input.\n", "`mat.target` is missing while benchmark mode is disabled.\n", "Please provide `mat.target` or enable benchmark mode.")

  }

  ###############################################
  # Deconvolve retina tissues using the prebuilt retina benchmark dataset and consensus reference.
  #
  # The pre-built reference matrix includes expression profiles of 10 distinct cell types:
  # - Seven major cell types: Amacrine cells (ACs), bipolar cells (BCs), Cone cells,
  #   horizontal cells (HCs), Müller glia cells (MGs), retinal ganglion cells (RGCs), and Rod cells.
  # - Three minor cell types: Astrocytes, microglia cells, and retinal pigment epithelium (RPE).
  #
  # The cell types were annotated using known marker genes from Liang et al. (2019).
  #
  # The unmatched AMD retinal samples used in our study were obtained from Ratnapriya et al. (2019).
  ###############################################

  if (option == "retina") {

    cat("Running DeMixSC in retina mode with the provided retina benchmark data and reference...\n")

    # Load the consensus reference matrix specific to retina tissues.
    # This reference contains the expression profiles of 10 different cell types as described above.
    reference = get("retina_consensus_reference", envir = asNamespace("DeMixSC"))

    # Load the retina benchmark dataset that contains bulk and pseudo-bulk expression matrices.
    mat.a = get("retina_benchmark_bulk_batch2", envir = asNamespace("DeMixSC"))
    mat.b = get("retina_benchmark_pseudobulk_batch2", envir = asNamespace("DeMixSC"))

    # Preprocess the input data using the DeMixSC preprocessing function.
    # This step prepares the target bulk along with the reference matrix for the deconvolution analysis.
    input = DeMixSC.preprocessing(mat.target = mat.target,
                                  mat.a = mat.a,
                                  mat.b = mat.b,
                                  reference = reference,
                                  min.expression = min.expression,
                                  scale.factor = scale.factor,
                                  adjp.cutoff = adjp.cutoff,
                                  log2fc.cutoff = log2fc.cutoff,
                                  top.ranked.genes = top.ranked.genes,
                                  adj.factor = adj.factor)

    # Run the deconvolution algorithm using DeMixSC with parallel processing.
    # The input includes the adjusted bulk matrix (mat.target.adj) and reference matrix (reference.adj).
    res = DeMixSC.wnnls.parallel(bulk.mat = input$mat.target.adj,
                                 ref = input$reference.adj,
                                 c = c,
                                 iter.max = iter.max,
                                 eps = eps,
                                 nthread = nthread,
                                 print.sample = print.sample,
                                 output.file = output.file)

    return(res)

  }

  ###############################################
  # Deconvolve high-grade serous ovarian carcinoma (HGSC) samples using the provided HGSC benchmark dataset and consensus reference.
  #
  # The pre-built reference matrix includes expression profiles of 13 cell types:
  # - Epithelial cells, endothelial cells, fibroblasts, B cells, plasma cells,
  #   natural killer (NK) cells, innate lymphoid cells (ILCs), monocytes, macrophages, mast cells,
  #   dendritic cells (DCs), plasmacytoid dendritic cells (pDCs), and T cells.
  #
  # These cell types were annotated from the original study Hippen et al. (2023).
  #
  # The unmatched HGSC samples used in our study were obtained from Lee et al. (2020).
  ###############################################

  if (option == "hgsc") {

    cat("Running DeMixSC in HGSC mode with the provided HGSC benchmark data and reference...\n")

    # Load the consensus reference matrix specific to HGSC tissues.
    # This reference contains the expression profiles of 13 different cell types as described above.
    reference = get("hgsc_consensus_reference", envir = asNamespace("DeMixSC"))

    # Load the HGSC benchmark dataset that contains bulk and pseudo-bulk expression matrices.
    mat.a = get("hgsc_benchmark_bulk", envir = asNamespace("DeMixSC"))
    mat.b = get("hgsc_benchmark_pseudobulk", envir = asNamespace("DeMixSC"))

    # Preprocess the input data using the DeMixSC preprocessing function.
    # This step prepares the target bulk along with the reference matrix for the deconvolution analysis.
    input = DeMixSC.preprocessing(mat.target = mat.target,
                                  mat.a = mat.a,
                                  mat.b = mat.b,
                                  reference = reference,
                                  min.expression = min.expression,
                                  scale.factor = scale.factor,
                                  adjp.cutoff = adjp.cutoff,
                                  log2fc.cutoff = log2fc.cutoff,
                                  top.ranked.genes = top.ranked.genes,
                                  adj.factor = adj.factor)

    # Run the deconvolution algorithm using DeMixSC with parallel processing.
    # The input includes the adjusted bulk matrix (mat.target.adj) and reference matrix (reference.adj).
    res = DeMixSC.wnnls.parallel(bulk.mat = input$mat.target.adj,
                                 ref = input$reference.adj,
                                 c = c,
                                 iter.max = iter.max,
                                 eps = eps,
                                 nthread = nthread,
                                 print.sample = print.sample,
                                 output.file = output.file)

    return(res)
  }

  ###############################################
  # DeMixSC deconvolution with user provided benchmark dataset and reference.
  #
  # If the unmatched bulk cohort (mat.target) is not provided (i.e., mat.target = NULL),
  # the DeMixSC.preprocessing function will run in benchmark mode.
  # In benchmark mode, it processes the provided reference and benchmark dataset,
  # and performs deconvolution for both benchmark bulk and pseudo-bulk samples.
  #
  # If mat.target is provided, DeMixSC.preprocessing runs in real data mode,
  # where the provided bulk cohort is deconvolved with descrepancy adjustment.
  ###############################################

  if (option == "user.defined") {

    cat("Running DeMixSC in user-defined mode...\n")

    # Preprocess the input data using the DeMixSC preprocessing function.
    # This step prepares the target bulk along with the reference matrix for the deconvolution analysis.

    input = DeMixSC.preprocessing(mat.target = mat.target,
                                  mat.a = mat.a,
                                  mat.b = mat.b,
                                  reference = reference,
                                  min.expression = min.expression,
                                  scale.factor = scale.factor,
                                  adjp.cutoff = adjp.cutoff,
                                  log2fc.cutoff = log2fc.cutoff,
                                  top.ranked.genes = top.ranked.genes,
                                  adj.factor = adj.factor)


    if (is.null(mat.target)) {

      cat("Start to deconvolve bulk benchmark data...\n")

      # Run the deconvolution algorithm using DeMixSC with parallel processing on benchmark bulk samples.
      res.bulk = DeMixSC.wnnls.parallel(bulk.mat = input$mat.a.adj,
                                        ref = input$reference.adj,
                                        c = c,
                                        iter.max = iter.max,
                                        eps = eps,
                                        nthread = nthread,
                                        print.sample = NULL,
                                        output.file = NULL)

      cat("Start to deconvolve pseudo-bulk benchmark data...\n")

      # Run the deconvolution algorithm using DeMixSC with parallel processing on benchmark pseudo-bulk samples.
      res.pseudobulk = DeMixSC.wnnls.parallel(bulk.mat = input$mat.b,
                                              ref = input$reference,
                                              c = c,
                                              iter.max = iter.max,
                                              eps = eps,
                                              nthread = nthread,
                                              print.sample = NULL,
                                              output.file = NULL)
      res = list(bulk.prop = res.bulk$cell.type.proportions,
                 pseudobulk.prop = res.pseudobulk$cell.type.proportions)

      return(res)

    } else {

      # Run the deconvolution algorithm using DeMixSC with parallel processing.
      # The input includes the adjusted bulk matrix (mat.target.adj) and reference matrix (reference.adj).
      res = DeMixSC.wnnls.parallel(bulk.mat = input$mat.target.adj,
                                   ref = input$reference.adj,
                                   c = c,
                                   iter.max = iter.max,
                                   eps = eps,
                                   nthread = nthread,
                                   print.sample = print.sample,
                                   output.file = output.file)

      return(res)
    }
  }
}
