#' @title Preprocess the input data for DeMixSC wNNLS framework
#'
#' @author Shuai Guo, Xiaoqian Liu, Wenyi Wang
#'
#' @description
#' This function preprocesses input data for DeMixSC by handling two modes:
#' 1) real data mode when target bulk data is provided, and
#' 2) benchmark mode when target bulk data is not provided.
#' The function performs filtering, normalization, and adjusts genes with technological discrepancies.
#'
#'
#' @param mat.target A gene expression matrix of target bulk cohort data to be deconvolved, where rows are genes and columns are samples.
#' If `NULL`, the function runs in benchmark mode (i.e., deconvolves the benchmark data `mat.a` and `mat.b`).
#'
#' @param mat.a A gene expression matrix 'a' from the benchmark dataset (e.g., bulk matrix), where rows are genes and columns are samples.
#'
#' @param mat.b A gene expression matrix 'b' from the benchmark dataset (e.g., pseudo-bulk matrix), where rows are genes and columns are samples.
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
#'
#' @return
#' **Benchmark Mode:**
#' When `mat.target` is not provided (i.e., benchmark mode), the function returns:
#'  \item{mat.a.adj}{Adjusted relative expression matrix 'a', with genes on the rows and samples on the columns.}
#'  \item{mat.b}{Relative expression matrix 'b', with genes on the rows and samples on the columns.}
#'  \item{reference.adj}{Adjusted reference matrix used for deconvolving mat.a.adj.}
#'  \item{reference}{Original reference matrix used for deconvolving mat.b.}
#'  \item{discrepancy.genes}{A list of genes highly affected with technological discrepancies.}
#'  \item{non.discrepancy.genes}{A list of genes without detected technological discrepancies.}
#'  \item{de.res}{Results of differential expression analysis between bulk and pseudobulk data, including p-values, fold changes, and mean expression levels.}
#'
#' **Real Data Mode:**
#' When `mat.target` is provided, the function returns:
#'  \item{mat.target.adj}{Adjusted expression matrix for the target bulk cohort after technological discrepancy correction.}
#'  \item{reference.adj}{Adjusted reference matrix after correcting for technological discrepancies.}
#'  \item{discrepancy.genes}{A list of genes highly affected with technological discrepancies.}
#'  \item{non.discrepancy.genes}{A list of genes without detected technological discrepancies.}
#'  \item{de.res}{Results of differential expression analysis between bulk and pseudobulk data, including p-values, fold changes, and mean expression levels.}
#'
#'
#' @export

DeMixSC.preprocessing = function(mat.target = NULL,
                                 mat.a = NULL,
                                 mat.b = NULL,
                                 reference = NULL,
                                 min.expression = 1,
                                 scale.factor = 1e5,
                                 adjp.cutoff = 0.05,
                                 log2fc.cutoff = 0.75,
                                 top.ranked.genes = 5000,
                                 adj.factor = 1000) {

  cat("DeMixSC preprocessing has started. \n")

  # Check if the target bulk cohort (mat.target) is provided.
  if (is.null(mat.target)) {

    ########################################################
    # If mat.target is NULL, the function will adjust inputs
    # for running in benchmark mode.
    ########################################################

    # Output message for benchmark mode.
    cat("No target bulk provided. Running the DeMixSC benchmark mode. \n")

    # Remove genes from mat.a and mat.b where all values are less than min.expression.
    if (any(rowSums(mat.a < min.expression) == ncol(mat.a))) {
      mat.a = mat.a[rowSums(mat.a < min.expression) != ncol(mat.a), ]
    }
    if (any(rowSums(mat.b < min.expression) == ncol(mat.b))) {
      mat.b = mat.b[rowSums(mat.b < min.expression) != ncol(mat.b), ]
    }

    # Determine sample sizes for mat.a and mat.b.
    n.mat.a = ncol(mat.a)
    n.mat.b = ncol(mat.b)
    n.total = n.mat.a + n.mat.b

    # Find common genes present in reference, mat.a, and mat.b.
    common.genes = Reduce(intersect, list(row.names(reference), row.names(mat.a), row.names(mat.b)))

    # Subset reference, mat.a, and mat.b to only include common genes.
    reference = reference[common.genes,]
    mat.a = mat.a[common.genes,]
    mat.b = mat.b[common.genes,]

    # Perform quantile normalization on the combined matrix of mat.a and mat.b.
    df_norm = preprocessCore::normalize.quantiles(as.matrix(cbind(mat.a, mat.b)))
    rownames(df_norm) = common.genes; colnames(df_norm) = c(colnames(mat.a), colnames(mat.b))
    df_norm = apply(df_norm, 2, function(x){ x/sum(x)*(scale.factor) })

    # Perform paired t-test to identify genes with significant differences between mat.a and mat.b.
    pval = suppressWarnings(apply(df_norm, 1, function(row) {
      t.test(row[1:n.mat.a], row[(n.mat.a+1):n.total], paired = T, exact = F)$p.value
    }))
    padj = p.adjust(pval, method = "fdr")
    baseMean.mat.a = abs(rowMeans(df_norm[, 1:n.mat.a]))
    baseMean.mat.b = abs(rowMeans(df_norm[, (n.mat.a+1):n.total]))
    baseMean = (baseMean.mat.a + baseMean.mat.b) / 2
    abs.dif = abs(baseMean.mat.a - baseMean.mat.b)
    log2fc = log2(baseMean.mat.a / baseMean.mat.b)

    # Compile differential expression results into a table.
    res = cbind(pval, padj, baseMean.mat.a, baseMean.mat.b, baseMean, abs.dif, log2fc)
    row.names(res) = row.names(df_norm)
    colnames(res) = c("pval", "padj", "baseMean.mat.a", "baseMean.mat.b", "baseMean", "abs.dif", "log2.fold.change")

    # Select differentially expressed (DE) genes based on log2 fold change and adjusted p-value.
    res.de = res[which(abs(res[, "log2.fold.change"]) >= log2fc.cutoff & res[, "padj"] <= adjp.cutoff),]

    # Rank DE genes by base mean expression level and identify discrepancy genes (top ranked DE genes).
    res.de.sorted = res.de[order(res.de[, "baseMean"], decreasing = TRUE), ]
    discrepancy.genes = rownames(res.de.sorted)[1:1:min(top.ranked.genes, nrow(res.de.sorted))]
    non.discrepancy.genes = setdiff(common.genes, discrepancy.genes)

    # Adjust mat.a and reference for technological discrepancies (scaling discrepancy genes).
    mat.a.final = apply(df_norm[, 1:n.mat.a], 2, function(x){ x/sum(x)*(scale.factor) })
    mat.b.final = apply(df_norm[, (n.mat.a+1):n.total], 2, function(x){ x/sum(x)*(scale.factor) })
    reference.adj = reference
    for (i in discrepancy.genes) {
      mat.a.final[i,] = mat.a.final[i,]/adj.factor
      reference.adj[i,] = reference.adj[i,]/adj.factor
    }

    # Return the processed matrices, references, and additional attributes.
    return(list(mat.a.adj = mat.a.final,
                mat.b = mat.b.final,
                reference.adj = reference.adj,
                reference = reference,
                discrepancy.genes = discrepancy.genes,
                non.discrepancy.genes = non.discrepancy.genes,
                de.res = res.de.sorted))

  } else {

    ###############################################
    # If mat.target is provided, function runs in real data mode.
    ###############################################

    # Output message for benchmark mode.
    cat("Target bulk provided. Running the DeMixSC real data mode. \n")

    # Remove genes from target bulk, mat.a and mat.b where all values are less than min.expression.
    if (any(rowSums(mat.a < min.expression) == ncol(mat.a))) {
      mat.a = mat.a[rowSums(mat.a < min.expression) != ncol(mat.a), ]
    }
    if (any(rowSums(mat.b < min.expression) == ncol(mat.b))) {
      mat.b = mat.b[rowSums(mat.b < min.expression) != ncol(mat.b), ]
    }
    if (any(rowSums(mat.target < min.expression) == ncol(mat.target))) {
      mat.target = mat.target[rowSums(mat.target < min.expression) != ncol(mat.target), ]
    }

    # Find common genes present in reference, target bulk, mat.a, and mat.b.
    common.genes = Reduce(intersect, list(row.names(reference), row.names(mat.target), row.names(mat.a), row.names(mat.b)))
    reference = reference[common.genes,]
    mat.target = mat.target[common.genes,]
    mat.a = mat.a[common.genes,]
    mat.b = mat.b[common.genes,]

    # Determine sample sizes for target bulk, mat.a, and mat.b.
    n.mat.target = ncol(mat.target)
    n.mat.a = ncol(mat.a)
    n.mat.b = ncol(mat.b)
    n.total = n.mat.a + n.mat.b

    # Align the target bulk with the mat.a using ComBat.
    bulk.input = preprocessCore::normalize.quantiles(as.matrix(cbind(mat.target, mat.a)), keep.names = T)
    mat.align = suppressMessages(sva::ComBat(dat = bulk.input, batch = c(rep("target", n.mat.target), rep("a", n.mat.a))))

    # Perform quantile normalization on the combined matrix of target bulk, mat.a, and mat.b.
    mat.align.norm = preprocessCore::normalize.quantiles(as.matrix(cbind(mat.align, mat.b)))
    rownames(mat.align.norm) = common.genes
    colnames(mat.align.norm) = colnames(cbind(mat.target, mat.a, mat.b))
    mat.target.align = mat.align.norm[,1:n.mat.target]
    mat.a.align = mat.align.norm[,(n.mat.target+1): (n.mat.target + n.mat.a) ]
    mat.b       = mat.align.norm[,(n.mat.target + n.mat.a + 1):ncol(mat.align.norm)]

    # Perform paired t-test to identify genes with significant differences between mat.a and mat.b.
    df_norm = as.matrix(cbind(mat.a.align, mat.b))
    df_norm = apply(df_norm, 2, function(x){ x/sum(x)*(scale.factor) })
    pval = suppressWarnings(apply(df_norm, 1, function(row) {
      t.test(row[1:n.mat.a], row[(n.mat.a+1):n.total], paired = T, exact = F)$p.value
    }))
    padj = p.adjust(pval, method = "fdr")
    baseMean.mat.a = abs(rowMeans(df_norm[, 1:n.mat.a]))
    baseMean.mat.b = abs(rowMeans(df_norm[, (n.mat.a+1):n.total]))
    baseMean = (baseMean.mat.a + baseMean.mat.b) / 2
    abs.dif = abs(baseMean.mat.a - baseMean.mat.b)
    log2fc = log2(baseMean.mat.a / baseMean.mat.b)

    # Compile differential expression results into a table.
    res = cbind(pval, padj, baseMean.mat.a, baseMean.mat.b, baseMean, abs.dif, log2fc)
    row.names(res) = row.names(df_norm)
    colnames(res) = c("pval", "padj", "baseMean.mat.a", "baseMean.mat.b", "baseMean", "abs.dif", "log2.fold.change")

    # Select differentially expressed (DE) genes based on log2 fold change and adjusted p-value.
    res.de = res[which(abs(res[, "log2.fold.change"]) >= log2fc.cutoff & res[, "padj"] <= adjp.cutoff),]

    # Rank DE genes by base mean expression level and identify discrepancy genes (top ranked DE genes).
    res.de.sorted = res.de[order(res.de[, "baseMean"], decreasing = TRUE), ]
    discrepancy.genes = rownames(res.de.sorted)[1:min(top.ranked.genes, nrow(res.de.sorted))]
    non.discrepancy.genes = setdiff(common.genes, discrepancy.genes)

    # Adjust mat.a and reference for technological discrepancies (scaling discrepancy genes).
    mat.target.align.final = apply(mat.target.align, 2, function(x){ x/sum(x)*(scale.factor) })
    for (i in discrepancy.genes) {
      mat.target.align.final[i,] = mat.target.align.final[i,]/adj.factor
      reference[i,] = reference[i,]/adj.factor
    }

    # Return the processed matrices, references, and additional attributes.
    return(list(mat.target.adj = mat.target.align.final,
                reference.adj = reference,
                discrepancy.genes = discrepancy.genes,
                non.discrepancy.genes = non.discrepancy.genes,
                de.res = res.de.sorted))
  }
}
