#' @title Functions used in the pre-processing steps implemented in DeMixSC

#' @author Shuai Guo, Xiaoqian Liu, Wenyi Wang

############################################################

#' @title Pre-process the single-cell Seurat Object
#'
#' @description Take the annotated Seurat object as input to build a reference matrix.
#'
#' @param Seurat.obj The pre-annotated single-cell Seurat object.
#' @param annotation A scalar or a character string to designate the column for cell-type information in the Seurat object.
#'
#' @return
#' \item{reference}{A matrix of the reference expression,
#'                  with genes on the rows and cell types on the columns.}
#' \item{theta}{A matrix of the expression abundance relative to total UMI count per cell type,
#'              with genes on the rows and cell types on the columns.}
#' \item{cell.size}{A vector of the mean cell size for each cell type.}
#'
#' @export

DeMixSC.ref = function(Seurat.obj, annotation){

  # Obtain the data from the Seurat object
  mat.norm  = Seurat.obj@assays$RNA@data
  cell.type = as.matrix(Seurat.obj@meta.data[,annotation]) ## maintain column and row names in the original matrix.

  # Filter genes with zero expression in all columns
  nz.gene = rownames(mat.norm)[(Matrix::rowSums(mat.norm) != 0 )]
  mat.norm = mat.norm[nz.gene, , drop = FALSE]

  # Calculate theta for each cell-type
  theta.ct = sapply(unique(cell.type), function(ct){
    x = mat.norm[,cell.type %in% ct, drop = FALSE]
    Matrix::rowSums(x)/sum(x)
  })

  # Calculate cell size for each cell-type
  size.ct = sapply(unique(cell.type), function(ct){
    x = mat.norm[, cell.type %in% ct, drop = FALSE]
    sum(x)/ncol(x)
  })

  # Construct the cell-type-specific reference matrix for each cell type
  ref.ct = t(t(theta.ct)*size.ct)

  return(list(reference = ref.ct,
              cell.size = size.ct,
              theta = theta.ct))
}

#' @title Identify genes less affected by technological discrepancies
#'
#' @description
#' Using the benchmark dataset as input, this function discerns genes that remain stable across technological platforms from those that exhibit high discrepancies.
#'
#' @param mat.a Expression matrix 'a' from the benchmark dataset, with genes as rows and cell types as columns (potentially the bulk matrix).
#' @param mat.b Expression matrix 'b' from the benchmark dataset, with genes as rows and cell types as columns (potentially the pseudo-bulk matrix).
#' @param cutoff A cutoff P-value threshold to classify genes influenced by technological discrepancies. Default is 0.05.
#'
#' @return
#'  \item{non.discrepancies.genes}{A list of genes minimally influenced by technological discrepancies.}
#'  \item{mat.a.normalized}{Normalized expression matrix 'a' from the benchmark dataset.}
#'  \item{mat.b.normalized}{Normalized expression matrix 'b' from the benchmark dataset.}
#'  \item{t.test.res}{Paired t-test results table for genes detailing various parameters, with genes on the rows and parameters on the columns.}
#'
#' @export

DeMixSC.discrepancy = function(mat.a, mat.b, cutoff = 0.05){

  # Determine sample sizes
  n.mat.a = ncol(mat.a)
  n.mat.b = ncol(mat.b)
  n.total = n.mat.a + n.mat.b

  # Normalize the expression matrix
  commom.genes = intersect(row.names(mat.a), row.names(mat.b))
  df_norm = preprocessCore::normalize.quantiles(as.matrix(cbind(mat.a[commom.genes,],mat.b[commom.genes,])))
  colnames(df_norm) = c(colnames(mat.a), colnames(mat.b))
  rownames(df_norm) = commom.genes

  # Remove rows with zero expression across all samples
  if(sum(df_norm == 0) > 0){
    df_norm = df_norm[df_norm !=0,]
  }

  # Perform differential expression analysis using paired t-test
  res = matrix(nrow = nrow(df_norm), ncol = 8)
  row.names(res) = row.names(df_norm)
  colnames(res) = c("pval", "baseMean.mat.a", "var.mat.a", "baseMean.mat.b", "ver.mat.b", "baseMean", "abs.dif", "log2FoldChange")
  for (i in 1:nrow(df_norm)) {
      p = t.test(df_norm[i, 1:n.mat.a], df_norm[i, (n.mat.a+1):n.total], paired = T, exact=FALSE)$p.value
      baseMean.mat.a = mean(df_norm[i,1:n.mat.a])
      var.mat.a = var(df_norm[i,1:n.mat.a])
      baseMean.mat.b = mean(df_norm[i,(n.mat.a+1):n.total])
      var.mat.b = var(df_norm[i,(n.mat.a+1):n.total])
      baseMean = mean(df_norm[i,])
      dif = abs(baseMean.mat.a - baseMean.mat.b)
      log2fc = log2(baseMean.mat.a/baseMean.mat.b)
      res[i,] = c(p, baseMean.mat.a, var.mat.a, baseMean.mat.b, var.mat.b, baseMean, dif, log2fc)
  }
  padj= p.adjust(res[,"pval"], method = "fdr") # Adjust p-values with FDR.
  res = cbind(res, padj)

  # Outputs
  return(list(non.discrepancy.genes = names(which(res[,"padj"] >= cutoff)),
              mat.a.normalized = df_norm[,1:n.mat.a],
              mat.b.normalized = df_norm[,(n.mat.a+1):n.total],
              t.test.res = res))
}

#' @title Preprocess a benchmark dataset for DeMixSC deconvolution.
#'
#' @description
#' Prepares and adjusts the benchmark dataset to handle technological discrepancies for the DeMixSC deconvolution.
#'
#' @param Seurat.obj The pre-annotated single-cell Seurat object.
#' @param annotation A column name or scalar indicating cell-type information within the Seurat object.
#' @param mat.a  mRNA expression matrix 'a' from the benchmark dataset (possibly bulk matrix), with genes on the rows and cell types on the columns.
#' @param mat.b  mRNA expression matrix 'b' from the benchmark dataset (possibly pseudo-bulk matrix), with genes on the rows and cell types on the columns.
#' @param cutoff Threshold to differentiate genes based on inter-platform discrepancies. Default is 0.05.
#' @param scale.factor Scale factor applied to the relative expression abundances matrix. Default is 1e6.
#'
#' @return
#'  \item{mat.a.adj}{Adjusted relative expression matrix 'a', with genes on the rows and samples on the columns.}
#'  \item{mat.b}{Relative expression matrix 'b', with genes on the rows and samples on the columns.}
#'  \item{reference.adj}{Adjusted reference matrix used for deconvolving mat.a.adj}
#'  \item{reference}{Original reference matrix used for deconvolving mat.b}
#'  \item{theta.adj}{Adjusted relative abundance matrix of reference expression.}
#'  \item{cell.size}{A vector of the mean cell size for each cell type.}
#'  \item{discrepancy.genes}{A list of genes highly affected with technological discrepancies.}
#'  \item{non.discrepancy.genes}{A list of genes hardly affected by technological discrepancies.}
#'
#' @export

DeMixSC.prep.benchmark = function(Seurat.obj, annotation, mat.a, mat.b, cutoff = 0.05, scale.factor = 1e6) {

  # Construct the reference expression matrix from Seurat object
  ref = DeMixSC.ref(Seurat.obj, annotation)

  # Determine genes with discrepancies between matrices
  discrepancy = DeMixSC.discrepancy(mat.a, mat.b, cutoff)
  commom.genes = intersect(row.names(ref$reference), row.names(discrepancy$mat.a.normalized))
  discrepancy.genes = setdiff(commom.genes, intersect(commom.genes, discrepancy$non.discrepancy.genes))

  # Extract relevant attributes
  cell.size = ref$cell.size
  theta     = ref$theta[commom.genes,]
  reference = ref$reference[commom.genes,]

  # Normalize the expression matrices for relative expression abundances
  mat.a.final = apply(discrepancy$mat.a.normalized[commom.genes,], 2, function(x){ x/sum(x)*(scale.factor) })
  mat.b.final = apply(discrepancy$mat.b.normalized[commom.genes,], 2, function(x){ x/sum(x)*(scale.factor) })

  # Adjust expression values of genes influenced by technological discrepancies
  for (i in discrepancy.genes) {
    log2_basemean = log2(discrepancy$t.test.res[i,"baseMean"])
    mat.a.final[i,] = mat.a.final[i,]/log2_basemean
    reference[i,] = reference[i,]/log2_basemean
    theta[i,] = theta[i,]/log2_basemean
  }

  # Return the processed and adjusted matrices, references, and additional attributes
  return(list(mat.a.adj = mat.a.final,
              mat.b = mat.b.final,
              reference.adj = reference,
              reference = ref$reference[commom.genes,],
              theta.adj = theta,
              cell.size = cell.size,
              discrepancy.genes = discrepancy.genes,
              non.discrepancy.genes = discrepancy$non.discrepancy.genes))
}

#' @title Preprocess a benchmark dataset for deconvolution of unmatched bulk cohort using DeMixSC
#'
#' @description
#' Adjusts and prepares the target bulk matrix and a reference matrix for deconvolution.
#'
#' @param Seurat.obj The pre-annotated single-cell Seurat object.
#' @param annotation A column name or scalar indicating cell-type information within the Seurat object.
#' @param mat.a mRNA expression matrix 'a' from the benchmark dataset (possibly bulk matrix), with genes on the rows and cell types on the columns.
#' @param mat.b mRNA expression matrix 'b' from the benchmark dataset (possibly pseudo-bulk matrix), with genes on the rows and cell types on the columns.
#' @param mat.target Targeted bulk cohort matrix for deconvolution, with genes on the rows and cell types on the columns.
#' @param cutoff Threshold to differentiate genes based on inter-platform discrepancies. Default is 0.05.
#' @param scale.factor Scale factor applied to the relative expression abundances matrix. Default is 1e6.
#' @param combat Logical indicator for using ComBat to align the unmatched bulk cohort with the benchmark dataset's bulk data. Default is TRUE.
#'
#' @return
#'  \item{mat.target.adj}{Adjusted target bulk cohort's relative expression matrix.}
#'  \item{reference.adj}{Adjusted reference matrix suitable for deconvolution}
#'  \item{theta.adj}{Adjusted relative abundance matrix of reference expression.}
#'  \item{cell.size}{A vector of the mean cell size for each cell type.}
#'  \item{discrepancy.genes}{A list of genes highly affected with technological discrepancies.}
#'  \item{non.discrepancy.genes}{A list of genes hardly affected by technological discrepancies.}
#'
#' @export

DeMixSC.prep = function(Seurat.obj, annotation,
                        mat.a, mat.b, mat.target,
                        cutoff = 0.05, scale.factor = 1e6, combat = TRUE) {

  # Remove zero genes across matrices
  if (sum(rowSums(mat.target) == 0) != 0) {
    mat.target = mat.target[rowSums(mat.target == 0) != ncol(mat.target),]
  }
  if (sum(rowSums(mat.a) == 0) != 0) {
    mat.a = mat.a[rowSums(mat.a == 0) != ncol(mat.a),]
  }
  if (sum(rowSums(mat.b) == 0) != 0) {
    mat.b = mat.b[rowSums(mat.b == 0) != ncol(mat.b),]
  }

  # Match genes across matrices
  bulk.genes = Reduce(intersect, list(row.names(mat.a), row.names(mat.b), row.names(mat.target)))
  mat.target = mat.target[bulk.genes,]
  mat.a = mat.a[bulk.genes,]
  mat.b = mat.b[bulk.genes,]

  # Normalize matrices and align the bulk RNA-seq data
  ncol.mat.target = ncol(mat.target) # number of samples in the targeted bulk to be deconvolved
  ncol.mat.a = ncol(mat.a) # number of samples in the bulk data from benchmark dataset
  bulk.input = preprocessCore::normalize.quantiles(as.matrix(cbind(mat.target, mat.a)))
  if (combat) {
    mat.align = sva::ComBat(dat = bulk.input, batch = c(rep("target", ncol.mat.target), rep("a", ncol.mat.a)))
  } else {
    mat.align = bulk.input
  }

  mat.align.norm = preprocessCore::normalize.quantiles(as.matrix(cbind(mat.align, mat.b)))
  rownames(mat.align.norm) = bulk.genes
  colnames(mat.align.norm) = colnames(cbind(mat.target, mat.a, mat.b))

  mat.target.align = mat.align.norm[,1:ncol.mat.target]
  mat.a.align = mat.align.norm[,(ncol.mat.target+1): (ncol.mat.target + ncol.mat.a) ]
  mat.b.align = mat.align.norm[,(ncol.mat.target + ncol.mat.a + 1):ncol(mat.align.norm)]

  # Identify genes with discrepancies post-alignment
  discrepancy = DeMixSC.discrepancy(mat.a.align, mat.b.align, cutoff)

  # Prepare the final input matrix
  ref = DeMixSC.ref(Seurat.obj, annotation)
  mat.target.align.final = apply(mat.target.align, 2, function(x){ x/sum(x)*(scale.factor) })

  commom.genes = intersect(row.names(ref$reference), row.names(discrepancy$mat.a.normalized))
  cell.size = ref$cell.size
  theta     = ref$theta[commom.genes,]
  reference = ref$reference[commom.genes,]
  mat.target.align.final = mat.target.align.final[commom.genes,]

  # Adjust for technological discrepancies in genes
  discrepancy.genes = setdiff(commom.genes, intersect(commom.genes, discrepancy$non.discrepancy.genes))
  for (i in discrepancy.genes) {
    log2_basemean = log2(discrepancy$t.test.res[i,"baseMean"])
    mat.target.align.final[i,] = mat.target.align.final[i,]/log2_basemean
    reference[i,] = reference[i,]/log2_basemean
    theta[i,] = theta[i,]/log2_basemean
  }

  return(list(mat.target.adj = mat.target.align.final,
              reference.adj = reference,
              theta.adj = theta,
              cell.size = cell.size,
              discrepancy.genes = discrepancy.genes,
              non.discrepancy.genes = discrepancy$non.discrepancy.genes))
}
