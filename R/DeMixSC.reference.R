#' @title Build the reference matrix from given Seurat object
#'
#' @author Shuai Guo, Xiaoqian Liu, Wenyi Wang
#'
#' @description Take the annotated Seurat object as input to build a reference matrix.
#'
#'
#' @param Seurat.obj The pre-annotated single-cell Seurat object.
#'
#' @param annotation A scalar or a character string to designate the column for cell-type information in the Seurat object.
#'
#'
#' @return
#' \item{reference}{A matrix of the reference expression,
#'                  with genes on the rows and cell types on the columns.}
#' \item{theta}{A matrix of the expression abundance relative to total UMI count per cell type,
#'              with genes on the rows and cell types on the columns.}
#' \item{cell.size}{A vector of the mean cell size for each cell type.}
#'
#'
#' @export

DeMixSC.ref = function(Seurat.obj, annotation){

  # Obtain the data from the Seurat object.
  mat.norm  = Seurat.obj@assays$RNA@counts # Raw UMI counts from the single-cell data
  cell.type = as.matrix(Seurat.obj@meta.data[,annotation]) # maintain column and row names in the original matrix

  # Filter genes with zero expression in all columns.
  nz.gene = rownames(mat.norm)[(Matrix::rowSums(mat.norm) != 0 )]
  mat.norm = mat.norm[nz.gene, , drop = FALSE]

  # Calculate theta for each cell-type.
  theta.ct = sapply(unique(cell.type), function(ct){
    x = mat.norm[,cell.type %in% ct, drop = FALSE]
    Matrix::rowSums(x)/sum(x)
  })

  # Calculate cell size for each cell-type.
  size.ct = sapply(unique(cell.type), function(ct){
    x = mat.norm[, cell.type %in% ct, drop = FALSE]
    sum(x)/ncol(x)
  })

  # Construct the cell-type-specific reference matrix for each cell type.
  ref.ct = t(t(theta.ct)*size.ct)

  return(list(reference = ref.ct,
              cell.size = size.ct,
              theta = theta.ct))
}
