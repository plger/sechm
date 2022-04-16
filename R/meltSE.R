#' meltSE
#'
#' Melts a SE object into a \code{\link[ggplot2]{ggplot}}-ready long data.frame.
#'
#' @param x An object of class
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @param features A vector of features (i.e. row.names) to include. Use 
#'   `features=NULL` to include all.
#' @param assayName The name(s) of the assay(s) to use. If NULL and the assays 
#'   are named, all of them will be included.
#' @param colDat.columns The colData columns to include (defaults includes all).
#'   Use `colDat.columns=NA` in order not to include any.
#' @param rowDat.columns The rowData columns to include (default all). Use 
#'   `rowData=NA` to not include any.
#'
#' @return A data.frame.
#'
#' @examples
#' data("Chen2017", package="sechm")
#' head(meltSE(Chen2017,"Fos"))
#'
#' @import SummarizedExperiment
#' @export
meltSE <- function(x, features, assayName=NULL, colDat.columns=NULL,
                   rowDat.columns=NULL){
  features <- intersect(features, row.names(x))
  if(is.null(colDat.columns)) colDat.columns <- colnames(colData(x))
  if(all(is.na(colDat.columns))) colDat.columns <- c()
  colDat.columns <- intersect(colDat.columns, colnames(colData(x)))
  if(is.null(rowDat.columns)) rowDat.columns <- colnames(rowData(x))
  if(all(is.na(rowDat.columns))) rowDat.columns <- c()
  rowDat.columns <- intersect(rowDat.columns, colnames(rowData(x)))
  if(is.null(assayName) && !is.null(assayNames(x)) ) assayName <- assayNames(x)
  if(is.null(assayName)){
    a <- list(value=assay(x))
  }else{
    a <- assays(x)[assayName]
  }
  if(is.numeric(assayName)) names(a) <- paste0("assay", assayName)
  a <- lapply(a, FUN=function(x) x[features,,drop=FALSE])
  df <- data.frame( feature=rep(features,ncol(x)),
                    sample=rep(colnames(x), each=length(features)) )
  for(f in colDat.columns) df[[f]] <- rep(colData(x)[[f]],each=length(features))
  for(f in rowDat.columns) df[[f]] <- rep(rowData(x)[features,f], ncol(x))
  for(v in names(a)) df[[v]] <- as.numeric(a[[v]])
  df
}
