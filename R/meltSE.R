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
#' @param flatten Logical, whether to flatten nested data.frames.
#' @param baseDF Logical, whether to return a base data.frame (removing columns
#'   containing other objects such as atomic lists). Filtering is applied after
#'   flattening.
#'
#' @return A data.frame (or a DataFrame).
#'
#' @examples
#' data("Chen2017", package="sechm")
#' head(meltSE(Chen2017,"Fos"))
#'
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
meltSE <- function(x, features, assayName=NULL, colDat.columns=NULL,
                   rowDat.columns=NULL, flatten=TRUE, baseDF=TRUE){
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
  df <- cbind(df, .flattenDF(colData(x)[rep(seq_len(ncol(x)),each=length(features)),
                                        colDat.columns, drop=FALSE],
                             remove=!flatten))
  df <- cbind(df, .flattenDF(rowData(x)[rep(features, ncol(x)),
                                        rowDat.columns, drop=FALSE],
                             remove=!flatten))
  if(baseDF){
    for(f in colnames(df)){
      if(!is.vector(df[[f]]) && !is.factor(df[[f]])) df[[f]] <- NULL
    }
    df <- as.data.frame(df)
  }
  for(v in names(a)) df[[v]] <- as.numeric(a[[v]])
  df
}


.flattenDF <- function(x, columns=colnames(x), remove=FALSE){
  x <- DataFrame(x)[,columns,drop=FALSE]
  isdf <- unlist(lapply(x, FUN=function(x) is.data.frame(x) || is(x, "DFrame")))
  for(f in names(x)[which(isdf)]){
    y <- x[[f]]
    colnames(y) <- paste(f, colnames(y), sep=".")
    x[[f]] <- NULL
    if(!remove) x <- cbind(x, y)
  }
  if(!any(sapply(x, FUN=function(x) !is.vector(x) && !is.factor(x))))
    x <- as.data.frame(x)
  x
}