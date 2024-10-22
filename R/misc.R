#' sortRows
#'
#' @param x A numeric matrix or data.frame.
#' @param z Whether to scale rows for the purpose of calculating order.
#' @param toporder Optional verctor of categories (length=nrow(x)) on which to
#' supra-order  when sorting rows.
#' @param na.rm Whether to remove missing values and invariant rows.
#' @param method Seriation method; 'MDS_angle' (default) or 'R2E' recommended.
#' @param toporder.meth Whether to perform higher-order sorting 'before'
#' (default) or 'after' the lower-order sorting.
#'
#' @return A reordered matrix or data.frame.
#'
#' @examples
#' # random data
#' m <- matrix( round(rnorm(100,mean=10, sd=2)), nrow=10,
#'              dimnames=list(LETTERS[1:10], letters[11:20]) )
#' m
#' sortRows(m)
#'
#' @importFrom seriation seriate get_order
#' @importFrom stats sd aggregate median dist
#' @export
sortRows <- function(x, z=FALSE, toporder=NULL, na.rm=FALSE, method="MDS_angle",
                     toporder.meth="before"){
  toporder.meth <- match.arg(toporder.meth, c("before","after"))
  if(is.numeric(toporder)) toporder <- as.character(toporder)
  if(na.rm){
    w <- which( apply(x, 1, FUN = function(y){ !any(is.na(y)) }) |
                  !(apply(x, 1, na.rm=TRUE, FUN=sd) > 0) )
    x <- x[w,]
    if(!is.null(toporder)) toporder <- toporder[w]
  }
  if(is.factor(toporder)) toporder <- droplevels(toporder)
  y <- x
  if(z) y <- safescale(x, byRow=TRUE)
  if(!is.null(toporder)){
    if(toporder.meth=="before"){
      ag <- aggregate(y, by=list(toporder), na.rm=TRUE, FUN=median)
      row.names(ag) <- ag[,1]
      ag <- ag[,-1]
      if(nrow(ag)>2){
        try( ag <- sortRows(ag, z=FALSE, na.rm=FALSE, method=method),
             silent=TRUE)
      }
      ll <- split(as.data.frame(y), toporder)
      ll <- lapply(ll, FUN=function(x){
          tryCatch(sortRows(x,method=method), error=function(e) return(x))
      })
      y <- unlist(lapply(ll[row.names(ag)], FUN=row.names))
      return(x[y,])
    }else{
      o1 <- get_order(seriate(dist(y), method=method))
      oa <- aggregate(o1,by=list(top=as.factor(toporder)),FUN=median)
      oa <- oa[order(oa[,2]),1]
      toporder <- factor(as.character(toporder), levels=as.character(oa))
      return(x[order(as.numeric(toporder),o1),])
    }
  }
  ss <- seriate(dist(y), method=method)
  x[get_order(ss),]
}


.chooseAssay <- function(se, assayName=NULL, returnName=FALSE){
  if(length(assays(se))==0) stop("The object has no assay!")
  if(is.null(assayName) &&
     !is.null(dan <- metadata(se)$default_view$assay) &&
     length(dan <- intersect(dan,assayNames(se)))>0)
    assayName <- dan[1]
  assayGiven <- !is.null(assayName)
  a <- .getDef("assayName")
  if(!assayGiven && !is.null(assayNames(se))){
    assayName <- intersect(a,assayNames(se))
    if(is.null(assayName)) assayName <- assayNames(se)[1]
  }else if(!assayGiven && is.null(assayNames(se))){
    assayName <- 1L
    if(length(assays(se))>1)
      message("Assay unspecified, and multiple assays present; ",
              "will use the first one.")
  }
  if(assayGiven && !is.numeric(assayName) &&
     !any(assayName %in% assayNames(se)))
      stop("Assay '", assayName, "' not found!")

  if(length(assayName)==0) assayName <- assayNames(se)
  if(length(assayName)>1) message("Using assay '", assayName[1],"'")
  assayName <- assayName[1]
  if(returnName) return(assayName)
  assays(se)[[assayName]]
}

#' @importFrom grDevices colorRampPalette
.getHMcols <- function(cols=NULL, n=101){
  if(is.null(cols)) cols <- .getDef("hmcols")
  if(is.function(cols)) return(cols)
  if(length(cols) %in% 2:3)  return(colorRampPalette(cols)(n))
  cols
}

.getBaseHMcols <- function(se, cols){
  if(!is.null(cols)) return(cols)
  if(!is.null(se) && !is.null(cols <- metadata(se)$hmcols)) return(cols)
  .getDef("hmcols")
}

#' getBreaks
#'
#' Produces symmetrical breaks for a color scale, with the scale steps
#' increasing for large values, which is useful to avoid outliers influencing
#' too much the color scale.
#'
#' @param x A matrix of log2FC (or any numerical values centered around 0)
#' @param n The desired number of breaks.
#' @param split.prop The proportion of the data points to plot on a linear
#' scale; the remaining will be plotted on a scale with regular frequency per
#' step (quantile).
#' @param symmetric Logical; whether breaks should be symmetric around 0
#'  (default TRUE)
#'
#' @return A vector of breaks of length = `n`
#' @export
#' @importFrom stats quantile
#'
#' @examples
#' dat <- rnorm(100,sd = 10)
#' getBreaks(dat, 10)
getBreaks <- function(x, n, split.prop=0.98, symmetric=TRUE){
  if(is.logical(split.prop)) split.prop <- ifelse(split.prop,0.98,1)
  if(symmetric){
    if((n %% 2) == 0)
      warning("For symmetrical colorscales, an uneven number of colors should",
              " be provided.")
    breaks <- getBreaks(as.numeric(c(0,abs(x))), n=ceiling((n+1)/2),
                        split.prop=split.prop, symmetric=FALSE)
    return(unique(c(-rev(breaks),breaks)))
  }
  minVal <- min(x, na.rm=TRUE)
  if(split.prop<1) q <- as.numeric(quantile(x,split.prop,na.rm=TRUE))
  if(split.prop==1 || !(q>minVal)){
    q <- max(x, na.rm=TRUE)
    split.prop <- 1
  }
  xr <- seq(from=minVal, to=q, length.out=floor(split.prop*n))
  n <- n-length(xr)
  if(n>0){
    q <- quantile(as.numeric(x)[which(x>q)],seq_len(n)/n, na.rm=TRUE)
    xr <- c(xr,as.numeric(q))
  }
  if(any(duplicated(xr))){
    ## duplicated breaks, probably because we have too few datapoints;
    ## we fall back onto a linear scale
    xr <- getBreaks(x, n, 1, symmetric=FALSE)
  }
  xr
}

.getAnnoCols <- function(se, given=list(), do.assign=FALSE){
    ll <- list( default=.getDef("anno_colors") )
    if(!is.null(metadata(se)$anno_colors))
      ll$object <- metadata(se)$anno_colors
    ll$given <- given
    ac <- .mergelists(ll)
    if(do.assign) ac <- .assignAnnoColors(se, ac)
    lapply(ac, unlist)
}

#' @importFrom randomcoloR distinctColorPalette
#' @importFrom SummarizedExperiment colData rowData
.assignAnnoColors <- function(x, anno_colors){
    fn <- function(x){
        if(is.factor(x)) return(levels(x))
        if(is.character(x)) return(unique(x))
        return(NULL)
    }
    if(is(x, "SummarizedExperiment")){
        df <- c( lapply(colData(x), fn), lapply(rowData(x), fn) )
    }else{
        df <- lapply( x, fn)
    }
    for(f in names(df)){
        if(!(f %in% names(anno_colors))) anno_colors[[f]] <- list()
        x <- setdiff(df[[f]], names(anno_colors[[f]]))
        if(length(x)>0) anno_colors[[f]][x] <- distinctColorPalette(length(x))
    }
    anno_colors
}


# non recursive, latest values win
.mergelists <- function(ll){
    names(ll) <- NULL
    names(nn) <- nn <- unique(unlist(lapply(ll,names)))
    lapply(nn, FUN=function(x){
        x <- lapply(ll, function(y) y[[x]])
        x <- x[!unlist(lapply(x,is.null))]
        if(length(x)==0) return(x)
        if(length(x)==1 || is.function(x[[1]])) return(x[[1]])
        x <- do.call(c,x)
        x[!duplicated(names(x))]
    })
}

.has_nan <- function(x){
    if(is(x,"SummarizedExperiment"))
        return(any( unlist(lapply(assays(x), .has_nan)) ))
    any(is.infinite(x) | is.na(x))
}

#' @importFrom grid gpar
#' @importFrom ComplexHeatmap Legend packLegend
.annoLegend <- function(anno_colors){
  leg <- lapply(names(anno_colors), FUN=function(n){
    x <- anno_colors[[n]]
    Legend(labels=names(x), title=n, legend_gp=gpar(fill=as.character(x)))
  })
  packLegend(list=leg)
}

.prepData <- function( se, genes=NULL, do.scale=FALSE,
                       assayName=.getDef("assayName"), includeMissing=FALSE ){
    x <- as.matrix(.chooseAssay(se, assayName))
    if(!is.null(genes)){
        genes <- unique(genes)
        x <- x[intersect(genes,row.names(x)),]
    }
    if(do.scale){
        x <- safescale(x, byRow=TRUE)
    }
    if(includeMissing && length(missg <- setdiff(genes, row.names(x)))>0){
        x2 <- matrix( NA_real_, ncol=ncol(x), nrow=length(missg),
                      dimnames=list(missg, colnames(x)) )
        x <- rbind(x,x2)[genes,]
    }
    as.matrix(x)
}

.parseToporder <- function(x, toporder=NULL){
    if(is(x, "SummarizedExperiment")) x <- rowData(x)
    if(is.null(toporder)) return(NULL)
    if(length(toporder)==1 && is.character(toporder)){
        if(toporder %in% colnames(x)){
            toporder <- x[[toporder]]
            names(toporder) <- row.names(x)
        }else{
            stop("Could not interpret `toporder`.")
        }
    }
    if(!is.null(names(toporder))){
        toporder <- toporder[row.names(x)]
    }else{
        names(toporder) <- row.names(x)
    }
    return(toporder)
}

.prepScale <- function(x, hmcols=NULL, breaks=.getDef("breaks")){
    hmcols <- .getHMcols(cols=hmcols)
    if(!is.null(breaks) && length(breaks)==1 && !is.na(breaks) &&
       (!is.logical(breaks) || breaks)){
      breaks <- getBreaks(x, length(hmcols), split.prop=breaks)
    }else if(is.null(breaks) || all(is.na(breaks)) ||
         (length(breaks)==1 && is.logical(breaks) && !breaks) ){
      breaks <- getBreaks(x, length(hmcols), 1, FALSE)
    }
    if(abs(ld <- length(breaks)-length(hmcols))>0){
      if(ld<0) stop("Too many breaks for the given colors!")
      if(ld>1) warning("Too few breaks - only part of the colors will be used.")
      hmcols <- hmcols[seq(from=1+floor(ld/2), to=length(breaks))]
    }
    list(breaks=breaks, hmcols=hmcols)
}

#' qualitativeColors
#'
#' @param names The names to which the colors are to be assigned, or an integer
#' indicating the desired number of colors
#' @param ... passed to `randomcoloR::distinctColorPalette`
#'
#' @return A vector (eventually named) of colors
#'
#' @importFrom randomcoloR distinctColorPalette
qualitativeColors <- function(names, ...){
    names <- unique(names)
    if(length(names)==1 && is.integer(names))
        return(distinctColorPalette(names))
    cols <- distinctColorPalette(length(names), ...)
    names(cols) <- names
    cols
}


#' safescale
#'
#' Equivalent to `base::scale`, but handling missing values and null variance
#' a bit more elegantly.
#'
#' @param x A matrix.
#' @param center Logical, whether to center values.
#' @param byRow Logical, whether to scale by rows instead of columns.
#'
#' @return A scaled matrix.
#' @export
#'
#' @importFrom matrixStats rowVars rowMaxs
#' @examples
#' m <- matrix(rnorm(100), nrow=10)
#' m.scaled <- safescale(m)
safescale <- function(x, center=TRUE, byRow=FALSE){
  if(!any(is.na(x)) && !byRow) return(base::scale(x, center=center))
  if(is.null(dim(x))) x <- matrix(x)
  if(!byRow) x <- t(x)
  y <- x
  if(center) y <- y - rowMeans(x, na.rm=TRUE)
  rv <- sqrt(matrixStats::rowVars(y, na.rm=TRUE))
  w <- which(rowSums(is.na(x))<ncol(y) & is.na(rv))
  rv[w] <-  matrixStats::rowMaxs(y[w,], na.rm=TRUE)
  y <- y/rv
  y[which(rowSums(is.na(y) | is.infinite(y))==ncol(y)),] <- 0
  if(!byRow) y <- t(y)
  y
}

.getGaps <- function(x, CD, silent=TRUE){
  if(is.null(x)) return(NULL)
  if(is.factor(x) && length(x)>1){
    if(!is.null(names(x))) x <- x[row.names(CD)]
    if(length(x)==nrow(CD)) return(x)
  }
  if(!all(x %in% colnames(CD))){
    if(!silent) warning("Gap field(s) not found in the object's data.")
    return(NULL)
  }
  x <- as.data.frame(CD)[,x,drop=FALSE]
  if(ncol(x)==0) return(NULL)
  x
}

#' @importFrom ComplexHeatmap HeatmapAnnotation
.prepareAnnoDF <- function(an, anno_colors, fields, whichComplex=NULL,
                           show_legend=TRUE, show_annotation_name=TRUE,
                           dropEmptyLevels=TRUE, anno_name_side=NULL,
                           highlight=NULL){
  if(!is.null(whichComplex))
    whichComplex <- match.arg(whichComplex, c("row","column"))
  an <- an[,intersect(fields, colnames(an)),drop=FALSE]
  an <- as.data.frame(an)
  if(ncol(an)==0){
    an <- NULL
  }else{
    for(i in colnames(an)){
      if(is.factor(an[[i]])){
        if(dropEmptyLevels) an[[i]] <- droplevels(an[[i]])
      }
      if(is.logical(an[[i]])){
        an[[i]] <- factor(as.character(an[[i]]),levels=c("FALSE","TRUE"))
        if(!(i %in% names(anno_colors))){
          anno_colors[[i]] <- c("FALSE"="white", "TRUE"="darkblue")
        }
      }else if(!is.null(anno_colors[[i]]) &&
               !is.function(anno_colors[[i]])){
        if(i %in% names(anno_colors)){
          w <- intersect(names(anno_colors[[i]]),unique(an[[i]]))
          if(length(w)==0){
            anno_colors[[i]] <- NULL
          }else{
            anno_colors[[i]] <- anno_colors[[i]][w]
          }
        }
      }
    }
  }
  if(is.null(whichComplex)) return(list(an=an, anno_colors=anno_colors))
  if(is.null(an)) return(NULL)
  anno_colors <- anno_colors[intersect(names(anno_colors),colnames(an))]
  if(length(anno_colors)==0){
    an <- HeatmapAnnotation(df=an, show_legend=show_legend, na_col="white",
                            which=whichComplex,
                            annotation_name_side=anno_name_side,
                            show_annotation_name=show_annotation_name,
                            highlight=highlight )
  }else{
    an <- HeatmapAnnotation(df=an, show_legend=show_legend, na_col="white",
                            which=whichComplex, col=anno_colors,
                            annotation_name_side=anno_name_side,
                            show_annotation_name=show_annotation_name,
                            highlight=highlight)
  }
  an
}


.defaultAnno <- function(se, type="left"){
  if(!is.null(dv <- metadata(se)$default_view) &&
     !is.null(dv <- dv[[paste0(type,"_annotation")]])) return(dv)
  if(type=="top" && !is.null(dv <- metadata(se)$default_view) &&
     sum(lengths(dv <- unique(c(dv$gridvar, dv$groupvar, dv$colvar))))>0)
    return(dv)
  if(!is.null(def <- .getDef(paste0(type,"_annotation")))) return(def)
  # for compatibility with older versions:
  if(type=="left") return(.getDef("anno_rows"))
  if(type=="top") return(.getDef("anno_columns"))
  return(NULL)
}


#' getDEA
#'
#' Extracts (standardized) DEA results from the rowData of an SE object.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}},
#'   with DEAs each saved as a rowData column of `se`, with the column name
#'   prefixed with "DEA."
#' @param dea The optional name of the DEA to extract
#' @param homogenize Logical; whether to homogenize the DEA
#' @param sort Logical; whether to return the table sorted by significance
#'
#' @return The DEA data.frame if `dea` is given, otherwise a named list of
#'   data.frames.
#' @export
#'
#' @examples
#' # loading example SE
#' data("Chen2017", package="sechm")
#' # this ones doesn't have saved DEAs in the standard format:
#' getDEA(Chen2017)
getDEA <- function(se, dea=NULL, homogenize=FALSE, sort=TRUE){
  stopifnot(is(se,"SummarizedExperiment"))
  deas <- grep("^DEA\\.", colnames(rowData(se)), value=TRUE)
  names(deas) <- gsub("^DEA\\.", "", deas)
  deas <- lapply(deas, FUN=function(x){
    x <- rowData(se)[[x]]
    if(homogenize) x <- homogenizeDEA(x)
    if(sort){
      cn <- grep("PValue|P\\.Value|pvalue|p_value|pval", colnames(x), value=TRUE)
      if(length(cn)>0) x <- x[order(x[,cn[1]]),]
    }
    x
  })
  if(!is.null(dea)){
    dea <- gsub("^DEA\\.","",dea)
    if(dea %in% names(deas)) return(deas[[dea]])
    matches <- grep(dea, names(deas), value=TRUE)
    if(length(matches)==0){
      message("No matching DEA found. Available DEAs:\n",
              paste(names(deas), collapse=", "))
    }else if(length(matches)==1){
      message("Returning DEA '", matches,"'")
      if(homogenize) return(homogenizeDEA(deas[[matches]]))
      return(deas[[matches]])
    }else{
      message("Multiple matches:\n",
              paste(matches, collapse=", "))
    }
    message("No DEA found!")
    return(NULL)
  }
  if(length(deas)==0) return(list())
  deas <- deas[!unlist(lapply(deas, is.null))]
  lapply(deas, FUN=function(x){
    x[!is.na(x[,1]),]
  })
}

#' Get DEGs from a SE or list of DEA results
#'
#' @param x A `SummarizedExperiment` object with DEA results in rowData, or a
#'   list of DEA result data.frames.
#' @param dea Which DEA(s) to use (default all). Used only if `x` is a
#'   `SummarizedExperiment`.
#' @param lfc.th Absolute log-foldchange threshold.
#' @param fdr.th FDR threshold.
#' @param direction If !=0, specifies whether to fetch only upregulated or
#'   downregulated features
#' @param merge Logical; whether to take the union of DEGs from the different
#'   DEAs (when more than one).
#'
#' @return A character vector with the significant features, or a list of such
#'   vectors.
#' @export
#'
#' @examples
#' # loading example SE
#' data("Chen2017", package="sechm")
#' # this ones doesn't have saved DEAs in the standard format:
#' getDEGs(Chen2017)
getDEGs <- function(x, dea=NULL, lfc.th=log2(1.3), fdr.th=0.05, direction=0,
                    merge=TRUE){
  if(is(x,"SummarizedExperiment")) x <- getDEA(x, dea=dea)
  stopifnot(is.list(x) || is(x,"DFrame"))
  if(is.data.frame(x) || is(x,"DFrame")){
    x <- x[which(abs(x$logFC)>=lfc.th & x$FDR<=fdr.th),]
    if(direction != 0) x <- x[which(sign(x$logFC)==sign(direction)),]
    return(row.names(x))
  }
  x <- lapply(x, lfc.th=lfc.th, fdr.th=fdr.th, direction=direction, FUN=getDEGs)
  if(length(x)==1) return(x[[1]])
  if(merge) return(unique(unlist(x)))
  return(x)
}


#' homogenizeDEA
#'
#' Standardizes the outputs of differential expression methods (to an
#'   edgeR-like style)
#'
#' @param x A data.frame containing the results of a differential expression
#'   analysis
#'
#' @return A standardized data.frame.
homogenizeDEA <- function(x){
  if(is(x,"data.table") || is(x, "DFrame")) x <- as.data.frame(x)
  if(!is.data.frame(x)) return(NULL)
  colnames(x) <- gsub("log2FoldChange|log2Fold|log2FC|log2\\(fold_change\\)|log2\\.fold_change\\.",
                      "logFC", colnames(x))

  abf <- head(intersect(colnames(x),
                        c("logCPM", "meanExpr", "AveExpr", "baseMean")), 1)
  if (length(abf)==1){
    x$meanExpr <- x[, abf]
    if(abf == "baseMean") x$meanExpr <- log1p(x$meanExpr)
  }else if(all(c("value_1","value_2") %in% colnames(x))){ # cufflinks
    x$meanExpr <- log(1+x$value_1+x$value_2)
  }
  colnames(x) <- gsub("P\\.Value|pvalue|p_value|pval", "PValue", colnames(x))
  colnames(x) <- gsub("padj|adj\\.P\\.Val|q_value|qval", "FDR", colnames(x))
  if (!("FDR" %in% colnames(x)))
    x$FDR <- p.adjust(x$PValue, method = "fdr")
  f <- grep("^logFC$",colnames(x),value=TRUE)
  if(length(f)==0) f <- grep("logFC",colnames(x),value=TRUE)
  if(length(f)==0) warning("No logFC found.")
  if(length(f)>1){
    message("Using ",f[1])
    x[["logFC"]] <- x[[f[1]]]
  }
  x$F <- NULL
  x$FDR[is.na(x$FDR)] <- 1
  x <- x[!is.na(x$logFC),]
  if(!is.null(x$PValue)){
    x <- x[!is.na(x$PValue),]
    x <- x[order(x$PValue),]
  }else{
    x <- x[order(x$FDR),]
  }
  x
}


#' log2FC
#'
#' Generates log2(foldchange) matrix/assay, eventually on a per-batch fashion.
#'
#' @param x A numeric matrix, or a `SummarizedExperiment` object
#' @param fromAssay The assay to use if `x` is a `SummarizedExperiment`
#' @param controls A vector of which samples should be used as controls for
#'   foldchange calculations.
#' @param by An optional vector indicating groups/batches by which the controls
#'   will be averaged to calculate per-group foldchanges.
#' @param isLog Logical; whether the data is log-transformed. If NULL, will
#'   attempt to figure it out from the data and/or assay name
#' @param agFun Aggregation function for the baseline (default rowMeans)
#' @param toAssay The name of the assay in which to save the output. If left to
#'   the default value, both a log2FC assay as well as a scaled log2FC assay
#'   (scaled by unit-variance, but not centered) will be saved in the object.
#' @param pseudocount If the origin assay is not log-transformed, `pseudocount`
#'   will be added to the values before calculating a log-transformation. This
#'   prevents infinite fold-changes and moderates them.
#' @param ndigits Number of digits after the decimal of the log2FC (and
#'   scaledLFC).
#'
#' @return An object of same class as `x`; if a `SummarizedExperiment`, will
#' have the additional assay named from `toAssay`.
#'
#' @examples
#' log2FC( matrix(rnorm(40), ncol=4), controls=1:2 )
#'
#' @import SummarizedExperiment
#' @export
log2FC <- function(x, fromAssay=NULL, controls, by=NULL, isLog=NULL,
                   agFun=rowMeans, toAssay="log2FC", pseudocount=1L, ndigits=2){
  if(is.null(colnames(x))) colnames(x) <- paste0("S",seq_len(ncol(x)))
  if(is(x, "SummarizedExperiment")){
    if(is.null(fromAssay))
      stop("If `x` is a SummarizedExperiment, specify the assay to use ",
           "using `fromAssay`")
    if(!(fromAssay %in% assayNames(x)))
      stop("`fromAssay` '", fromAssay, "' not found.")
    if(!is.null(by) && length(by)==1 && by %in% colnames(colData(x)))
      by <- colData(x)[[by]]
    a <- assays(x)[[fromAssay]]
  }else{
    if(!is.matrix(x))
      stop("`x` should either be a SummarizedExperiment or a numeric matrix.")
    a <- x
  }
  if(is.null(isLog)){
    if(!is.null(fromAssay) && grepl("^log",fromAssay, ignore.case=TRUE)){
      isLog <- TRUE
    }else{
      isLog <- any(a<0)
      if(isLog) message("Assuming assay values to be on a log-scale")
    }
  }
  if(!isLog) log2(a+pseudocount)
  if(is.logical(controls)) controls <- which(controls)
  if(!all(controls %in% seq_len(ncol(a))))
    stop("Some control indexes are out of range.")
  if(is.null(by)) by <- rep(1,ncol(a))
  i <- split(1:ncol(a),by)
  lfc <- do.call(cbind, lapply(i, FUN=function(x){
    c2 <- intersect(x,controls)
    if(length(c2)==0) stop("Some groups of `by` have no controls.")
    a[,x,drop=FALSE]-agFun(a[,c2,drop=FALSE],na.rm=TRUE)
  }))
  lfc <- lfc[,colnames(x)]
  if(!is.null(ndigits)) lfc <- round(lfc, ndigits)
  if(is(x, "SummarizedExperiment")){
    assays(x)[[toAssay]] <- lfc
    if(toAssay=="log2FC"){
      slfc <- safescale(assays(x)$log2FC, center=FALSE, byRow=TRUE)
      if(!is.null(ndigits)) slfc <- round(slfc, ndigits)
      assays(x)$scaledLFC <- slfc
    }
    return(x)
  }
  lfc
}

#' Set rowData attribute of given rows
#'
#' @param se A `SummarizedExperiment` object
#' @param values A named vector of values, where the names correspond to rows of
#'   `se`
#' @param name The name of the rowData column in which to store the attribute.
#' @param clear Logical; whether to clear out any pre-existing such column.
#' @param other The value for unspecified rows (default NA)
#'
#' @return The modified `se` object.
#' @export
#'
#' @examples
#' data("Chen2017", package="sechm")
#' Chen2017 <- setRowAttr(Chen2017, c("Arc"=1,"Junb"=1,"Npas4"=2))
setRowAttr <- function(se, values, name="cluster", clear=TRUE, other=NA){
  if(is.data.frame(values)){
    stopifnot(!is.null(row.names(values)))
    values <- lapply(values, FUN=function(x) setNames(x,row.names(values)))
    for(v in names(values)){
      se <- setRowAttr(se, values[[v]], name=v, clear=clear)
    }
    return(se)
  }
  stopifnot(is.vector(values))
  stopifnot(!is.null(names(values)))
  stopifnot(any(names(values) %in% row.names(se)))
  if(clear || is.null(rowData(se)[[name]])) rowData(se)[[name]] <- other
  rowData(se)[names(values),name] <- as.vector(values)
  se
}
