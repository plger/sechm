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
  if(z) y <- t(.safescale(t(x)))
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
  if(!is.null(dan <- metadata(se)$default_view$assay) && 
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

  if(length(assayName)>1) message("Using assay ", assayName[1])
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
        #x <- x[apply(x,1,FUN=sd)>0,]
        x <- t(.safescale(t(x)))
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


.safescale <- function(x){
    if(!any(is.na(x))) base::scale(x)
    if(!is.null(dim(x))){
        y <- apply(x,2,.safescale)
        row.names(y) <- row.names(x)
        return(y)
    }
    if(all(is.na(x))) return(x)
    if(sd(x,na.rm=TRUE)>0) return(base::scale(x))
    if(sum(!is.na(x))==0) return(base::scale(as.numeric(!is.na(x))))
    rep(0,length(x))
}


#' scale2
#'
#' A wrapper for non-centered unit-variance scaling
#' @param x A matrix whose rows are to be scaled.
#'
#' @return A matrix of dimensions like x
#' @export
#'
#' @examples
#' scale2(matrix(1:9,nrow=3))
scale2 <- function(x){
  y <- t(scale(t(x),center=FALSE))
  y[is.nan(y)] <- 0
  y
}

.getGaps <- function(x, CD, silent=TRUE){
  if(is.null(x)) return(NULL)
  if(!all(x %in% colnames(CD))){
    if(!silent) warning("Gap field(s) not found in the object's data.")
    return(NULL)
  }
  as.data.frame(CD)[,x,drop=FALSE]
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