#' crossHm
#'
#' Plot a multi-panel heatmap from a list of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#' @param ses A (named) list of
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} objects,
#'  with some matching row.names between them.
#' @param features A vector of features (i.e. row.names) to plot.
#' @param assayName The name of the assay to use; if multiple names are given,
#' the first available will be used. Defaults to "logcpm", "lognorm".
#' @param only.common Logical; whether to plot only rows common to all SEs
#' (default TRUE).
#' @param do.scale Logical; whether to scale rows in each SE (default TRUE).
#' @param uniqueScale Logical; whether to force the same colorscale for
#' each heatmap.
#' @param sortBy Names or indexes of `ses` to use for sorting rows (default all)
#' @param cluster_cols Logical; whether to cluster columns (default FALSE).
#' @param cluster_rows Logical; whether to cluster rows (default TRUE if
#' `do.sortRows=FALSE`, FALSE otherwise).
#' @param toporder Optional verctor of categories on which to supra-order when
#' sorting rows, or name of a `rowData` column to use for this purpose.
#' @param hmcols Colors for the heatmap.
#' @param breaks Breaks for the heatmap colors. Alternatively, symmetrical
#' breaks can be generated automatically by setting `breaks` to a numerical
#' value between 0 and 1. The value is passed as the `split.prop` argument to
#' the \code{\link{getBreaks}} function, and indicates the proportion of the
#' points to map to a linear scale, while the more extreme values will be
#' plotted on a quantile scale. `breaks=FALSE` will disable symmetrical scale
#' and quantile capping, while retaining automatic breaks. `breaks=1` will
#' produce a symmetrical scale without quantile capping.
#' @param gaps_at Columns of `colData` to use to establish gaps between columns.
#' @param gaps_row A named vector according to which rows will be split.
#' @param left_annotation Columns of `rowData` to use for left annotation.
#' @param top_annotation Columns of `colData` to use for top annotation.
#' @param name The title of the heatmap key.
#' @param anno_colors List of colors to use for annotation.
#' @param show_rownames Whether to show row names (default TRUE if 50 rows or
#' less).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param rel.width Relative width of the heatmaps
#' @param merge_legends Logical; passed to
#' \code{\link[ComplexHeatmap]{draw-HeatmapList-method}}
#' @param ... Any other parameter passed to each call of
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A Heatmap list.
#'
#' @examples
#' data("Chen2017", package="sechm")
#' se1 <- Chen2017[,1:6]
#' se2 <- Chen2017[,7:15]
#' se3 <- crossHm(list(se1=se1, se2=se2), row.names(se1)[1:10] )
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @importFrom S4Vectors metadata metadata<- SimpleList
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
crossHm <- function(ses, features, do.scale=TRUE, uniqueScale=FALSE,
                    assayName=.getDef("assayName"), sortBy=seq_along(ses),
                    only.common=TRUE, cluster_cols=FALSE,
                    cluster_rows=is.null(sortBy), toporder=NULL, hmcols=NULL,
                    breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                    gaps_row=NULL,  name=NULL,
                    top_annotation=.getDef("anno_columns"), 
                    left_annotation=.getDef("anno_rows"), 
                    anno_colors=list(), show_rownames=NULL, merge_legends=FALSE,
                    show_colnames=FALSE, rel.width=NULL, ... ){

  stopifnot(is.vector(features))
  if(is.factor(features)) features <- as.character(features)
  stopifnot(is.character(features))
  
  if(is(ses,"SummarizedExperiment")) ses <- list(ses)
  if(is.null(names(ses))) names(ses) <- paste("SE", seq_along(ses))
  if(!is.null(rel.width) && length(rel.width)!=length(ses))
      stop("If given, `rel.width` should have the same length as `ses`.")
  if(is.null(rel.width)) rel.width <- rep(1,length(ses))

  tt <- table(unlist(lapply(ses,row.names)))
  if(only.common){
      features <- intersect(features, names(tt)[which(tt==length(ses))])
  }else{
    features <- intersect(features, names(tt))
  }
  if(length(features)==0)
      stop("There appears to be not feature in common across `ses`")
  if(length(features)<=2){
    sortBy <- NULL
    cluster_rows <- FALSE
  }
  if(!is.null(toporder)){
      if(is.null(names(toporder)))
          stop("`toporder` should be a vector named by feature")
      if(!all(features %in% names(toporder)))
          warning("Some features are missing from `toporder`")
      toporder <- toporder[features]
      names(toporder) <- features
  }

  dats <- lapply( ses, FUN=.prepData, genes=features, assayName=assayName,
                  do.scale=do.scale && !uniqueScale, includeMissing=TRUE )
  if(do.scale && uniqueScale){
      x <- do.call(cbind, dats)
      x <- t(.safescale(t(x)))
      dl <- lapply(dats, FUN=function(x) seq_len(ncol(x)))
      dl2 <- c(0,cumsum(lengths(dl)[-length(dl)]))
      dats <- lapply( seq_along(dl), FUN=function(i)
          x[,dl[[i]]+dl2[[i]],drop=FALSE] )
      names(dats) <- names(dl)
  }else{
      x <- NULL
  }

  if(uniqueScale){
      if(is.null(x)) x <- do.call(cbind, dats)
      if(is.null(breaks) && do.scale) breaks <- 0.995
      cscale <- .prepScale(x, hmcols=.getHMcols(hmcols), breaks=breaks)
      breaks <- cscale$breaks
      hmcols <- cscale$hmcols
  }

  if(!is.null(sortBy) && length(sortBy)>0){
      xs <- dats
      if(do.scale && !uniqueScale)
          xs <- lapply(xs, FUN=function(x){ t(.safescale(t(x))) })
      xs <- do.call(cbind, xs[sortBy])
      features <- row.names(sortRows(xs,toporder=toporder))
      dats <- lapply(dats, FUN=function(x) x[features,,drop=FALSE])
  }

  CDs <- lapply(ses, ac=top_annotation, FUN=function(x,ac){
      ac <- intersect(colnames(colData(x)),ac)
      as.data.frame(colData(x)[,ac,drop=FALSE])
  })

  # make sure factors share the levels across datasets
  facts <- unique(unlist(lapply(CDs, FUN=function(x){
      x <- vapply(x, class, character(1))
      names(x)[x=="factor"]
  })))
  for(v in facts){
      lvls <- unique(unlist(lapply(CDs, FUN=function(x){
          if(v %in% colnames(x)) return(levels(droplevels(as.factor(x[[v]]))))
          NULL
      })))
      CDs <- lapply(CDs, FUN=function(x){
          if(!(v %in% colnames(x))) return(x)
          x[[v]] <- factor(as.character(x[[v]]), levels=lvls)
          x
      })
  }

  ses <- lapply(seq_along(ses), FUN=function(i){
      RD <- rowData(ses[[i]])[features,,drop=FALSE]
      se <- SummarizedExperiment( list(a=dats[[i]]), colData=CDs[[i]],
                                  rowData=RD, metadata=metadata(ses[[i]]) )
      row.names(se) <- features
      se
  })
  names(ses) <- names(dats)
  tmp <- unlist(lapply(ses, FUN=function(x) colnames(rowData(x))))
  afields <- intersect(c(top_annotation, left_annotation),
                       c(tmp, colnames(CDs[[1]])))
  anno_colors <- .getAnnoCols(ses[[1]], anno_colors, do.assign=TRUE)
  anno_colors <- anno_colors[afields]

  if(is.null(show_rownames)) show_rownames <- length(features)<50

  hlp <- list()
  if(uniqueScale){
      if(!is.null(assayName) && length(assayName)==1 &&
         !is.numeric(assayName)){
          hlp$title <- ifelse(do.scale, paste0("scaled\n",assayName),
                              assayName)
      }else{
          hlp$title <- ifelse(do.scale, "z-scores", "")
      }
  }

  htlist <- lapply(seq_along(ses), FUN=function(i){
      sechm(ses[[i]], features=features, do.scale=(do.scale && !uniqueScale),
            assayName="a", name=names(ses)[i], toporder=toporder,
            hmcols=hmcols, breaks=breaks, left_annotation=left_annotation,
            top_annotation=top_annotation, anno_colors=anno_colors,
            cluster_rows=cluster_rows, cluster_cols=cluster_cols,
            show_rownames=(show_rownames && i==length(ses)), sortRowsOn=NULL,
            show_colnames=show_colnames, isMult=i!=length(ses),
            show_heatmap_legend=(!uniqueScale || i==length(ses)),
            heatmap_legend_param=hlp, column_title=names(ses)[i],
            show_annotation_legend=FALSE, includeMissing=!only.common,
            width=rel.width[i], ...)
  })

  ht <- NULL
  for(f in htlist) ht <- ht + f
  if(length(anno_colors)>0)
      return(ComplexHeatmap::draw(ht,
              annotation_legend_list=.annoLegend(anno_colors),
              merge_legends=merge_legends))
  ht
}
