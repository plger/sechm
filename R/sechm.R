#' sechm
#'
#' ComplexHeatmap wrapper for
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @param features A vector of features (i.e. row names of `se`). Alternatively,
#'   can be a list of feature sets, in which case these will be plotted as
#'   different row chunks.
#' @param do.scale Logical; whether to scale rows (default FALSE).
#' @param assayName An optional vector of assayNames to use. The first available
#'  will be used, or the first assay if NULL.
#' @param cluster_cols Whether to cluster columns (default F)
#' @param cluster_rows Whether to cluster rows; default FALSE if
#' `do.sortRows=TRUE`.
#' @param sortRowsOn Sort rows by MDS polar order using the specified columns
#' (default all)
#' @param toporder Optional vector of categories on which to supra-order when
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
#' @param gaps_row Passed to the heatmap function; if missing, will
#' be set automatically according to toporder.
#' @param left_annotation Columns of `rowData` to use for left annotation.
#' Alternatively, an `HeatmapAnnotation` object.
#' @param right_annotation Columns of `rowData` to use for left annotation.
#' Alternatively, an `HeatmapAnnotation` object.
#' @param top_annotation Columns of `colData` to use for top annotation.
#' Alternatively, an `HeatmapAnnotation` object. To disable (overriding
#'   defaults), use `top_annotation=character()`.
#' @param bottom_annotation Columns of `colData` to use for bottom annotation.
#' Alternatively, an `HeatmapAnnotation` object.
#' @param anno_colors List of colors to use for annotation.
#' @param annorow_title_side Side (top or bottom) of row annotation names
#' @param annocol_title_side Side (left or right) of column annotation names
#' @param name The name of the heatmap, eventually appearing as title of the
#' color scale.
#' @param show_rownames Whether to show row names (default TRUE if less than
#' 50 rows to plot).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param show_annotation_legend Logical; whether to show the annotation legend.
#' @param includeMissing Logical; whether to include missing features (default
#' FALSE)
#' @param mark An optional vector of gene names to highlight.
#' @param isMult Logical; used to silence labels when plotting multiple heatmaps
#' @param show_heatmap_legend Logical; whether to show heatmap legend
#' @param sort.method Row sorting method (see \code{\link{sortRows}})
#' @param na_col Color of NA values
#' @param ... Further arguments passed to `Heatmap`
#'
#' @return A a \code{\link[ComplexHeatmap]{Heatmap-class}}.
#'
#' @examples
#' data("Chen2017", package="sechm")
#' sechm(Chen2017, row.names(Chen2017)[1:10], do.scale=TRUE)
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
#' @name sechm
#' @rdname sechm
sechm <- function(se, features, do.scale=FALSE, assayName=NULL, name=NULL,
                  sortRowsOn=NULL, cluster_cols=FALSE, cluster_rows=NULL, 
                  toporder=NULL, hmcols=NULL, breaks=.getDef("breaks"), 
                  gaps_at=NULL, gaps_row=NULL, left_annotation=NULL,
                  right_annotation=NULL, top_annotation=NULL,
                  bottom_annotation=NULL, anno_colors=list(),
                  show_rownames=NULL, show_colnames=FALSE,
                  isMult=FALSE, show_heatmap_legend=!isMult,
                  show_annotation_legend=TRUE, mark=NULL, na_col="white",
                  annorow_title_side=ifelse(show_colnames,"bottom","top"),
                  annocol_title_side="right", includeMissing=FALSE,
                  sort.method="MDS_angle", ...){

  if(is.list(features)){
    if(!is.null(gaps_row))
      warning("`gaps_row` ignored when passing a list of features")
    gaps_row <- rep(factor(names(features),names(features)),lengths(features))
    names(gaps_row) <- features <- unlist(features)
  }
  stopifnot(is.vector(features))
  if(is.numeric(features) && all(features==round(features)) && all(features>=1))
      features <- row.names(se)[features]
  if(is.factor(features)) features <- as.character(features)
  stopifnot(is.character(features))
  if(length(intersect(features, row.names(se)))==0)
    stop("None of the features were found in the object.")
  
  assayName <- .chooseAssay(se, assayName, returnName = TRUE)
  if(is.null(name)){
      if(is.numeric(assayName)){
          name <- ifelse(do.scale, "z-scores", "assay")
      }else{
          name <- ifelse(do.scale, paste0("scaled\n",assayName), assayName)
      }
  }

  x <- .prepData(se, genes=features, do.scale=do.scale, assayName=assayName,
                 includeMissing=includeMissing )
  if(nrow(x)<3){
    sortRowsOn <- cluster_rows <- FALSE
  }
  
  if(is.null(sortRowsOn) && is.null(cluster_rows)){
    if(nrow(x)>1500){
      message("Row sorting/clustering disable by default due to heatmap size.")
      sortRowsOn <- cluster_rows <- FALSE
    }else{
      sortRowsOn <- seq_len(ncol(se))
    }
  }
  if(!is.null(sortRowsOn)) cluster_rows <- FALSE

  toporder <- .parseToporder(rowData(se)[row.names(x),,drop=FALSE], toporder)
  row_order <- NULL
  if(!is.null(sortRowsOn) && !isFALSE(sortRowsOn) && length(sortRowsOn)>0 && 
     nrow(x)>2){
    x2 <- sortRows(x[,sortRowsOn,drop=FALSE], toporder=toporder,
                   na.rm=TRUE, method=sort.method)
    row_order <- row.names(x2)
  }

  if( is.null(breaks) ){
      if( (!is.null(assayName) && grepl("^log[2]?FC$|scaledLFC",assayName))
          || do.scale)  breaks <- 0.995
  }
  hmcols <- .getBaseHMcols(se, hmcols)
  cscale <- .prepScale(x, hmcols=hmcols, breaks=breaks)

  breaks <- cscale$breaks
  hmcols <- circlize::colorRamp2(breaks, cscale$hmcols)

  anno_colors <- .getAnnoCols(se, anno_colors)

  if(!is.null(mark) && length(mark <- which(row.names(x) %in% mark))>0){
    mark <- anno_mark(mark, row.names(x)[mark], which="row")
  }else{
    mark <- NULL
  }

  if(is.null(left_annotation)) left_annotation <- .defaultAnno(se, "left")
  if(is.null(right_annotation) && is.null(mark))
    right_annotation <- .defaultAnno(se, "right")
  if(is.null(top_annotation)) top_annotation <- .defaultAnno(se, "top")
  if(is.null(bottom_annotation)) bottom_annotation <- .defaultAnno(se, "bottom")

  if(!is(left_annotation,"HeatmapAnnotation")){
    if(length(left_annotation)>0 && is.character(left_annotation)){
      left_annotation <- .prepareAnnoDF(
        rowData(se)[row.names(x),,drop=FALSE], anno_colors,
        left_annotation, whichComplex="row", show_legend=show_annotation_legend,
        show_annotation_name=!is.na(annorow_title_side),
        anno_name_side=ifelse(is.na(annorow_title_side), "top",
                              annorow_title_side))
    }else{
      left_annotation <- NULL
    }
  }

  if(!is(right_annotation,"HeatmapAnnotation")){
    if(length(right_annotation)>0 && is.character(right_annotation)){
      right_annotation <- .prepareAnnoDF(
        rowData(se)[row.names(x),,drop=FALSE], anno_colors,
        right_annotation, whichComplex="row",
        show_legend=show_annotation_legend,
        show_annotation_name=!is.na(annorow_title_side), highlight=mark)
    }else{
      right_annotation <- NULL
    }
  }

  if(is.null(right_annotation) && !is.null(mark)){
    right_annotation <- rowAnnotation(highlight=mark)
  }

  if(!is(top_annotation,"HeatmapAnnotation")){
    if(length(top_annotation)>0 && is.character(top_annotation)){
      top_annotation <- .prepareAnnoDF(
        colData(se), anno_colors, top_annotation, whichComplex="column",
        show_legend=(show_annotation_legend && !isMult),
        show_annotation_name=(!is.na(annocol_title_side) && !isMult),
        anno_name_side=ifelse(is.na(annocol_title_side), "right",
                              annocol_title_side) )
    }else{
      top_annotation <- NULL
    }
  }
  if(!is(bottom_annotation,"HeatmapAnnotation")){
    if(length(bottom_annotation)>0 && is.character(bottom_annotation)){
      bottom_annotation <- .prepareAnnoDF(
        colData(se), anno_colors, bottom_annotation, whichComplex="column",
        show_legend=(show_annotation_legend && !isMult),
        show_annotation_name=(!is.na(annocol_title_side) && !isMult),
        anno_name_side=ifelse(is.na(annocol_title_side), "right",
                              annocol_title_side) )
    }else{
      bottom_annotation <- NULL
    }
  }

  if(is.null(gaps_at) && !is.null(metadata(se)$default_view))
    gaps_at <- metadata(se)$default_view$gridvar
  if(is.null(gaps_at)) gaps_at <- .getDef("gaps_at")
  if(!is.null(gaps_at) && length(gaps_at)!=ncol(se))
    gaps_at <- intersect(gaps_at, colnames(colData(se)))
  if(isFALSE(gaps_at)) gaps_at <- NULL
  gaps_col <- .getGaps(gaps_at, colData(se), silent=TRUE)
  gaps_row <- .getGaps(gaps_row, rowData(se)[row.names(x),,drop=FALSE])

  if(is.null(show_rownames)) show_rownames <- nrow(x)<50 && is.null(mark)
  if(nrow(x)<=2) cluster_rows <- FALSE

  Heatmap( x, col=hmcols, na_col=na_col, name=name, row_order=row_order,
    show_row_names=show_rownames, show_column_names=show_colnames,
    row_split=gaps_row, column_split=gaps_col, cluster_rows=cluster_rows,
    show_heatmap_legend=show_heatmap_legend, cluster_columns=cluster_cols,
    top_annotation=top_annotation, left_annotation=left_annotation,
    bottom_annotation=bottom_annotation, right_annotation=right_annotation,
    ...)
}
