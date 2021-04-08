#' Heatmap wrapper for
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @name SE-heatmap
#' @rdname SE-heatmap
#' @aliases sehm sechm
#'
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#' @param genes An optional vector of genes (i.e. row names of `se`)
#' @param do.scale Logical; whether to scale rows (default FALSE).
#' @param assayName An optional vector of assayNames to use. The first available
#'  will be used, or the first assay if NULL.
#' @param sortRowsOn Sort rows by MDS polar order using the specified columns
#' (default all)
#' @param cluster_cols Whether to cluster columns (default F)
#' @param cluster_rows Whether to cluster rows; default FALSE if
#' `do.sortRows=TRUE`.
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
#' @param gaps_row Passed to the heatmap function; if missing, will
#' be set automatically according to toporder.
#' @param anno_rows Columns of `rowData` to use for annotation.
#' @param anno_columns Columns of `colData` to use for annotation.
#' @param anno_colors List of colors to use for annotation.
#' @param anno_rows_title_side Side (top or bottom) of row annotation names
#' @param name The name of the heatmap, eventually appearing as title of the
#' color scale.
#' @param show_rownames Whether to show row names (default TRUE if 50 rows or
#' less).
#' @param show_colnames Whether to show column names (default FALSE).
#' @param show_annotation_legend Logical; whether to show the annotation legend.
#' @param includeMissing Logical; whether to include missing genes (default
#' FALSE)
#' @param mark An optional vector of gene names to highlight.
#' @param right_annotation Passed to `ComplexHeatmap::Heatmap`
#' @param ... Further arguments passed to `pheatmap` (`sehm`) or `Heatmap`
#' (`sechm`).
#'
#' @return A a \code{\link[ComplexHeatmap]{Heatmap-class}}.
#'
#' @examples
#' data("SE", package="sechm")
#' sehm(SE, row.names(SE)[1:10], do.scale=TRUE)
#'
#' @param isMult Logical; used to silence labels when plotting mulitple heatmaps
#' @param show_heatmap_legend Logical; whether to show heatmap legend
#'
#' @importFrom circlize colorRamp2
#' @importFrom methods is
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @export
sechm <- function(se, genes, do.scale=FALSE, assayName=.getDef("assayName"),
                  sortRowsOn=seq_len(ncol(se)), cluster_cols=FALSE,
                  cluster_rows=is.null(sortRowsOn), toporder=NULL, hmcols=NULL,
                  breaks=.getDef("breaks"), gaps_at=.getDef("gaps_at"),
                  gaps_row=NULL, anno_rows=.getDef("anno_rows"),
                  anno_columns=.getDef("anno_columns"), name=NULL,
                  anno_colors=list(), show_rownames=NULL, show_colnames=FALSE,
                  isMult=FALSE, show_heatmap_legend=!isMult,
                  show_annotation_legend=TRUE, mark=NULL,
                  annorow_title_side=ifelse(show_colnames,"bottom","top"),
                  right_annotation=NULL, includeMissing=FALSE,
                  sort.method="MDS_angle", ...){

  assayName <- .chooseAssay(se, assayName, returnName = TRUE)
  if(is.null(name)){
      if(is.numeric(assayName)){
          name <- ifelse(do.scale, "z-scores", "assay")
      }else{
          name <- ifelse(do.scale, paste0("scaled\n",assayName), assayName)
      }
  }

  x <- .prepData(se, genes=genes, do.scale=do.scale, assayName=assayName,
                 includeMissing=includeMissing )

  toporder <- .parseToporder(rowData(se)[row.names(x),,drop=FALSE], toporder)
  if(!is.null(sortRowsOn) && length(sortRowsOn)>0 && nrow(x)>2){
      x2 <- sortRows(x[,sortRowsOn,drop=FALSE], toporder=toporder,
                     na.rm=TRUE, method=sort.method)
      x <- x[row.names(x2),]
  }

  if( is.null(breaks) ){
      if( (!is.null(assayName) && grepl("^log[2]?FC$",assayName)) || do.scale)
          breaks <- 0.995
  }
  hmcols <- .getBaseHMcols(se, hmcols)
  cscale <- .prepScale(x, hmcols=hmcols, breaks=breaks)

  breaks <- cscale$breaks
  hmcols <- circlize::colorRamp2(breaks, cscale$hmcols)


  anno_colors <- .getAnnoCols(se, anno_colors)

  anr <- .prepareAnnoDF( rowData(se)[row.names(x),,drop=FALSE], anno_colors,
                         anno_rows, whichComplex="row",
                         show_legend=show_annotation_legend,
                         show_annotation_name=!is.na(annorow_title_side),
                         anno_name_side=ifelse(is.na(annorow_title_side), "top",
                                               annorow_title_side))

  an <- .prepareAnnoDF( colData(se), anno_colors, anno_columns,
                        whichComplex="column",
                        show_legend=(show_annotation_legend && !isMult),
                        show_annotation_name=!isMult,
                        anno_name_side="right" )

  gaps_col <- .getGaps(gaps_at, colData(se), silent=TRUE)
  gaps_row <- .getGaps(gaps_row, rowData(se)[row.names(x),,drop=FALSE])

  if(is.null(show_rownames)) show_rownames <- nrow(x)<50 && is.null(mark)
  if(nrow(x)<=2) cluster_rows <- FALSE
  if(!is.null(mark) && is.null(right_annotation)){
      mark <- which(row.names(x) %in% mark)
      right_annotation <- rowAnnotation(mark=anno_mark(mark,
                                                       row.names(x)[mark]))
  }
  Heatmap(x, col=hmcols, na_col="white", name=name,
          show_row_names=show_rownames, show_column_names=show_colnames,
          top_annotation=an, left_annotation=anr, row_split=gaps_row,
          column_split=gaps_col, show_heatmap_legend=show_heatmap_legend,
          cluster_rows=cluster_rows, cluster_columns=cluster_cols,
          right_annotation=right_annotation, ...)
}
