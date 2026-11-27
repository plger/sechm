#' commonScale
#'
#' Computes a single color scale encompassing several heatmaps, so that these
#' can be re-generated with exactly the same color mapping (and a single legend
#' shown).
#'
#' @param ... Any number of \code{\link[ComplexHeatmap]{Heatmap-class}} or
#'   \code{\link[ComplexHeatmap]{HeatmapList-class}} objects (e.g. as produced
#'   by \code{\link{sechm}}), or lists thereof.
#' @param hmcols Colors to use. If missing, the colors of the first heatmap are
#'   used.
#' @param breaks Either a number between 0 and 1 (symmetrical scale, with the
#'   given proportion of the points on a linear scale and the rest on a
#'   quantile scale), `1` (symmetrical, no quantile capping), `FALSE` (linear,
#'   non-symmetrical), or a vector of breaks, which is then used as is. If
#'   `NULL` (default), the scale will be symmetrical (with `0.995`) if all the
#'   heatmaps' original scales were symmetrical around 0, and linear otherwise.
#'
#' @return A list with the elements `breaks` and `hmcols`, which can be passed
#'   as the arguments of the same names to \code{\link{sechm}} (to re-generate
#'   the heatmaps with identical scales).
#'
#' @examples
#' data("Chen2017", package="sechm")
#' g <- row.names(Chen2017)[1:10]
#' h1 <- sechm(Chen2017[,1:6], g, do.scale=TRUE)
#' h2 <- sechm(Chen2017[,7:15], g, do.scale=TRUE)
#' cs <- commonScale(h1, h2)
#' h1 <- sechm(Chen2017[,1:6], g, do.scale=TRUE, breaks=cs$breaks,
#'             hmcols=cs$hmcols)
#' h2 <- sechm(Chen2017[,7:15], g, do.scale=TRUE, breaks=cs$breaks,
#'             hmcols=cs$hmcols, show_heatmap_legend=FALSE)
#' h1 + h2
#'
#' @importFrom circlize colorRamp2
#' @export
commonScale <- function(..., hmcols=NULL, breaks=NULL){
  hts <- .flattenHeatmaps(list(...))
  hts <- Filter(function(x) is.numeric(x@matrix), hts)
  if(length(hts)==0) stop("No heatmap with a continuous color scale found.")
  x <- unlist(lapply(hts, function(x) as.numeric(x@matrix)))
  x <- x[is.finite(x)]
  if(length(x)==0) stop("No finite values in the heatmaps.")
  cf <- lapply(hts, function(x) x@matrix_color_mapping@col_fun)
  ob <- lapply(cf, attr, which="breaks")
  oc <- lapply(cf, attr, which="colors")
  if(is.null(hmcols)){
    if(is.null(oc[[1]])) stop("Could not retrieve colors; please set `hmcols`.")
    hmcols <- oc[[1]]
  }
  hmcols <- .getHMcols(hmcols)
  if(is.null(breaks)){
    sym <- all(vapply(ob, FUN.VALUE=logical(1), FUN=function(b){
      !is.null(b) && isTRUE(all.equal(b, -rev(b), check.attributes=FALSE))
    }))
    breaks <- if(sym) 0.995 else FALSE
  }
  cs <- .prepScale(x, hmcols=hmcols, breaks=breaks)
  cs$col_fun <- colorRamp2(cs$breaks, cs$hmcols)
  cs
}

.flattenHeatmaps <- function(x){
  if(is(x, "Heatmap")) return(list(x))
  if(is(x, "HeatmapList")) return(.flattenHeatmaps(x@ht_list))
  if(is.list(x)) return(unlist(lapply(x, .flattenHeatmaps), recursive=FALSE))
  stop("Expected Heatmap or HeatmapList objects.")
}
