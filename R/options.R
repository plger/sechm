.options <- local({
  options <- list("assayName"=c("logFC", "log2FC", "logcpm", "lognorm"),
                  "anno_colors"=list(), hmcols=c("blue","black","yellow"),
                  "top_annotation"=c(
                    "Batch", "batch", "Condition","condition", "Group", "group",
                    "Dataset", "Genotype", "genotype", "cluster_id", "group_id",
                    "celltype"), gaps_at=c("Dataset","cluster_id"),
                  breaks=NULL)
  names(emptyAnnos) <- emptyAnnos <- c("bottom_annotation", "left_annotation",
                                       "right_annotation", "anno_rows",
                                       "anno_columns")
  options <- c(options, lapply(emptyAnnos, FUN=function(x) character(0)))

  env <- new.env(parent=emptyenv())

  list(set=function(variable, value) {
    stopifnot(is.character(variable), length(variable) == 1L, !is.na(variable))
    stopifnot(variable %in% names(options))
    env[[variable]] <- value
  }, get=function(variable, default=NULL) {
    if(variable %in% names(env)) return(env[[variable]])
    if(variable %in% names(options)) return(options[[variable]])
    warning("Unknown option ", variable)
    return(NULL)
  }, reset=function(){
    for(f in names(options)) env[[f]] <- options[[f]]
  })
})

.getDef <- function(x, ...) .options$get(x)

#' resetAllSechmOptions
#'
#' Resents all package options
#'
#' @return None
#'
#' @examples
#' resetAllSechmOptions()
#'
#' @export
resetAllSechmOptions <- function() .options$reset()

#' setSechmOption
#'
#' Sets a package-wide option for `sechm`
#'
#' @param variable The name of the variable to set
#' @param value The parameter value to save
#'
#' @return None
#'
#' @examples
#' setSechmOption("hmcols", value=c("blue","black","yellow"))
#'
#' @export
setSechmOption <- function(variable, value){
  .options$set(variable,value)
}