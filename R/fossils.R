#' Fossils object
#'
#' Create a fossil record object. The input is taken to be a dataframe or list.
#'
#' The fossil record object contains 4 fields for each fossil with the following information:
#' \itemize{
#'  \item \code{sp} the label of the corresponding species. This label matches the edge labels in the corresponding phylo object or the species labels in the corresponding taxonomy object if additional taxonomic
#'  information was provided
#'  \item \code{edge} the label of the sampled node or tip in the phylogeny, i.e the node at the end of the edge along which the fossil was sampled
#'  \item \code{hmin} the age of the fossil or the youngest bound of the time interval in which the fossil was sampled
#'  \item \code{hmax} the oldest bound of the time interval in which the fossil was sampled.
#'  This is equal to \code{hmin} if exact sampling times are known
#' }
#'
#' @param data Dataframe or list of sampled fossils. See Details for the list of required fields. If NULL, the function creates an empty fossils object.
#' @param from.taxonomy Boolean indicating whether the fossils were sampled using a taxonomy object, as opposed to a tree object. Default = FALSE.
#'
#' @export
fossils <- function(data = NULL, from.taxonomy = FALSE){
  if(is.null(data)) {
    data = data.frame(sp = numeric(), edge = numeric(), hmin = numeric(), hmax = numeric(), stringsAsFactors = F)
  }
  else {
    if(is.list(data)) data <- as.data.frame(data)

    # check for required fields
    required_fields = c("sp", "edge", "hmin", "hmax")
    missing = !required_fields %in% colnames(data)
    if(any(missing)) stop(paste0("Missing required fields: ", paste(required_fields[missing], collapse = ", ")))
    
    # check for mistyped fields (e.g. species as string rather than index)
    mistyped = sapply(required_fields, function(f) mode(data[[f]]) != "numeric")
    if(any(mistyped)) stop(paste0("All required fields should be of numeric type. Mistyped fields: ", paste(colnames(data)[mistyped], collapse = ", ")))

    # check for additional fields
    additional = !colnames(data) %in% required_fields
    if(any(additional)) {
      warning(paste0("These fields will be discarded: ", paste(colnames(data)[additional], collapse = ", ")))
      data[,additional] = NULL
    }
  }

  attr(data, "class") <- c("fossils", class(data))
  attr(data, "from.taxonomy") <- from.taxonomy

  return(data)
}

#' @export
#' @aliases fossils
print.fossils <- function(x, max.length = 10, ...){
  if(length(x$sp) > 0){
    if(length(x$sp) < max.length)
      max.length = length(x$sp)
    print(as.data.frame(x)[1:max.length,])
    if(length(x$sp) > max.length)
      cat("...\n")
  }
  cat("Fossil record with", length(x$sp), "occurrences representing", length(unique(x$sp)), "species\n")
  if(!is.null(attr(x,"from.taxonomy"))){
    if(attr(x,"from.taxonomy")) cat("Fossils record simulated from or assigned to an existing taxonomy \n")
    else cat("Fossil record not simulated using taxonomy: all speciation events are assumed to be symmetric \n")
  }
}

#' @export
#' @aliases fossils
summary.fossils <- function(object, max.length = 10, ...){
  print(object, max.length = max.length)
}

#' @export
#' @rdname fossils
as.fossils<-function(data, from.taxonomy = FALSE) UseMethod("as.fossils")

#' @export
as.fossils.default<-function(data, ...){
  fossils(data, ...)
}

#' @export
#' @rdname fossils
is.fossils<-function(data){
  inherits(data, "fossils")
}
