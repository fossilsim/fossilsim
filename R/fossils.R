#' Fossils object
#'
#' Create a fossil record object from a dataframe or a list. The input is taken to be a dataframe or list.
#'
#' The fossil record object contains 5 fields containing for each fossil the following data:
#' \itemize{
#'  \item \code{sp} the label of the corresponding species. This label matches with tip labels in the corresponding phylogeny
#'  and with species labels in the taxonomy record.
#'  \item \code{edge} the label of the sampling node in the phylogeny, i.e the node at the end of the edge on which the fossil was sampled
#'  \item \code{origin} the label of the node where the corresponding species originated in the phylogeny
#'  \item \code{hmin} the lowest bound of the time interval in which that fossil was sampled
#'  \item \code{hmax} the highest bound of the time interval in which that fossil was sampled. Is equal to \code{hmin} if sampling times are exact.
#' }
#'
#' @param data Dataframe or list of sampled fossils. See Details for the list of required fields.
#' @param from.taxonomy Boolean indicating whether the fossils were sampled using a taxonomy object. Defaults to FALSE.
#'
#' @export
fossils<-function(data, from.taxonomy = FALSE){
  if(is.list(data)) data <- as.data.frame(data)

  # check for required fields
  required_fields = c("origin", "sp", "edge", "hmin", "hmax")
  missing = which(! required_fields %in% colnames(data))
  if(length(missing) > 0) stop(paste0("Missing required fields: ", paste(required_fields[missing], collapse = ", ")))

  # check for additional fields
  additional = which(! colnames(data) %in% required_fields)
  if(length(additional) > 0) {
    warning(paste0("These fields will be discarded: ", paste(colnames(data)[additional], collapse = ", ")))
    data[,additional] = NULL
  }

  attr(data, "class") <- c("fossils", class(data))
  attr(data, "from.taxonomy") <- from.taxonomy

  return(data)
}

#' @export
#' @aliases fossils
print.fossils<-function(f, max.length = 10, ...){
  if(length(f$sp) > 0){
    if(length(f$sp) < max.length)
      max.length = length(f$sp)
    print(as.data.frame(f)[1:max.length,])
    if(length(f$sp) > max.length)
      cat("...\n")
  }
  cat("Fossil record with", length(f$sp), "occurrences representing", length(unique(f$sp)), "species\n")
  if(attr(f,"from.taxonomy")) cat("Fossils record simulated from a taxonomy")
  else cat("Fossils record not simulated from a taxonomy: all speciation was assumed symmetric")
}

#' @export
#' @aliases fossils
summary.fossils<-function(f, max.length = 10, ...){
  print(f, max.length = max.length)
}

#' @export
#' @rdname  fossils
as.fossils<-function(data, ...) UseMethod("as.fossils")

as.fossils.default<-function(data, ...){
  fossils(data, ...)
}

#' @export
#' @rdname  fossils
is.fossils<-function(fossils){
  if(inherits(fossils, "fossils"))
    TRUE
  else
    FALSE
}
