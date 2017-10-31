#' Fossils object
#'
#' Create a fossil record object from a dataframe or a list. The input is taken to be a dataframe or list.
#' 
#' The fossil record object contains 6 fields containing for each fossil the following data:
#' \itemize{
#'  \item \code{sp} the label of the corresponding species. This label matches with tip labels in the corresponding phylogeny
#'  and with species labels in the taxonomy record.
#'  \item \code{node} the label of the sampling node in the phylogeny, i.e the node at the end of the edge on which the fossil was sampled
#'  \item \code{origin} the label of the node where the corresponding species originated in the phylogeny
#'  \item \code{mode} the speciation mode that gave rise to the corresponding species. Can be one of \code{b} (budding/asymmetric), \code{s} (symmetric/bifurcation) or \code{a} (anagenic).
#'  \item \code{hmin} the lowest bound of the time interval in which that fossil was sampled
#'  \item \code{hmax} the highest bound of the time interval in which that fossil was sampled. Is equal to \code{hmin} if sampling times are exact.
#' }
#'
#' @param data Dataframe or list of sampled fossils. See Details for the list of required fields.
#' @param speciation.mode String indicating the speciation concept fossils were generated under (e.g. asymmetric/budding, symmetric/bifurcating, mixed). Defaults to mixed.
#'
#' @export
fossils<-function(data, speciation.mode = c("asymmetric", "symmetric", "mixed")){
  if(is.list(data)) data <- as.data.frame(data)
  
  # check for required fields
  required_fields = c("origin", "sp", "node", "mode", "hmin", "hmax")
  missing = which(! required_fields %in% colnames(data))
  if(length(missing) > 0) stop(paste0("Missing required fields: ", paste(required_fields[missing], collapse = ", ")))
  
  # check for additional fields
  additional = which(! colnames(data) %in% required_fields)
  if(length(additional) > 0) {
    warning(paste0("These fields will be discarded: ", paste(colnames(data)[additional], collapse = ", ")))
    data[,additional] = NULL
  }
  
  attr(data, "class") <- c("fossils", class(data))
  
  # handling speciation.mode attribute
  if(length(speciation.mode) > 1) speciation.mode = "mixed"
  if(!speciation.mode %in% c("asymmetric", "symmetric", "mixed")) {
    warning("Unrecognized speciation mode, defaulting to mixed")
    speciation.mode = "mixed"
  }
  attr(data, "speciation") <- speciation.mode
  
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
  if(!is.null(attr(f,"speciation")))
    cat("Speciation mode:", attr(f, "speciation"))
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

# head(f)
# tail(f)
