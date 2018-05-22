#' Taxonomy object
#'
#' @description
#' Create a taxonomy object relating species identity to a phylo object.
#'
#' @details
#' The taxonomy object includes the following 7 fields for each edge in the corresponding phylo object:
#' \itemize{
#' \item{\code{sp} true species identity label.
#' If all species originated via budding or bifurcation this will always correspond to the terminal-most edge label (i.e. the youngest) associated with each species.
#' This may not be case if the data set also contains anagenic species, when multiple species may be associated with a single edge}
#' \item{\code{edge} edge label of the branch in the corresponding phylo object.
#' Note some species may be associated with multiple edges}
#' \item{\code{parent} = ancestor of species \code{sp}. Parent labels follow the same convention as species.
#' The label assigned to the parent of the origin or root will be zero}
#' \item{\code{start} = origin time of species}
#' \item{\code{end} = end time of species}
#' \item{\code{mode} = speciation mode. "o" = origin or "r" = root (the edge/species that began the process).
#' "b" = asymmetric or budding speciation. "s" = symmetric or bifurcating speciation. "a" = anagenic speciation}
#' \item{\code{origin} = edge beginning the species}
#' }
#'
#' Optional fields:
#' \itemize{
#' \item{\code{cryptic = TRUE} if speciation event was cryptic otherwise the function assumes \code{cryptic = FALSE}}
#' \item{\code{cryptic.id} = cryptic species identity. If cryptic = TRUE \code{cryptic.id} will differ from the true species identity \code{sp}}
#' }
#'
#' @param data Dataframe of species taxonomy. See Details for the list of required fields.
#'
#' @export
taxonomy<-function(data){
  if(is.null(data$sp) || is.null(data$edge) || is.null(data$mode) || is.null(data$end) || is.null(data$start) || is.null(data$end))
    stop("Species identity, edge labels, mode, start and end times must be specified using 'sp', 'edge', 'mode', 'start' and 'end'")

  if(is.null(data$cryptic))
    data <- cbind(data, cryptic = 0, cryptic.id = data$sp)

  me <- data
  attr(me, "class") <- c("taxonomy", class(me))
  return(me)
}

#' @export
#' @rdname summary.taxonomy
print.taxonomy<-function(x, max.length = 50, round.x = 12, ...){
  summary(x, max.length = max.length, details = FALSE)
}

#' Display taxonomy object
#'
#' @param x Taxonomy object.
#' @param max.length Max number of rows to print out.
#' @param round.x Number of decimal places to be used for species and edge ages.
#' @param details If TRUE include summary statistics.
#'
#' @export
summary.taxonomy<-function(object, max.length = 50, round.x = 12, details = TRUE, ...){

  x = data.frame(lapply(object, function(y) if(is.numeric(y)) round(y, round.x) else y))

  if(length(x$sp) > 0){
    if(length(x$sp) < max.length)
      max.length = length(x$sp)
    print(as.data.frame(x)[1:max.length,])
    if(length(x$sp) > max.length)
      cat("...\n")
  }
  cat("Taxonomy representing", length(unique(x$sp)), "species across", length(unique(x$edge)), "edges.\n")
  if(details){
    cat("\t", length(unique(x$sp[which(x$mode == "b")])), "budding species\n\t",
    length(unique(x$sp[which(x$mode == "s")])), "bifurcating species\n\t",
    length(unique(x$sp[which(x$mode == "a")])), "anagenic species\n\t",
    length(unique(x$sp[which(x$mode == "o")])), "origin species\n\t",
    length(unique(x$sp[which(x$mode == "r")])), "root species\n\t",
    length(unique(x$sp[which(x$cryptic == 1)])), "cryptic speciation events\n")
  }
}

#' @export
#' @rdname taxonomy
as.taxonomy<-function(data) UseMethod("as.taxonomy")

as.taxonomy.default<-function(data, ...){
  taxonomy(data, ...)
}

#' @export
#' @rdname taxonomy
is.taxonomy<-function(data){
  inherits(taxonomy, "taxonomy")
}
