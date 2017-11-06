#' Taxonomy object
#'
#' Create a taxonomy object related species identity to a phylo object.
#'
#' @description
#' The input is a dataframe that must include the following columns: \cr \cr
#' \code{sp} = true species identity label.
#' If all species originated via budding or bifurcation this will always correspond to the terminal-most edge label (i.e. the youngest) associated with this species.
#' This may not be case if the data set also contains anagenic species. UPDATE THIS \cr \cr
#' \code{edge} = edge label of the branch in the corresponding phylo object.
#' Note some species may be associated with multiple branches. \cr \cr
#' \code{mode} = speciation mode. "o" = origin or "r" = root, edge/species that began the process. "b" = asymmetric or budding speciation. "s" = symmetric or bifurcating speciation. "a" = anagenic speciation.
#' "NA" = no speciation event is associated with the edge label or speciation mode is unknown. \cr \cr
#' Additional information: \cr \cr
#' \code{parent} = ancestor of sp. This will correspond to the terminal-most edge label of the parent species. \cr \cr
#' \code{start} = origin time of species. \cr \cr
#' \code{end} = end time of species. \cr \cr
#' \code{cryptic} = TRUE if speciation event was cryptic. \cr \cr
#' \code{cryptic.id} = cryptic species identity. If cryptic = TRUE \code{cryptic.id} will differ from the true species identity (\code{sp}). \cr
#'
#' @param data Dataframe of species taxonomy.
#'
#' @export
taxonomy<-function(data){
  if(is.null(data$sp) || is.null(data$edge) || is.null(data$mode))
    stop("Species identity, edge labels and mode must be specified using 'sp', 'edge' and 'mode'")

  me <- data
  attr(me, "class") <- c("taxonomy", class(me))
  return(me)
}

#' @export
#' @aliases taxonomy
print.taxonomy<-function(t, max.length = 50){
  summary(t, max.length = max.length, details = FALSE)
}

#' @export
#' @aliases taxomy
summary.taxonomy<-function(t, max.length = 50, details = TRUE, ...){
  if(length(t$sp) > 0){
    if(length(t$sp) < max.length)
      max.length = length(t$sp)
    print(as.data.frame(t)[1:max.length,])
    if(length(t$sp) > max.length)
      cat("...\n")
  }
  cat("Taxonomy representing", length(unique(t$sp)), "species across", length(unique(t$edge)), "edges.\n")
  if(details){
    cat("\t", length(which(t$mode == "b")), "budding species\n\t",
    length(which(t$mode == "s")), "bifurcating species\n\t",
    length(which(t$mode == "a")), "anagenic species\n\t",
    length(which(t$mode == "o")), "origin species\n\t",
    length(which(t$mode == "NA" & t$cryptic == 1)), "cryptic speciation events\n")
  }
}

#' @export
#' @rdname taxonomy
as.taxonomy<-function(data, ...) UseMethod("as.taxonomy")

as.taxonomy.default<-function(data, ...){
  taxonomy(data, ...)
}

#' @export
#' @rdname taxonomy
is.taxonomy<-function(taxonomy){
  if(inherits(taxonomy, "taxonomy"))
    TRUE
  else
    FALSE
}

#eof
