#' Fossils object
#'
#' Create a fossil record object. The input is taken to be a dataframe or list with phylo object edge labels (sp) and sampled fossil ages (h).
#'
#' @param data Dataframe or list of sampled fossils (sp = edge labels, h = ages).
#' @param ages String indicating how fossil ages are defined relative to interval ages (e.g. continuous, interval.max)
#' @param speciation.mode String indicating the speciation concept fossils were generated under (e.g. asymmetric/budding, symmetric/bifurcating, mixed)
#'
#' @export
fossils<-function(data, ages = NULL, speciation.mode = NULL){
  if(is.null(data$sp) || is.null(data$h))
    stop("Edge labels and fossil ages must be specified using 'sp' and 'h'")
  if(is.list(data))
    data <- as.data.frame(data)
  me <- data
  attr(me, "class") <- c("fossils", class(me))
  attr(me, "ages") <- ages
  attr(me, "speciation") <- speciation.mode
  return(me)
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
  cat("Fossil record with", length(f$h), "occurrences repesenting", length(unique(f$sp)), "species\n")
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
