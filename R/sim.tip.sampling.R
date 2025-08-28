#' Include extant samples in the fossil object, with optional rho sampling.
#'
#' @param fossils Fossils object.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param rho Extant species sampling probability. Can be a single value or a vector. Vector entries will be applied to extant tips in the order in which they appear in the taxonomy object.
#' @param tol Rounding error tolerance for tip ages.
#'
#' @return An object of class fossils containing extant tip samples equal to the age of the tips (i.e. 0.0).
#'
#' @examples
#' # simulate tree
#' lambda = 0.1
#' mu = 0.05
#' tips = 8
#' t = TreeSim::sim.bd.taxa(tips, 1, lambda, mu)[[1]]
#'
#' # simulate fossils
#' f = sim.fossils.poisson(0.5, t)
#'
#' # simulate extant samples
#' f = sim.extant.samples(f, t, rho = 0.5)
#' plot(f, t)
#'
#' @export
sim.extant.samples = function(fossils, tree = NULL, taxonomy = NULL, rho = 1, tol = NULL){

  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(is.null(tree) && is.null(taxonomy))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(taxonomy))
    warning("tree and taxonomy both defined, using taxonomy")

  if(!is.null(attr(fossils, "from.taxonomy"))){
    from.taxonomy = attr(fossils, "from.taxonomy")
    if(!is.null(taxonomy) & !from.taxonomy)
      stop("species taxonomy defined but fossils not based on taxonomy")
  }

  if(any(rho < 0 | rho > 1))
    stop("rho must be a probability between 0 and 1")

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1)
    from.taxonomy = FALSE
  } else from.taxonomy = TRUE

  if(is.null(tol))
    tol = if(is.null(tree)) 1e-8 else min(min(tree$edge.length)/100, 1e-8)

  nextant = length(which(taxonomy$end > (0 - tol) & taxonomy$end < (0 + tol)))
  if(length(rho) > 1){
    if(nextant != length(rho)){
      stop("rho must be a single value or a vector with length equal to the number of extant tips")
    }
  } else{
    rho = rep(rho, nextant)
  }

  j = 0
  for(i in unique(taxonomy$sp)){

    end = min(taxonomy$end[which(taxonomy$sp == i)])

    if(!(end > (0 - tol) & end < (0 + tol))) next
    j = j + 1

    if(runif(1) < rho[j]){
      # identify the edge ending zero
      edge = taxonomy$edge[which(abs(taxonomy$end) < tol & taxonomy$sp == i)]

      fossils<-rbind_fill(fossils, data.frame(sp = i, edge = edge, hmin = 0, hmax = 0))
    }

  }

  if(!is.fossils(fossils))
    fossils = as.fossils(fossils, from.taxonomy = from.taxonomy)
  return(fossils)
  #eof
}

#' Include extant and extinct tip samples in the fossil object, with optional rho sampling.
#'
#' @param fossils Fossils object.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param rho Tip sampling probability. Can be a single value or a vector. Vector entries will be applied to extant tips in the order in which they appear in the taxonomy object.
#'
#' @return An object of class fossils containing extant or extinct tip samples equal to the age of the tips.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(2, t)
#'
#' # simulate tip samples
#' f = sim.tip.samples(f, t, rho = 0.5)
#' plot(f, t)
#'
#' @export
# The tree is required to identify which edges are terminal.
sim.tip.samples<-function(fossils, tree, taxonomy = NULL, rho = 1){

  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(taxonomy))
    warning("tree and taxonomy both defined, using taxonomy")

  if(!is.null(attr(fossils, "from.taxonomy"))){
    from.taxonomy = attr(fossils, "from.taxonomy")
    if(!is.null(taxonomy) & !from.taxonomy)
      stop("species taxonomy defined but fossils not based on taxonomy")
  }

  if(any(rho < 0 | rho > 1))
    stop("rho must be a probability between 0 and 1")

  ntips = length(tree$tip.label)
  if(length(rho) > 1){
    if(ntips != length(rho)){
      stop("rho must be a single value or a vector with length equal to the number of tips")
    }
  } else{
    rho = rep(rho, ntips)
  }

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1)
    from.taxonomy = FALSE
  } else from.taxonomy = TRUE

  j = 0
  for (i in unique(taxonomy$sp)){

    # identify the terminal most edge
    end = min(taxonomy$end[which(taxonomy$sp == i)])
    edge = taxonomy$edge[which(taxonomy$sp == i & taxonomy$end == end)]

    if(is.tip(edge,tree)){

      j = j + 1
      if(runif(1) < rho[j]){
        fossils<-rbind_fill(fossils, data.frame(sp = i, edge = edge, hmin = end, hmax = end))
      }
    }

  }
  if(!is.fossils(fossils))
    fossils = as.fossils(fossils, from.taxonomy = from.taxonomy)
  return(fossils)
}
