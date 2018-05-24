#' Include extant samples in the fossil object, with optional rho sampling.
#'
#' @param fossils Fossils object.
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param rho Extant species sampling probability.
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
sim.extant.samples = function(fossils, tree = NULL, species = NULL, rho = 1, tol = NULL){

  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(is.null(tree) && is.null(species))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  if(!is.null(attr(fossils, "from.taxonomy"))){
    from.taxonomy = attr(fossils, "from.taxonomy")
    if(!is.null(species) & !from.taxonomy)
      stop("species taxonomy defined but fossils not based on taxonomy")
  }

  if(!(rho >= 0 && rho <= 1))
    stop("rho must be a probability between 0 and 1")

  if(is.null(species)){
    species = sim.taxonomy(tree, beta = 1)
    from.taxonomy = FALSE
  } else from.taxonomy = TRUE

  tol = min(min(tree$edge.length)/100, 1e-8)

  for (i in unique(species$sp)){

    end = species$end[which(species$sp == i)][1]

    if(!(end > (0 - tol) & end < (0 + tol))) next

    if(runif(1) < rho){
      # identify the edge ending zero
      edge = species$edge[which(species$sp == i & species$edge.end == end)]
      origin = species$origin[which(species$sp == i)][1]

      fossils<-rbind(fossils, data.frame(sp = i, edge = edge, origin = origin, hmin = 0, hmax = 0))
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
#' @param species Taxonomy object.
#' @param rho Tip sampling probability.
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
sim.tip.samples<-function(fossils, tree, species = NULL, rho = 1){

  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  if(!is.null(attr(fossils, "from.taxonomy"))){
    from.taxonomy = attr(fossils, "from.taxonomy")
    if(!is.null(species) & !from.taxonomy)
      stop("species taxonomy defined but fossils not based on taxonomy")
  }

  if(!(rho >= 0 && rho <= 1))
    stop("rho must be a probability between 0 and 1")

  if(is.null(species)){
    species = sim.taxonomy(tree, beta = 1)
    from.taxonomy = FALSE
  } else from.taxonomy = TRUE

  for (i in unique(species$sp)){

    # identify the terminal most edge
    end = species$end[which(species$sp == i)][1]
    edge = species$edge[which(species$sp == i & species$edge.end == end)][1]

    if(is.tip(edge,tree)){

      if(runif(1) < rho){
        # identify the edge ending zero
        origin = species$origin[which(species$sp == i)][1]

        fossils<-rbind(fossils, data.frame(sp = i, edge = edge, origin = origin, hmin = end, hmax = end))
      }
    }

  }
  if(!is.fossils(fossils))
    fossils = as.fossils(fossils, from.taxonomy = from.taxonomy)
  return(fossils)
}
