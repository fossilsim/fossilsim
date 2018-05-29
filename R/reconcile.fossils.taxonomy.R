#' Reconcile existing fossil and taxonomy objects
#'
#' This function uses edge identifiers (\code{edge}) and fossil sampling times (\code{hmin}) to reassign fossil species identifiers (\code{sp, origin}) using an existing taxonomy object.
#' It can only be used if exact fossil sampling times are known (i.e. \code{hmin = hmax}), otherwise edges containing multiple species may be indistinguishable.
#'
#' @param fossils Fossils object.
#' @param taxonomy Taxonomy object.
#'
#' @return An object of class fossils.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils using the tree
#' rate = 2
#' f = sim.fossils.poisson(rate, tree = t)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s = sim.taxonomy(t, 0.5, 1, 0.5)
#' f = reconcile.fossils.taxonomy(f, s)
#' plot(f, t)
#'
#' @export
reconcile.fossils.taxonomy = function(fossils, taxonomy){
  
  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")
  
  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")
  
  # if fossils contain edges not in taxonomy
  if(!all(fossils$edge %in% taxonomy$edge))
    stop("incompatible fossils and taxonomy: not all fossil edges found in taxonomy")
  
  if(!identical(fossils$hmin, fossils$hmax))
    stop("exact fossil sampling times must be specified to use this function (i.e. hmin = hmax)")
  
  if(attr(fossils,"from.taxonomy"))
    warning("fossils already assigned based on taxonomy")
  
  # for each fossil identify the edge
  for(i in 1:length(fossils$edge)){
    edge = fossils$edge[i]
    # identify the edges in the corresponding taxonomy obj
    j = which(taxonomy$edge == edge)
    if(length(j) == 1){
      # reassign species
      fossils$sp[i] = taxonomy$sp[j]
      fossils$origin[i] = taxonomy$origin[j]
    } else { # {more than one species is associated with the edge }
      age = fossils$hmin[i]
      edges = taxonomy$edge[j]
      j = which(taxonomy$edge %in% edges & taxonomy$start > age & taxonomy$end < age)
      # reassign species
      fossils$sp[i] = taxonomy$sp[j]
      fossils$origin[i] = taxonomy$origin[j]
    }
  }
  fossils = as.fossils(fossils, from.taxonomy = TRUE)
  return(fossils)
}