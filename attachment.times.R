#' Identify asymmetric attachment ages and extinction times in an incompletely sampled tree
#'
#' The age at which a species attaches to a tree may not be equivalent to the time of origin of a species
#' if sampling is incomplete.
#' This function takes an object of class fossils and the corresponding phylo object and calculates
#' the speciation (= attachment) times taking into account incomplete sampling.
#' The function assumes all speciation events are asymmetric (budding).
#' If the fossil object does not represent asymmetric species, the function converts species to asymmetric using the corresponding tree.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#'
#' @return Dataframe containing the speciation & extinction times in an incompletely sampled tree.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(2, t)
#'
#' # add extant samples?
#' f = add.extant.occ(f, t, rho = 0.5)
#'
#' # calculate attachment times
#' attachment.times(t, f)
#'
#' @export
attachment.times<-function(tree,fossils){

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  # asymmetric sampling
  species = create.taxonomy(tree)

  if(attr(fossils,"from.taxonomy"))
    warning("fossils already assigned based on taxonomy")
  else fossils = reconcile.fossils.taxonomy(fossils, species)

  # identify the root or origin
  if("r" %in% species$mode) root = species$sp[which(species$mode == "r")][1]
  else if ("o" %in% species$mode) root = species$sp[which(species$mode == "o")][1]

  # define the root edge
  root.edge = length(tree$tip.label) + 1

  asym.d = function(edge, tree, species){
    d = fetch.descendants(edge, tree, return.edge.labels = TRUE)
    unique(species$sp[which(species$edge %in% d)])
  }

  p <- data.frame(sp = numeric(), lineage.starts = numeric(), lineage.ends = numeric(), first.appearance = numeric(), last.appearance = numeric())

  root.decs = NULL

  for(i in unique(fossils$sp)){

    # if i is the root or origin
    if(i == root)
      attaches = i
    else {
      # identify the asymmetric parent
      parent = species$parent[which(species$sp == i)][1]
      j = i

      attaches = NULL
      while(is.null(attaches)){

        if(parent == root){
          # if the root has been sampled
          if(root %in% fossils$sp){
            attaches = j
          } else {

            # fetch the root descendants
            if(is.null(root.decs))
              root.decs = asym.d(root.edge, tree, species)

            # if i is not the only other sampled descendant
            if(length(which(root.decs[!root.decs == i] %in% fossils$sp)) > 0){

              # identify other sampled decs
              decs = root.decs[which(root.decs %in% fossils$sp)]

              # identify right most sampled dec
              right = tree$edge[,2][min(which(tree$edge[,2] %in% decs))]

              # if i/j is not the right most sample in the sym tree
              if(i != right) attaches = j
              else attaches = root

            } else attaches = root

          }
        }
        # if parent is sampled
        # -> attachment identity = self (i) or nearest ancestor (j)
        else if(parent %in% fossils$sp) attaches = j
        else {
          # identify parent edge & fetch descendants
          edge = species$edge[which(species$sp == parent)][1]
          decs = asym.d(edge, tree, species)

          # if i is not the only other sampled descendant
          if ( length(which(decs[!decs==i] %in% fossils$sp)) > 0 ){

            # identify other sampled decs
            decs = decs[which(decs %in% fossils$sp)]

            # identify right most sampled dec
            right = tree$edge[,2][min(which(tree$edge[,2] %in% decs))]

            # if i/j is not the right most sample in the sym tree
            # -> attachment identity = self i/j
            if(i != right) attaches = j
            else {
              j = parent
              parent = species$parent[which(species$sp == j)][1]
            }
          } else{
            j = parent
            parent = species$parent[which(species$sp == j)][1]
          }
        }
      }
    }
    start = species$start[which(species$sp == attaches)][1]
    end = species$end[which(species$sp == i)][1]
    fa = max(fossils$hmax[which(fossils$sp == i)])
    la = min(fossils$hmin[which(fossils$sp == i)])
    p <- rbind(p, data.frame(sp = i, lineage.starts = start, lineage.end = end, first.appearance = fa, last.appearance = la))
    #eol
  }
  return(p)
  #eof
}

#' Map asymmetric fossil lineages
#'
#' Map fossils onto a tree assuming asymmetric (budding) speciation.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @return An object of class fossils.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(tree = t, 2)
#'
#' # add extant samples
#' f = add.extant.occ(f, tree = t, rho = 0.5)
#'
#' # asymmetric mapping
#' f = asymmetric.fossil.mapping(t, f)
#'
#' @export
asymmetric.fossil.mapping<-function(tree,fossils){

  if(!"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")
  if(!"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  species = create.taxonomy(tree)

  fossils = reconcile.fossils.taxonomy(fossils, species)

  return(fossils)
}
