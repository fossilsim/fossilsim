#' Tree with sampled ancestors represented as zero-length edges
#'
#' @description
#' Converts a phylo object to SAtree
#'
#' @param tree Phylo object.
#'
#' @export
SAtree = function(tree) {
  if(! "phylo" %in% class(tree)) stop("SAtree must be a valid phylo object")
  attr(tree, "class") <- c("SAtree", class(tree))
  tree
}

#' Transforms a tree and fossils dataframe to a combined SA tree.
#' Sampled ancestors are represented as tips on zero-length edges to maintain compatibility with the ape format.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @return A tree integrating the fossils.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 2, tree = t)
#'
#' # transform format
#' t2 = SAtree.from.fossils(t,f)
#' plot(t2)
#' @export
SAtree.from.fossils = function(tree, fossils) {
  if(length(fossils[,1])==0) return(tree)

  fossils$h = (fossils$hmin + fossils$hmax)/2
  fossils = fossils[order(fossils$edge, -fossils$h),]

  ntips = length(tree$tip.label)
  totalnodes = ntips + tree$Nnode

  depths = ape::node.depth.edgelength(tree)
  times = max(depths) - depths

  # adding root edge in case fossils appear on it
  if(!is.null(tree$root.edge)) {
    root = (ntips + length(fossils[,1]))*2
    tree$edge = rbind(tree$edge, c(root, ntips +1))
    tree$edge.length = c(tree$edge.length, tree$root.edge)
    times[root] = max(times) + tree$root.edge
  }

  #renaming all species not in fossils
  for(i in 1:ntips) {
    if(!i %in% fossils$sp) {
      tree$tip.label[i] = paste0(tree$tip.label[i], "_", 1)
    }
  }

  current_spec = 0
  count_spec = 1
  for(i in 1:length(fossils[,1])) {
    if(fossils$sp[i] !=  current_spec) {
      if(current_spec <= ntips) tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)
      current_spec = fossils$sp[i]
      count_spec = 1
    }
    #adding new speciation node
    edge = which(tree$edge[,2] == fossils$edge[i])
    tree$edge.length[edge] = times[tree$edge[edge,1]]-fossils$h[i]
    tree$edge = rbind(tree$edge,c(totalnodes+1,tree$edge[edge,2]))
    tree$edge.length = c(tree$edge.length,fossils$h[i]-times[tree$edge[edge,2]])
    tree$edge[edge,2]=totalnodes+1
    times[totalnodes+1] = fossils$h[i]
    totalnodes=totalnodes+1
    tree$Nnode=tree$Nnode+1

    #adding fossil tip
    tree$edge = rbind(tree$edge,c(totalnodes,-i))
    tree$edge.length = c(tree$edge.length,0)
    if(current_spec <= ntips) tree$tip.label = c(tree$tip.label, paste0(tree$tip.label[current_spec], "_", count_spec))
    else tree$tip.label = c(tree$tip.label, paste0("t", current_spec, "_", count_spec))
    count_spec = count_spec +1
  }
  if(current_spec <= ntips) tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)

  #handling root edge again, mrca may have been modified by the inclusion of fossils
  if(!is.null(tree$root.edge)) {
    rootedge = which(tree$edge[,1] == root)
    newroot = tree$edge[rootedge,2]

    rootidx = which(tree$edge == newroot)
    tree$edge[which(tree$edge == ntips + 1)] = newroot
    tree$edge[rootidx] = ntips + 1

    tree$root.edge = tree$edge.length[rootedge]
    tree$edge = tree$edge[-rootedge,]
    tree$edge.length = tree$edge.length[-rootedge]
  }

  #renumbering all nodes to maintain ape format
  for(n in totalnodes:(ntips+1)) {
    tree$edge[which(tree$edge==n)] = n + length(fossils[,1])
  }
  for(i in 1:length(fossils[,1])) {
    tree$edge[which(tree$edge==-i)] = ntips + i
  }

  #force reordering for nice plotting
  attr(tree,"order")=NULL
  tree = ape::reorder.phylo(tree)

  tree
}
