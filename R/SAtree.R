#' Tree with sampled ancestors represented as zero-length edges
#'
#' @description
#' Converts a phylo object to SAtree, without modification of tip labels.
#'
#' @param tree Phylo object.
#' @param complete Whether the tree is complete. Default TRUE. If the tree is not complete, then all fossil tips correspond to fossil samples, otherwise only sampled ancestors are considered samples.
#'
#' @export
SAtree = function(tree, complete = TRUE) {
  if(! "phylo" %in% class(tree)) stop("SAtree must be a valid phylo object")
  tree$complete = complete
  attr(tree, "class") <- c("SAtree", class(tree))
  tree
}

#' Transforms a tree and fossils dataframe to a combined SA tree.
#'
#' Sampled ancestors are represented as tips on zero-length edges to maintain compatibility with the ape format.
#' Tip labels are set to "species id"_"index". The order of the indexes is given by `tip_order`: either the oldest tip of a given species
#' receives index 1 and indexes increase towards the present (default) or the reverse.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @param tip_order Order of indexes to assign to the tips, either `oldest_first` (by default, indexes increase towards the present) or `youngest_first`
#' (indexes increase towards the past).
#' @return A list of `tree`, the SA tree integrating the fossils, and `fossils`, the fossils object updated with the tip label of each sample.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 2, tree = t)
#'
#' # transform format
#' t2 = SAtree.from.fossils(t,f)
#' plot(t2$tree)
#' @export
SAtree.from.fossils = function(tree, fossils, tip_order = c("oldest_first", "youngest_first")) {
  if(!is.fossils(fossils)) stop("Argument fossils must be a valid fossils object")

  if(length(tip_order) > 1) tip_order = tip_order[1]

  if(length(fossils[,1])==0) {
    tree$tip.label = paste0(tree$tip.label, "_", 1)
    return(list(tree = SAtree(tree, TRUE), fossils = fossils))
  }

  fossils$h = (fossils$hmin + fossils$hmax)/2
  fossils = fossils[order(fossils$sp, -fossils$h),]

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
  extant_tips = tree$tip.label
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
      if(current_spec > ntips) { #avoiding duplicates with existing tip labels
        edge_label = paste0("t", current_spec)
        while(edge_label %in% extant_tips) edge_label = paste0(edge_label, "b")
      }
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
    else tree$tip.label = c(tree$tip.label, paste0(edge_label, "_", count_spec))
    fossils$tip.label[i] = tree$tip.label[length(tree$tip.label)]
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

  if(tip_order == "youngest_first") {
    new_labels = tree$tip.label
    split_tip_labels = strsplit(tree$tip.label, split = "_", fixed = T)
    sp_labels = sapply(split_tip_labels, function(t) t[1])
    idxs = sapply(split_tip_labels, function(t) t[2])

    for(sp in unique(sp_labels)) {
      sp_positions = which(sp_labels == sp)
      sp_idxes = as.numeric(idxs[sp_positions])
      sp_idxes = max(sp_idxes) + 1 - sp_idxes
      new_labels[sp_positions] = paste0(sp, "_", sp_idxes)
    }
    names(new_labels) = tree$tip.label
    fossils$tip.label = new_labels[fossils$tip.label]
    tree$tip.label = new_labels
  }

  #force reordering for nice plotting
  attr(tree,"order")=NULL
  tree = ape::reorder.phylo(tree)

  list(tree = SAtree(tree, TRUE), fossils = fossils)
}

#' Removes all unsampled lineages from a combined tree.
#' Extinct tips are only sampled if they are fossils. With default settings all extant tips are sampled.
#'
#' @param tree Combined tree with fossils.
#' @param rho Sampling probability of extant tips. Default 1, will be disregarded if sampled_tips is not null.
#' @param sampled_tips List of tip labels corresponding to sampled extant tips.
#' @return Sampled tree with fossils.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 2, tree = t)
#'
#' # transform format
#' t2 = SAtree.from.fossils(t,f)$tree
#'
#' # transform to sampled tree
#' t3 = sampled.tree.from.combined(t2)
#' plot(t3)
#' @export
sampled.tree.from.combined = function(tree, rho = 1, sampled_tips = NULL) {
  if(!("SAtree" %in% class(tree)) ){
    if("phylo" %in% class(tree)) tree = SAtree(tree)
    else stop(paste('object "',class(tree),'" is not of class "SAtree"',sep=""))
  }
  if(!tree$complete && rho == 1 && is.null(sampled_tips)) stop("Tree is already sampled")

  remove_tips = c()

  depths = ape::node.depth.edgelength(tree)
  times = max(depths) - depths

  for(i in 1:length(tree$tip.label)) {
    if(times[i] < 1e-5) { #extant tip
      if((!is.null(sampled_tips) && !tree$tip.label[i] %in% sampled_tips) || #tip not sampled from sampled_tips
         (is.null(sampled_tips) && runif(1) > rho)) { #tip not sampled from rho
        remove_tips = c(remove_tips, i)
      }
    }
    else if(tree$complete) { #extinct tip
      edge = which(tree$edge[,2]==i)
      if(tree$edge.length[edge] > 1e-5) { #not on zero-length edge = not a fossil
        remove_tips = c(remove_tips, i)
      }
    }
  }

  tree = ape::drop.tip(tree, remove_tips)
  tree$complete = FALSE
  tree
}

#' Removes all intermediate fossils from a combined tree and labels the first and last fossils of each lineage.
#'
#' First and last are based on the order used in the `\link{SAtree.from.fossils}` function, i.e. youngest first or oldest first.
#' Can be used with sampled or complete trees. If only one fossil is present for a particular species it is labelled as first.
#'
#' @param tree Combined SAtree with fossils.
#' @return Tree with pruned fossils.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 2, tree = t)
#'
#' # transform format
#' t2 = SAtree.from.fossils(t,f)$tree
#'
#' # prune the tree to keep only ranges
#' t4 = prune.SAtree.to.ranges(t2)
#'
#' # or transform to sampled tree first
#' t3 = sampled.tree.from.combined(t2)
#' t4 = prune.SAtree.to.ranges(t3)
#' plot(t4)
#' @export
prune.SAtree.to.ranges = function(tree) {
  if(!("SAtree" %in% class(tree)) ){
    if("phylo" %in% class(tree)) tree = SAtree(tree)
    else stop(paste('object "',class(tree),'" is not of class "SAtree"',sep=""))
  }

  remove_tips = c()

  split_names = cbind(sub("_[^_]*$","",tree$tip.label),sub("^.+_","",tree$tip.label))
  for(name in unique(split_names[,1])) {
    idx = which(split_names[,1] == name)
    mn = min(split_names[idx,2])
    mx = max(split_names[idx,2])
    for(id in idx) {
      if(split_names[id,2] == mn) tree$tip.label[id] = paste0(name,"_first")
      else if(mx > mn && split_names[id,2] == mx) tree$tip.label[id] = paste0(name,"_last")
      else remove_tips = c(remove_tips, id) # intermediate sample, to remove
    }
  }

  tree = ape::drop.tip(tree, remove_tips)
  tree
}
