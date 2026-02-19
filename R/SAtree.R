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
#' The fossils object is updated with the corresponding tip labels for each sampled and returned. If a `taxonomy` object is provided, 
#' it will be updated with the added edges and also returned.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @param taxonomy Optional taxonomy map for the tree.
#' @param tip_order Order of indexes to assign to the tips, either `oldest_first` (by default, indexes increase towards the present) or `youngest_first` 
#' (indexes increase towards the past).
#' @return A list of `tree`, the SA tree integrating the fossils, `fossils`, the fossils object updated with the tip label of each sample, 
#' and updated `taxonomy` object if provided.
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
SAtree.from.fossils = function(tree, fossils, taxonomy = NULL, tip_order = c("oldest_first", "youngest_first")) {
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
    
    # changing to a new species
    if(fossils$sp[i] !=  current_spec) {
      if(current_spec <= ntips) tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)
      current_spec = fossils$sp[i]
      if(current_spec > ntips) { #avoiding duplicates with existing tip labels
        edge_label = paste0("t", current_spec)
        while(edge_label %in% extant_tips) edge_label = paste0(edge_label, "b")
      }
      count_spec = 1
    }
    
    #adding the new speciation node
    edge = which(tree$edge[,2] == fossils$edge[i])
    end_node = tree$edge[edge,2]
    tree$edge.length[edge] = times[tree$edge[edge,1]] - fossils$h[i]
    tree$edge = rbind(tree$edge, c(totalnodes+1, end_node))
    tree$edge.length = c(tree$edge.length,fossils$h[i] - times[end_node])
    tree$edge[edge,2] = totalnodes+1
    times[totalnodes+1] = fossils$h[i]
    totalnodes=totalnodes+1
    tree$Nnode=tree$Nnode+1
    
    if(!is.null(taxonomy)) {
      tax_row = which(taxonomy$edge == end_node & taxonomy$start > fossils$h[i] & taxonomy$end < fossils$h[i])
      taxonomy = rbind(taxonomy, taxonomy[tax_row, ])
      end_row = nrow(taxonomy)
      taxonomy$start[end_row] = taxonomy$end[tax_row] = fossils$h[i]
      taxonomy$edge[tax_row] = totalnodes
      taxonomy$mode[end_row] = "f"
    }
    
    #adding the fossil tip on a zero-length edge
    tree$edge = rbind(tree$edge, c(totalnodes, -i))
    tree$edge.length = c(tree$edge.length, 0)
    if(current_spec <= ntips) tree$tip.label = c(tree$tip.label, paste0(tree$tip.label[current_spec], "_", count_spec))
    else tree$tip.label = c(tree$tip.label, paste0(edge_label, "_", count_spec))
    fossils$tip.label[i] = tree$tip.label[length(tree$tip.label)]
    count_spec = count_spec +1
    
    if(!is.null(taxonomy)) {
      taxonomy = rbind(taxonomy, taxonomy[end_row, ])
      taxonomy$end[end_row + 1] = fossils$h[i]
      taxonomy$mode[end_row + 1] = "f"
      taxonomy$edge[end_row + 1] = -i
    }
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
    tree$edge[tree$edge==n] = n + length(fossils[,1])
    if(!is.null(taxonomy)) taxonomy$edge[taxonomy$edge==n] = n + length(fossils[,1])
  }
  for(i in 1:length(fossils[,1])) {
    tree$edge[tree$edge==-i] = ntips + i
    if(!is.null(taxonomy)) taxonomy$edge[taxonomy$edge==-i] = ntips + i
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
  
  list(tree = SAtree(tree, TRUE), fossils = fossils, taxonomy = taxonomy)
}
