#' Transforms a tree and fossils dataframe to combined format.
#' Sampled ancestors are represented as tips on zero-length edges to keep compatibility with the ape format.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @return A tree integrating the fossils
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # transform format
#' t2 = combined.tree.with.fossils(t,f)
#' plot(t2)
#' @export
combined.tree.with.fossils = function(tree, fossils) {
  if(length(fossils[,1])==0) return(tree)
  
  species = asymmetric.species.identities(tree)
  fossils$species = species[fossils$sp]
  fossils = fossils[order(fossils$species, -fossils$h),]
  
  #renaming all species not in fossils
  for(i in 1:length(tree$tip.label)) {
    if(!i %in% fossils$species) {
      tree$tip.label[i] = paste0(tree$tip.label[i], "_", 1)
    }
  }
  
  depths = node.depth.edgelength(tree)
  times = max(depths) - depths
  
  current_spec = 0
  count_spec = 1
  totalnodes = length(tree$tip.label) + tree$Nnode
  ntips = length(tree$tip.label)
  for(i in 1:length(fossils[,1])) {
    if(fossils$species[i] != current_spec) {
      tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)
      current_spec = fossils$species[i]
      count_spec = 1
    }
    #adding new speciation node
    edge = which(tree$edge[,2] == fossils$sp[i])
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
    tree$tip.label = c(tree$tip.label, paste0(tree$tip.label[current_spec], "_", count_spec))
    count_spec = count_spec +1
  }
  
  #renumbering all nodes to keep ape format
  for(n in totalnodes:(ntips+1)) {
    tree$edge[which(tree$edge==n)] = n + length(fossils[,1])
  }
  for(i in 1:length(fossils[,1])) {
    tree$edge[which(tree$edge==-i)] = ntips + i
  }
  
  #force reordering for nice plotting
  attr(tree,"order") =NULL
  tree = reorder.phylo(tree)
  
  tree
}

#' Removes all unsampled lineages from a combined tree.
#' Extinct tips are only sampled if they are fossils. With default settings all extant tips are sampled.
#'
#' @param tree Combined tree with fossils.
#' @param rho Sampling probability of extant tips. Default 1, will be disregarded is sampled_tips is not null.
#' @param sampled_tips List of tip labels corresponding to sampled extant tips.
#' @return Sampled tree with fossils
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # transform format
#' t2 = combined.tree.with.fossils(t,f)
#' # transform to sampled tree
#' t3 = sampled.tree.from.combined(t2)
#' plot(t3)
#' @export
sampled.tree.from.combined = function(tree, rho = 1, sampled_tips = NULL) {
  remove_tips = c()
  
  depths = node.depth.edgelength(tree)
  times = max(depths) - depths
  
  for(i in 1:length(tree$tip.label)) {
    if(times[i] < 1e-5 && #extant tip
       ((!is.null(sampled_tips) && !tree$tip.label[i] %in% sampled_tips) || #tip not sampled from sampled_tips
        (is.null(sampled_tips) && runif(1) > rho))) { #tip not sampled from rho
      remove_tips = c(remove_tips, i)
    }
    if(times[i] > 1e-5) { #extinct tip
      edge = which(tree$edge[,2]==i)
      if(tree$edge.length[edge] > 1e-5) { #not on zero-length edge = not a fossil
        remove_tips = c(remove_tips, i)
      }
    }
  }
  
  tree = drop.tip(tree, remove_tips)
  tree
}

#' Removes all intermediate fossils from a combined tree and labels the first and last fossils of each lineage. 
#' Can be used with sampled or complete trees. If only one fossil is present for a particular species it is labeled as first.
#'
#' @param tree Combined tree with fossils.
#' @return Tree with pruned fossils
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # transform format
#' t2 = combined.tree.with.fossils(t,f)
#' # prune fossils
#' t4 = prune.fossils(t2)
#' plot(t4)
#' @export
prune.fossils = function(tree) {
  remove_tips = c()
  
  split_names = cbind(sub("_[^_]*$","",tree$tip.label),sub("^.+_","",tree$tip.label))
  for(name in unique(split_names[,1])) {
    idx = which(split_names[,1] == name)
    mx = max(split_names[idx,2])
    for(id in idx) {
      if(split_names[id,2] == 1) tree$tip.label[id] = paste0(name,"_first") # 1 corresponds to oldest sample in that lineage
      else if(mx >1 && split_names[id,2] == mx) tree$tip.label[id] = paste0(name,"_last") # earliest sample
      else remove_tips = c(remove_tips, id) # intermediate sample, to remove 
    }
  }
  
  tree = drop.tip(tree, remove_tips)
  tree
}

#' Transforms a tree and fossils into a sampled tree in beast-usable format and writes it in Newick format.
#' Designed to work with FBD.
#'
#' @param tree Complete tree
#' @param fossils fossils dataframe
#' @param rho Sampling probability of extant tips. Default 1, will be disregarded is sampled_tips is not null.
#' @param sampled_tips List of tip labels corresponding to sampled extant tips.
#' @param ... Additional parameters will be passed to ape::write.tree
#' @return Output of write.tree
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # output for BEAST
#' format.for.beast(t,f) # output on the console
#' format.for.beast(t,f, file="example.tre") # output in file
#' @export
format.for.beast = function(tree, fossils, rho = 1, sampled_tips = NULL, ...) {
  proc_tree = prune.fossils(sampled.tree.from.combined(combined.tree.with.fossils(tree,fossils), rho = rho, sampled_tips = sampled_tips))
  write.tree(proc_tree, ...)
}