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
  
  depths = node.depth.edgelength(tree)
  times = max(depths) - depths
  
  current_spec = 0
  count_spec = 1
  totalnodes = length(tree$tip.label) + tree$Nnode
  ntips = length(tree$tip.label)
  for(i in 1:length(fossils[,1])) {
    if(fossils$species[i] != current_spec) {
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