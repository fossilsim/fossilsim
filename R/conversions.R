#' Transforms a tree and fossils dataframe to a combined format.
#' Sampled ancestors are represented as tips on zero-length edges to maintain compatibility with the ape format.
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @return A tree integrating the fossils.
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
  
  fossils = fossils[order(fossils$species, -fossils$h),]
  
  #renaming all species not in fossils
  for(i in 1:length(tree$tip.label)) {
    if(!i %in% fossils$species) {
      tree$tip.label[i] = paste0(tree$tip.label[i], "_", 1)
    }
  }
  
  depths = ape::node.depth.edgelength(tree)
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
    edge = which(tree$edge[,2] == fossils$node[i])
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
  tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)
  
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

#' Removes all unsampled lineages from a combined tree.
#' Extinct tips are only sampled if they are fossils. With default settings all extant tips are sampled.
#'
#' @param tree Combined tree with fossils.
#' @param rho Sampling probability of extant tips. Default 1, will be disregarded if sampled_tips is not null.
#' @param sampled_tips List of tip labels corresponding to sampled extant tips.
#' @return Sampled tree with fossils.
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
  
  depths = ape::node.depth.edgelength(tree)
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
  
  tree = ape::drop.tip(tree, remove_tips)
  tree
}

#' Removes all intermediate fossils from a combined tree and labels the first and last fossils of each lineage.
#' Can be used with sampled or complete trees. If only one fossil is present for a particular species it is labeled as first.
#'
#' @param tree Combined tree with fossils.
#' @return Tree with pruned fossils.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # transform format
#' t2 = combined.tree.with.fossils(t,f)
#' # prune fossils
#' t4 = prune.fossils(t2)
#' # or transform to sampled tree first
#' t3 = sampled.tree.from.combined(t2)
#' t4 = prune.fossils(t3)
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
  
  tree = ape::drop.tip(tree, remove_tips)
  tree
}

#' Transforms a tree and fossils into a sampled tree in beast-usable format and writes it in Newick format.
#' Designed to work with FBD.
#'
#' @param tree Complete tree.
#' @param fossils fossils dataframe.
#' @param rho Sampling probability of extant tips. Default 1, will be disregarded if sampled_tips is not null.
#' @param sampled_tips List of tip labels corresponding to sampled extant tips.
#' @param ... Additional parameters will be passed to ape::write.tree
#' @return Output of write.tree.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # simulate fossils
#' f<-sim.fossils.poisson(t, 2)
#' # output for BEAST
#' format.for.beast(t, f) # output on the console
#' format.for.beast(t, f, file="example.tre") # output in file
#' @export
format.for.beast = function(tree, fossils, rho = 1, sampled_tips = NULL, ...) {
  proc_tree = prune.fossils(sampled.tree.from.combined(combined.tree.with.fossils(tree,fossils), rho = rho, sampled_tips = sampled_tips))
  ape::write.tree(proc_tree, ...)
}

#' Transforms a fossilRecordSimulation object from package paleotree to a tree and taxonomy and fossils dataframes.
#'
#' @param record fossilRecordSimulation object.
#' @return A list containing the converted tree, taxonomy and fossils
#' @examples
#' # simulate record
#' record <- paleotree::simFossilRecord(p=0.1, q=0.1,r=0.1, nruns=1,nTotalTaxa=c(30,40), nExtant=0, nSamp = c(5,25))
#' # transform format
#' l_tf = paleotreeRecordToFossils(record)
#' l_tf$tree
#' l_tf$taxonomy
#' l_tf$fossils
#' @export
# NB: taxonomy times not branch specific
# NB: modes not branch specific
# NB: not tested with cryptic speciation
# NB: origin for anagenic species set to node directly above appearance
paleotreeRecordToFossils = function(record) {
  tree = paleotree::taxa2phylo(paleotree::fossilRecord2fossilTaxa(record))
  # recording node labels to keep track after changing the phylogeny
  tree$node.label = (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)
  fossils = data.frame(hmin=numeric(),hmax = numeric(), sp=numeric(),node=numeric(),origin=numeric(),mode = character(), stringsAsFactors = F)
  taxonomy = data.frame(edge=numeric(),sp=numeric(),start=numeric(),end=numeric(),mode = character(), 
                        cryptic = numeric(), cryptic.id = numeric(), parent = numeric(),stringsAsFactors = F)
  ages = n.ages(tree)
  
  #calculate age of record - paleotree allows for fully extinct trees so the youngest sample may not be at 0
  youngest = tree$tip.label[which(ages < 1e-5)[1]]
  offset = record[[youngest]]$taxa.data['ext.time']
  root_time = 0
  
  # this is a hack
  # basically paleotree outputs bifurcation nodes as 2 consecutive zero-length edges, which is fine
  # except I need the edge leading to the tip with the dead species above the other
  # so I'm reordering all bifurcation nodes before proceeding
  for(e in 1:length(tree$edge.length)) {
    if(tree$edge.length[e] == 0 && tree$edge[e,2] <= length(record)) {
      parent = which(tree$edge[,2] == tree$edge[e,1])
      if(length(parent) > 0 && tree$edge.length[parent] == 0) {
        other_edge = which(tree$edge[,1] == tree$edge[parent,1])
        other_edge = other_edge[which(other_edge != parent)]
        tip = tree$edge[e,2]
        tree$edge[e,2] = tree$edge[other_edge,2]
        tree$edge.length[e] = tree$edge.length[other_edge]
        tree$edge[other_edge,2] = tip
        tree$edge.length[other_edge] = 0
      }
    }
  }
  
  for(i in 1:length(record)) {
    if(record[[i]]$taxa.data['orig.time'] > root_time) root_time = record[[i]]$taxa.data[['orig.time']]
    
    # recording positions of species on phylogeny
    node_idx = which(tree$tip.label == names(record)[i])
    sampled_nodes = c()
    if(tree$edge.length[which(tree$edge[,2] == node_idx)] == 0) {
      node_idx = tree$edge[which(tree$edge[,2] == node_idx),1]
    }
    all_nodes = node_idx
    age = offset + ages[node_idx] + tree$edge.length[which(tree$edge[,2] == node_idx)] #age of the parent
    
    for(t in sort(record[[i]]$sampling.times)) {
      while(node_idx != length(record) + 1 && t > age) {
        node_idx = tree$edge[which(tree$edge[,2] == node_idx),1]
        age = age + tree$edge.length[which(tree$edge[,2] == node_idx)]
        all_nodes = c(all_nodes, node_idx)
      }
      sampled_nodes = c(sampled_nodes, node_idx)
    }
    
    while(node_idx != length(record) + 1 && record[[i]]$taxa.data[["orig.time"]] > age) {
      node_idx = tree$edge[which(tree$edge[,2] == node_idx),1]
      age = age + tree$edge.length[which(tree$edge[,2] == node_idx)]
      all_nodes = c(all_nodes, node_idx)
    }
    
    if(is.na(record[[i]]$taxa.data[["ancestor.id"]])) {
      anc_node = NA
      mode = 'root'
    }
    else {
      anc_node = tree$edge[which(tree$edge[,2] == node_idx),1]
      # checking for mode
      if(anc_node != length(record) +1 && tree$edge.length[which(tree$edge[,2] == anc_node)] == 0) {
        mode = 's' # bifurcation event
      } else {
        desc_anc = which(tree$edge[,1] == anc_node)
        other_desc = desc_anc[which(tree$edge[desc_anc,2] != node_idx)]
        if(tree$edge.length[other_desc] == 0) mode = 'a' # anagenic event
        else mode = 'b' # budding
      }
    }
    
    taxonomy = rbind(taxonomy, data.frame(edge = all_nodes, sp = names(record)[i],start=record[[i]]$taxa.data[['orig.time']],end=record[[i]]$taxa.data[['ext.time']],
                                          mode = mode, cryptic = !(record[[i]]$taxa.data[['taxon.id']] == record[[i]]$taxa.data[['looks.like']]), 
                                          cryptic.id = names(record)[record[[i]]$taxa.data[['looks.like']]], 
                                          parent = names(record)[record[[i]]$taxa.data[['ancestor.id']]], stringsAsFactors = F))
    if(length(record[[i]]$sampling.times)>0)
      fossils = rbind(fossils, data.frame(hmin = sort(record[[i]]$sampling.times), hmax = sort(record[[i]]$sampling.times), 
                                          sp = names(record)[i], node = sampled_nodes, origin = anc_node, mode = mode, stringsAsFactors = F))
  }
  row.names(taxonomy) = NULL
  row.names(fossils) = NULL
  
  # drop zero-edges from phylogeny produced by paleotree for anagenic/bifurcation events
  remove_tips = c()
  for(i in 1:length(tree$tip.label)) {
    e = which(tree$edge[,2] == i)
    if(tree$edge.length[e] == 0) {
      remove_tips = c(remove_tips, i)
      anc_edges = which(tree$edge[,1] == tree$edge[e,1])
      if(length(anc_edges) == 2) { # anagenic event: node will be removed as well so node labels need to be handled
        other_node = tree$edge[anc_edges[which(anc_edges!=e)],2]
        if(other_node > length(record)) 
          tree$node.label[which(tree$node.label == tree$node.label[tree$edge[e,1] - length(record)])] = tree$node.label[other_node - length(record)]
        else 
          tree$node.label[which(tree$node.label == tree$node.label[tree$edge[e,1] - length(record)])] = tree$tip.label[other_node]
      }
    }
  }
  # assign now all proper labels
  taxonomy$edge = c(tree$tip.label, tree$node.label)[taxonomy$edge]
  fossils$node = c(tree$tip.label, tree$node.label)[fossils$node]
  fossils$origin = c(tree$tip.label, tree$node.label)[fossils$origin]
  
  tree = ape::drop.tip(tree, remove_tips)
  # assign now all proper indexes
  taxonomy$edge = sapply(taxonomy$edge, function(x) {
    which(c(tree$tip.label, tree$node.label) == x)
  })
  fossils$node = sapply(fossils$node, function(x) {
    which(c(tree$tip.label, tree$node.label) == x)
  })
  fossils$origin = sapply(fossils$origin, function(x) {
    if(is.na(x)) return(NA)
    which(c(tree$tip.label, tree$node.label) == x)
  })
  
  tree$root.edge = root_time - tree$root.time
  tree$origin.time = root_time
  
  return(list(tree = tree, fossils = fossils, taxonomy = taxonomy))
}

#' Transforms a tree, fossils dataframe and taxonomy (optional) into a fossilRecordSimulation object from package paleotree.
#'
#' @param tree phylo object containing the tree
#' @param fossils fossils object
#' @param taxonomy optional taxonomy object. If NULL, all speciation is assumed symmetric
#' @param merge.cryptic whether cryptic species should be kept separate or merged, default FALSE
#' @return The converted paleotree record
#' @export
# NB: not tested with taxonomy or cryptic species
# NB: assumes taxonomy sorted from oldest to youngest on same branch
# NB: assumes taxonomy end and start are branch-dependent, some code can be removed if not
# NB: TODO assumes same species at branching if speciation mode not 's', probably not correct
# NB: TODO check ids in cryptic species
# NB: TODO assumes species ancestor is always species just above (except cryptic), check ?
# NB: TODO no mode in the fossils at the moment, fix
fossilsToPaleotreeRecord = function(tree, fossils, taxonomy = NULL, merge.cryptic = F) {
  
  .convert_one_species = function(current_node, ancestor, record) {
    # handling first species on branch
    if(is.null(taxonomy) || taxonomy$mode[which(taxonomy$edge == current_node)[1]] == 's') { #new species
      if(!is.null(taxonomy)) {
        id = which(taxonomy$edge == current_node)[1]
        current_species = taxonomy$sp[id]
        record[[current_species]] = list(taxa.data = c(length(record) +1, ancestor, taxonomy$start[id], taxonomy$end[id], (taxonomy$end[anaid] < 1e-3), length(record) +1), sampling.times = c())
        names(record[[current_species]]$taxa.data) = rec_names
        
        f = which(fossils$sp == taxonomy$sp[id])
      }
      else {
        if(current_node <= length(tree$tip.label)) current_species = tree$tip.label[current_node]
        else if(!is.null(tree$node.label)) current_species = tree$node.label[current_node - length(tree$tip.label)]
        else current_species = paste0("t",current_node)
        
        above_node = tree$edge[which(tree$edge[,2] == current_node),1]
        record[[current_species]] = list(taxa.data = c(length(record) +1, ancestor, node.ages[above_node], node.ages[current_node], (node.ages[current_node] < 1e-3), length(record) +1), sampling.times = c())
        names(record[[current_species]]$taxa.data) = rec_names
        
        f = which(fossils$node == current_node) #no taxonomy, all fossils on this branch belong to one species
      }
      
      # add all fossils
      for(fid in f) {
        record[[current_species]]$sampling.times = c(record[[current_species]]$sampling.times, (fossils$min[fid]+fossils$max[fid])/2)
      }
    }
    
    else { #same species
      id = which(taxonomy$edge == current_node)[1]
      current_species = names(record)[ancestor]
      record[[current_species]]$taxa.data[["ext.time"]] = taxonomy$end[id]
      # fossils already handled when species was started
    }
    
    # handling following species
    if(!is.null(taxonomy) && length(which(taxonomy$edge == current)) > 1) { #anagenic/cryptic species
      ancestor = which(names(record) == current_species)
      for(anaid in which(taxonomy$edge == current)[-1]) {
        
        if(taxonomy$cryptic[anaid] && merge.cryptic) { #in this case we attribute it to previous species
          record[[current_species]]$taxa.data[["ext.time"]] = taxonomy$end[anaid]
          
          #if fossils also merged, then already handled with the previous species
          if(!attr(fossils, "cryptic.merged")) {
            f = which(fossils$sp == taxonomy$sp[anaid])
            for(fid in f) {
              record[[current_species]]$sampling.times = c(record[[current_species]]$sampling.times, (fossils$min[fid]+fossils$max[fid])/2)
            }
          }
        }
        else {
          if(taxonomy$cryptic[anaid]) { # cryptic species, not merged TODO check proper ids here 
            current_species = taxonomy$cryptic.id[anaid]
            looks.like = which(names(record) == taxonomy$sp[id])
          }
          else {
            current_species = taxonomy$sp[anaid]
            looks.like = length(record) +1
          }
          
          record[[current_species]] = list(taxa.data = c(length(record) +1, ancestor, taxonomy$start[anaid], taxonomy$end[anaid], (taxonomy$end[anaid] < 1e-3), looks.like), sampling.times = c())
          names(record[[current_species]]$taxa.data) = rec_names
          
          f = which(fossils$sp == current_species)
          for(fid in f) {
            record[[current_species]]$sampling.times = c(record[[current_species]]$sampling.times, (fossils$min[fid]+fossils$max[fid])/2)
          }
        }
      }
    }
    
    desc = tree$edge[which(tree$edge[,1] == current_node),2]
    for(d in desc) {
      record = .convert_one_species(d, which(names(record) == current_species), record)
    }
    record
  }
  
  node.ages = ape::node.depth.edgelength(tree)
  node.ages = max(node.ages) - node.ages
  rec_names = c("taxon.id","ancestor.id","orig.time","ext.time", "still.alive","looks.like")
  
  current_node = length(tree$tip.label) + 1
  if(is.null(tree$root.edge)) {
    if(!is.null(tree$origin.time)) node.ages = node.ages + (tree$origin.time - max(node.ages))
    record = list()
    desc = tree$edge[which(tree$edge[,1] == current_node),2]
    for(d in desc) {
      record = .convert_one_species(d, NA, record)
    }
  } else {
    newn = length(tree$tip.label) + tree$Nnode
    tree$edge = rbind(c(newn, current_node), tree$edge)
    tree$edge.length = c(tree$root.edge, tree$edge.length)
    node.ages = c(node.ages, max(node.ages) + tree$root.edge)
    if(!is.null(tree$origin.time)) node.ages = node.ages + (tree$origin.time - max(node.ages))
    record = .convert_one_species(current_node, NA, list())
  }
  class(record) = "fossilRecordSimulation"
  record
}
