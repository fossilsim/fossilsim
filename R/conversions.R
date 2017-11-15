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

  fossils$h = (fossils$hmin + fossils$hmax)/2
  fossils = fossils[order(fossils$sp, -fossils$h),]

  ntips = length(tree$tip.label)
  totalnodes = ntips + tree$Nnode

  #renaming all species not in fossils
  for(i in 1:ntips) {
    if(!i %in% fossils$sp) {
      tree$tip.label[i] = paste0(tree$tip.label[i], "_", 1)
    }
  }

  depths = ape::node.depth.edgelength(tree)
  times = max(depths) - depths

  # adding root edge in case fossils appear on it
  if(!is.null(tree$root.edge)) {
    root = (ntips + length(fossils[,1]))*2
    tree$edge = rbind(tree$edge, c(root, ntips +1))
    tree$edge.length = c(tree$edge.length, tree$root.edge)
    times[root] = max(times) + tree$root.edge
  }

  current_spec = 0
  count_spec = 1
  for(i in 1:length(fossils[,1])) {
    if(fossils$sp[i] != current_spec) {
      tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)
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
    tree$tip.label = c(tree$tip.label, paste0(tree$tip.label[current_spec], "_", count_spec))
    count_spec = count_spec +1
  }
  tree$tip.label[current_spec] = paste0(tree$tip.label[current_spec], "_", count_spec)

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
    if(times[i] < 1e-5) { #extant tip
      if((!is.null(sampled_tips) && !tree$tip.label[i] %in% sampled_tips) || #tip not sampled from sampled_tips
         (is.null(sampled_tips) && runif(1) > rho)) { #tip not sampled from rho
        remove_tips = c(remove_tips, i)
      }
    }
    else { #extinct tip
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
#' beast.fbd.format(t, f) # output on the console
#' beast.fbd.format(t, f, file="example.tre") # output in file
#' @export
beast.fbd.format = function(tree, fossils, rho = 1, sampled_tips = NULL, ...) {
  proc_tree = prune.fossils(sampled.tree.from.combined(combined.tree.with.fossils(tree,fossils), rho = rho, sampled_tips = sampled_tips))
  ape::write.tree(proc_tree, ...)
}

#' Transforms a fossilRecordSimulation object from package paleotree to a tree and taxonomy and fossils dataframes.
#'
#' The returned tree is in paleotree format, with zero-length edges leading to tips at bifurcation and anagenic events.
#' Fossils and taxonomy are only specified on non-zero-length edges.
#'
#' @param record fossilRecordSimulation object.
#' @return A list containing the converted tree, taxonomy and fossils
#' @examples
#' # simulate record
#' record <- paleotree::simFossilRecord(p=0.1, q=0.1,r=0.1, nruns=1,nTotalTaxa=c(30,40), nExtant=0, nSamp = c(5,25))
#' # transform format
#' l_tf = paleotree.record.to.fossils(record)
#' l_tf$tree
#' l_tf$taxonomy
#' l_tf$fossils
#' @export
# NB: taxonomy times not branch specific
# NB: modes not branch specific
# NB: cryptic speciation: parent id is true parent, ie can be = cryptic id
paleotree.record.to.fossils = function(record) {
  # check that paleotree is installed - should be but you never know
  if (!requireNamespace("paleotree", quietly = TRUE)) {
    stop("Paleotree needed for this function to work. Please install it.", call. = FALSE)
  }

  tree = paleotree::taxa2phylo(paleotree::fossilRecord2fossilTaxa(record))
  # recording node labels to keep track after changing the phylogeny
  tree$node.label = (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)
  fossildf = fossils()
  taxonomy = data.frame(edge=numeric(),sp=numeric(),start=numeric(),end=numeric(),mode = character(),
                        cryptic = numeric(), cryptic.id = numeric(), parent = numeric(),stringsAsFactors = F)
  ages = n.ages(tree)

  #calculate age of record - paleotree allows for fully extinct trees so the youngest sample may not be at 0
  youngest = tree$tip.label[which(ages < 1e-5)[1]]
  offset = record[[youngest]]$taxa.data['ext.time']
  root_time = 0

  # this is a hack so the mode checking will work
  # basically paleotree outputs bifurcation nodes as 2 zero-length edges with the tip randomly assigned
  # here I'm reordering so the 2 edges always have the same ancestor
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
      prev = tree$edge[which(tree$edge[,2] == node_idx),1]
      if(prev %in% taxonomy$edge) break()
      lgth = tree$edge.length[which(tree$edge[,2] == prev)]
      if(length(lgth) > 0 && lgth == 0) break()

      node_idx = prev
      age = age + lgth
      all_nodes = c(all_nodes, node_idx)
    }

    if(is.na(record[[i]]$taxa.data[["ancestor.id"]])) { #root species
      origin = NA
      mode = 'root'
    } else {
      origin = tree$edge[which(tree$edge[,2] == node_idx),1]
      orig_edge = which(tree$edge[,2] == origin)
      if(length(orig_edge) == 0 || tree$edge.length[orig_edge] > 0) {
        desc_anc = which(tree$edge[,1] == origin)
        if(any(tree$edge.length[desc_anc] == 0)) mode = 'a' # anagenic event
        else mode = 'b' # budding
      } else {
        mode = 's' # bifurcation event
      }
    }

    taxonomy = rbind(taxonomy, data.frame(edge = all_nodes, sp = names(record)[i],start=record[[i]]$taxa.data[['orig.time']],end=record[[i]]$taxa.data[['ext.time']],
                                          mode = mode, cryptic = !(record[[i]]$taxa.data[['taxon.id']] == record[[i]]$taxa.data[['looks.like']]),
                                          cryptic.id = names(record)[record[[i]]$taxa.data[['looks.like']]],
                                          parent = names(record)[record[[i]]$taxa.data[['ancestor.id']]], stringsAsFactors = F))
    if(length(record[[i]]$sampling.times)>0)
      fossildf = rbind(fossildf, data.frame(hmin = sort(record[[i]]$sampling.times), hmax = sort(record[[i]]$sampling.times),
                                          sp = names(record)[i], edge = sampled_nodes, origin = origin))
  }

  row.names(taxonomy) = NULL
  row.names(fossildf) = NULL
  fossildf = as.fossils(fossildf, TRUE)
  taxonomy = as.taxonomy(taxonomy)

  tree$root.edge = root_time - tree$root.time
  tree$origin.time = root_time

  return(list(tree = tree, fossils = fossildf, taxonomy = taxonomy))
}

#' Transforms a fossils dataframe and either taxonomy or tree into a fossilRecordSimulation object from package paleotree.
#'
#' @param fossils fossils object
#' @param tree phylo object containing the tree. If provided and taxonomy = NULL, all speciation is assumed symmetric
#' @param taxonomy taxonomy object. If both tree and taxonomy are provided, only taxonomy will be used.
#' @return The converted paleotree record
#' @export
fossils.to.paleotree.record = function(fossils, tree = NULL, taxonomy = NULL) {
  if(is.null(taxonomy) && is.null(tree)) stop("Either tree or taxonomy needs to be provided")

  rec_names = c("taxon.id","ancestor.id","orig.time","ext.time", "still.alive","looks.like")

  if(!is.null(taxonomy)) {
    # then record based purely on taxonomy
    species = species.record.from.taxonomy(taxonomy)
    record = vector("list", length = length(species$sp))
    names(record) = species$sp
    for(i in 1:length(species$sp)) {
      if(is.na(species$parent[i])) anc = NA
      else anc = which(names(record) == species$parent[i])
      record[[i]] = list(taxa.data = c(i, anc, species$start[i], species$end[i], (species$end[i] < 1e-5), which(names(record) == species$cryptic.id[i])),
                         sampling.times = numeric())
      names(record[[i]]$taxa.data) = rec_names

      f = which(fossils$sp == species$sp[i])
      for(fid in f) {
        record[[i]]$sampling.times = c(record[[i]]$sampling.times, (fossils$hmin[fid]+fossils$hmax[fid])/2)
      }
      record[[i]]$sampling.times = sort(record[[i]]$sampling.times, decreasing = T)
    }
  } else {
    #auxiliary function for handling one edge
    .convert_one_species = function(current_node, ancestor, record) {
      if(current_node <= length(tree$tip.label)) current_species = as.character(tree$tip.label[current_node])
      else if(!is.null(tree$node.label)) current_species = as.character(tree$node.label[current_node - length(tree$tip.label)])
      else current_species = paste0("t",current_node)

      above_node = tree$edge[which(tree$edge[,2] == current_node),1]
      record[[current_species]] = list(taxa.data = c(length(record) +1, ancestor, node.ages[above_node], node.ages[current_node], (node.ages[current_node] < 1e-3), length(record) +1),
                                       sampling.times = numeric())
      names(record[[current_species]]$taxa.data) = rec_names

      f = which(fossils$edge == current_node) #no taxonomy, all fossils on this branch belong to one species
      for(fid in f) {
        record[[current_species]]$sampling.times = c(record[[current_species]]$sampling.times, (fossils$hmin[fid]+fossils$hmax[fid])/2)
      }

      desc = tree$edge[which(tree$edge[,1] == current_node),2]
      for(d in desc) {
        record = .convert_one_species(d, which(names(record) == current_species), record)
      }
      record
    }

    # no taxonomy => parse tree assuming symmetric speciation everywhere
    node.ages = ape::node.depth.edgelength(tree)
    node.ages = max(node.ages) - node.ages

    current_node = length(tree$tip.label) + 1
    if(is.null(tree$root.edge)) {
      # no root edge => assuming that the tree started with 2 new species created at the root time
      if(!is.null(tree$origin.time)) node.ages = node.ages + (tree$origin.time - max(node.ages))
      record = list()
      desc = tree$edge[which(tree$edge[,1] == current_node),2]
      for(d in desc) {
        record = .convert_one_species(d, NA, record)
      }
    } else {
      # root edge present => assuming that the tree started with one species which existed only on the root edge
      newn = length(tree$tip.label) + tree$Nnode + 1
      tree$edge = rbind(c(newn, current_node), tree$edge)
      tree$edge.length = c(tree$root.edge, tree$edge.length)
      node.ages = c(node.ages, max(node.ages) + tree$root.edge)
      if(!is.null(tree$origin.time)) node.ages = node.ages + (tree$origin.time - max(node.ages))
      record = .convert_one_species(current_node, NA, list())
    }
  }

  class(record) = "fossilRecordSimulation"
  record
}
