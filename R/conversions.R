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
#' t2 = SAtree.from.fossils(t,f)
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

#' Returns tree and fossil objects that you can use to plot the reconstructed tree.
#'
#' Note that for datasets containing extinct only samples (& rho = 0) the ages output are scaled so that the youngest sample = 0.
#'
#' @param tree Tree object.
#' @param fossils Fossils object.
#' @param rho Extant species sampling probability. Default = 1, will be disregarded if fossils object already contains extant samples.
#'
#' @return A list containing the tree and fossil objects.
#'
#' @examples
#' # simulate tree
#' birth = 0.1
#' death = 0.05
#' tips = 10
#' t = TreeSim::sim.bd.taxa(tips, 1, birth, death)[[1]]
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 0.3, tree = t)
#'
#' # simulate extant samples
#' f = sim.extant.samples(f, tree = t, rho = 0.5)
#'
#' # plot the complete tree
#' plot(f,t)
#'
#' # generate tree & fossil objects corresponding to the reconstructed tree
#' out = reconstructed.tree.fossils.objects(f, t)
#' f.reconst = out$fossils
#' t.reconst = out$tree
#'
#' # plot the reconstructed tree
#' plot(f.reconst, t.reconst)
#'
#' @export
# function to generate tree and corresponding fossil object for the reconstructed tree
reconstructed.tree.fossils.objects = function(fossils, tree, rho = 1){

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(!(rho >= 0 && rho <= 1))
    stop("rho must be a probability between 0 and 1")

  tol = min((min(tree$edge.length)/100),1e-8)
  samp_tips = NULL

  if(any( abs(fossils$hmax) < tol )){
    # identify extant sampled tips
    samp_tips = tree$tip.label[fossils$edge[which(abs(fossils$hmax) < tol)]]
    # eliminate them from the dataframe (required for SAtree.from.fossils)
    fossils = fossils[-c(which(abs(fossils$hmax) < tol)),]
  }

  # create SAtree object
  sa.tree = SAtree.from.fossils(tree, fossils)

  # match samp_tips and sa.tree tip labels
  if(!is.null(samp_tips)){
    for(i in 1:length(samp_tips)){
      samp_tips[i] = sa.tree$tip.label [ grep(paste0(samp_tips[i],"_"), sa.tree$tip.label)[1] ]
    }
  }

  # removes all unsampled lineages
  samp.tree = sampled.tree.from.combined(sa.tree, rho = rho, sampled_tips = samp_tips)

  # identify sampled ancestors and drop from the sampled tree
  sa = c()
  for(i in 1:length(samp.tree$tip.label)){
    blength = samp.tree$edge.length[which(samp.tree$edge[,2]==i)]
    if(abs(blength) < tol) sa = c(sa, samp.tree$tip.label[i])
  }
  no.sa.tree = ape::drop.tip(samp.tree, sa)

  # create new fossils object based on the reconstructed tree
  # & deal with sampled ancestors
  if(length(fossils$sp) > 0){
    f.new = data.frame()
    nages = n.ages(samp.tree)
    for(i in sa){
      anc = ancestor(which(samp.tree$tip.label==i), samp.tree)

      s1 = ape::extract.clade(samp.tree,anc)$tip.label
      s2 = no.sa.tree$tip.label[which(no.sa.tree$tip.label %in% s1)]

      # assign fossils to nearest descendent node in the tree
      if(length(s2) > 1){
        j = ape::getMRCA(no.sa.tree,s2)
      } else { # i is SA on terminal branch
        j = which(no.sa.tree$tip.label==s2)
      }
      h = nages[[which(samp.tree$tip.label==i)]]
      f.new = rbind(f.new, data.frame(sp = j,
                                      edge = j,
                                      hmin = h,
                                      hmax = h))
    }
    f.new = fossils(f.new)
  } else {
    f.new = fossils()
  }

  # deal with non SA samples
  f.new = FossilSim::sim.tip.samples(f.new, no.sa.tree)

  # add a root edge for any stem fossils
  if( any(f.new$edge == (length(no.sa.tree$tip.label)+1) )){
    h = f.new$hmin[which(f.new$edge == (length(no.sa.tree$tip.label)+1))]
    no.sa.tree$root.edge = max(h) - max(n.ages(no.sa.tree))
  }

  return(list(tree = no.sa.tree, fossils = f.new))

}

#' Removes all intermediate fossils from a combined tree and labels the first and last fossils of each lineage.
#' Can be used with sampled or complete trees. If only one fossil is present for a particular species it is labelled as first.
#'
#' @param tree Combined tree with fossils.
#' @return Tree with pruned fossils.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 2, tree = t)
#'
#' # transform format
#' t2 = SAtree.from.fossils(t,f)
#'
#' # prune fossils
#' t4 = prune.fossils(t2)
#'
#' # or transform to sampled tree first
#' t3 = sampled.tree.from.combined(t2)
#' t4 = prune.fossils(t3)
#' plot(t4)
#' @export
prune.fossils = function(tree) {
  if(!("SAtree" %in% class(tree)) ){
    if("phylo" %in% class(tree)) tree = SAtree(tree)
    else stop(paste('object "',class(tree),'" is not of class "SAtree"',sep=""))
  }

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
#' t = ape::rtree(6)
#'
#' # simulate fossils
#' f = sim.fossils.poisson(rate = 2, tree = t)
#'
#' # output for BEAST
#' beast.fbd.format(t, f) # output on the console
#' beast.fbd.format(t, f, file="example.tre") # output in file
#' @export
beast.fbd.format = function(tree, fossils, rho = 1, sampled_tips = NULL, ...) {
  proc_tree = prune.fossils(sampled.tree.from.combined(SAtree.from.fossils(tree,fossils), rho = rho, sampled_tips = sampled_tips))
  ape::write.tree(proc_tree, ...)
}

#' Transforms a fossilRecordSimulation object from package paleotree to a tree and taxonomy and fossils objects.
#'
#' The returned tree is in paleotree format, with zero-length edges leading to tips at bifurcation and anagenetic events.
#' Fossils and taxonomy are only specified on non-zero-length edges.
#' The label assigned to the parent of the origin or root will be zero.
#'
#' @param record fossilRecordSimulation object.
#' @param alphanumeric If TRUE function will return alphanumeric species labels (i.e. species labels contain the "t" prefix) (default). If FALSE function will return numeric only species labels.
#' @return A list containing the converted tree, taxonomy and fossils
#' @examples
#' if (requireNamespace("paleotree", quietly = TRUE)) {
#' # simulate record
#' record = paleotree::simFossilRecord(p=0.1, q=0.1,r=0.1, nruns=1, nTotalTaxa=c(30,40),
#'     nExtant=0, nSamp = c(5,25))
#'
#' # transform format
#' l_tf = paleotree.record.to.fossils(record)
#' l_tf$tree
#' l_tf$taxonomy
#' l_tf$fossils
#' }
#'
#' @export
#' @seealso \code{\link{taxonomy}}, \code{\link{fossils}}, \code{\link{fossils.to.paleotree.record}}
# NB: modes not branch specific
# NB: cryptic speciation: parent id is true parent, ie can be = cryptic id
paleotree.record.to.fossils = function(record, alphanumeric = TRUE) {
  # check that paleotree is installed - should be but you never know
  if (!requireNamespace("paleotree", quietly = TRUE)) {
    stop("Paleotree needed for this function to work. Please install it.", call. = FALSE)
  }

  tree = paleotree::taxa2phylo(paleotree::fossilRecord2fossilTaxa(record))
  # recording node labels to keep track after changing the phylogeny
  tree$node.label = (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)
  fossildf = fossils()
  taxonomy = data.frame(edge = numeric(), sp = numeric(), start = numeric(), end = numeric(), mode = character(),
                        cryptic = numeric(), cryptic.id = numeric(), parent = numeric(), stringsAsFactors = F)
  ages = n.ages(tree)

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
    age = ages[node_idx] + tree$edge.length[which(tree$edge[,2] == node_idx)] #age of the parent
    ends = record[[i]]$taxa.data[['ext.time']]
    starts = min(age, record[[i]]$taxa.data[["orig.time"]])

    for(t in sort(record[[i]]$sampling.times)) {
      while(node_idx != length(record) + 1 && t > age) {
        ends = c(ends, age)
        node_idx = tree$edge[which(tree$edge[,2] == node_idx),1]
        age = age + tree$edge.length[which(tree$edge[,2] == node_idx)]
        all_nodes = c(all_nodes, node_idx)
        starts = c(starts, min(age, record[[i]]$taxa.data[["orig.time"]]))
      }
      sampled_nodes = c(sampled_nodes, node_idx)
    }

    while(node_idx != length(record) + 1 && record[[i]]$taxa.data[["orig.time"]] > age) {
      prev = tree$edge[which(tree$edge[,2] == node_idx),1]
      if(prev %in% taxonomy$edge) break()
      lgth = tree$edge.length[which(tree$edge[,2] == prev)]
      if(length(lgth) > 0 && lgth == 0) break()

      ends = c(ends, age)
      node_idx = prev
      age = age + lgth
      all_nodes = c(all_nodes, node_idx)
      starts = c(starts, min(age, record[[i]]$taxa.data[["orig.time"]]))
    }

    if(is.na(record[[i]]$taxa.data[["ancestor.id"]])) { #root species
      parent = "t0"
      mode = 'r'
    } else {
      parent = names(record)[record[[i]]$taxa.data[['ancestor.id']]]
      origin = tree$edge[which(tree$edge[,2] == node_idx),1]
      orig_edge = which(tree$edge[,2] == origin)
      if(length(orig_edge) == 0 || tree$edge.length[orig_edge] > 0) {
        desc_anc = which(tree$edge[,1] == origin)
        if(any(tree$edge.length[desc_anc] == 0)) mode = 'a' # anagenetic event
        else mode = 'b' # budding
      } else {
        mode = 's' # bifurcation event
      }
    }

    taxonomy = rbind(taxonomy, data.frame(sp = names(record)[i], edge = all_nodes, parent = parent,
                                          start = starts, end = ends,
                                          mode = mode, cryptic = !(record[[i]]$taxa.data[['taxon.id']] == record[[i]]$taxa.data[['looks.like']]),
                                          cryptic.id = names(record)[record[[i]]$taxa.data[['looks.like']]],
                                          stringsAsFactors = F))
    if(length(record[[i]]$sampling.times)>0)
      fossildf = rbind(fossildf, data.frame(sp = names(record)[i], edge = sampled_nodes,
                                            hmin = sort(record[[i]]$sampling.times), hmax = sort(record[[i]]$sampling.times),
                                            stringsAsFactors = F))
  }

  row.names(taxonomy) = NULL
  row.names(fossildf) = NULL
  fossildf = as.fossils(fossildf, TRUE)
  taxonomy = as.taxonomy(taxonomy)

  if(!alphanumeric){
    fossildf$sp = gsub("t", "", fossildf$sp)
    taxonomy$sp = gsub("t", "", taxonomy$sp)
    taxonomy$cryptic.id = gsub("t", "", taxonomy$cryptic.id)
    taxonomy$parent = gsub("t", "", taxonomy$parent)
  }

  tree$root.edge = root_time - tree$root.time
  tree$origin.time = root_time

  #removing extant samples (if present) from fossils
  ext = which(fossildf$hmax < 1e-8)
  if(length(ext > 0)) fossildf = fossildf[-ext,]

  return(list(tree = tree, fossils = fossildf, taxonomy = taxonomy))
}

#' Transforms a fossils dataframe and either taxonomy or tree into a fossilRecordSimulation object from package paleotree.
#'
#' @param fossils fossils object
#' @param tree phylo object containing the tree. If provided and taxonomy = NULL, all speciation is assumed symmetric
#' @param taxonomy taxonomy object. If both tree and taxonomy are provided, only taxonomy will be used.
#' @return The converted paleotree record
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#' # simulate fossils using taxonomy
#' s = sim.taxonomy(t, 0.5, 1, 0.5)
#' f = sim.fossils.poisson(2, taxonomy = s)
#' # transform format
#' record = fossils.to.paleotree.record(f, taxonomy = s)
#' @seealso \code{\link{taxonomy}}, \code{\link{fossils}}, \code{\link{paleotree.record.to.fossils}}
#' @export
fossils.to.paleotree.record = function(fossils, tree = NULL, taxonomy = NULL) {
  if(is.null(taxonomy) && is.null(tree)) stop("Either tree or taxonomy needs to be provided")
  fossils = sim.extant.samples(fossils, tree = tree, taxonomy = taxonomy)

  rec_names = c("taxon.id","ancestor.id","orig.time","ext.time", "still.alive","looks.like")

  if(length(fossils$sp) > 0 & !any(grepl("t", fossils$sp))){
    fossils$sp = paste0("t", fossils$sp)
  }

  if(!is.null(taxonomy)) {
    # then record based purely on taxonomy

    # add "t" prefix if missing
    if(length(taxonomy$sp) > 0 & !any(grepl("t", taxonomy$sp))){
      taxonomy$sp = paste0("t", taxonomy$sp)
      taxonomy$parent = paste0("t", taxonomy$parent)
      taxonomy$cryptic.id = paste0("t", taxonomy$cryptic.id)
    }

    # order by species id
    taxonomy = taxonomy[order(as.numeric(gsub("t", "", taxonomy$sp))),]

    species = species.record.from.taxonomy(taxonomy)
    record = vector("list", length = length(species$sp))
    names(record) = species$sp
    for(i in 1:length(species$sp)) {
      if(is.na(species$parent[i]) | species$parent[i] == 0 | species$parent[i] == "t0") anc = NA else anc = which(names(record) == species$parent[i])
      record[[i]] = list(taxa.data = c(i, anc, species$species.start[i], species$species.end[i], (species$species.end[i] < 1e-5), which(names(record) == species$cryptic.id[i])),
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

#' Generate output in the format used by the program PyRate
#'
#' @param fossils Fossils object.
#' @param python If TRUE the function outputs the data in the python format used by PyRate (default).
#' If FALSE the function outputs a tab-delimited table used by tools associated with PyRate.
#' @param traits Vector of trait values equal to the number of unique species in the fossils dataframe.
#' The order should correspond to the order in which they appear in \code{unique(fossils$sp)}.
#' @param cutoff Exclude occurrences with age uncertainty greater than this value i.e. \code{hmax - hmin > cutoff}.
#' @param random If TRUE use a random number from within the interval U(hmin, hmax) for specimen ages,
#' otherwise use the midpoint of this interval (default). Applicable only when \code{python = TRUE} and for specimens with \code{hmin != hmax}.
#' @param min Value used to represent the minimum possible interval age of extinct specimens with \code{hmin = 0}. By default \code{min = NULL} and the function will use the sampling times in the fossils dataframe.
#' @param exclude.extant.singletons If TRUE exclude species that have extant samples only (default = TRUE).
#' @param file Output file name.
#' @param use.sp.names If TRUE use the value in fossils$sp as the complete taxon name, otherwise the function adds the prefix "taxa" (default = FALSE).
#'
#' @examples
#'
#' set.seed(123)
#'
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # assign a max age based on tree height
#' max.age = tree.max(t)
#'
#' # define a set of non-uniform length intervals
#' times = c(0, sort(runif(3, min = 0, max = max.age)), max.age)
#' rates = c(1,2,3,4)
#'
#' # simulate fossils reflect age uncertainty
#' f = sim.fossils.intervals(tree = t, interval.ages = times, rates = rates,
#'     use.exact.times = FALSE)
#'
#' # simulate extant samples
#' rho = 1
#' f = sim.extant.samples(f, t, rho = 1)
#'
#' plot(f, t)
#'
#' # generate input files for pyrate
#' fossils.to.pyrate(f)
#' fossils.to.pyrate(f, python = FALSE)
#'
#' # add trait values
#' traits = runif(length(unique(f$sp)))
#' fossils.to.pyrate(f, traits = traits)
#'
#' @export
fossils.to.pyrate = function(fossils, python = TRUE, traits = NULL, cutoff = NULL, random = FALSE, min = NULL,
                             exclude.extant.singletons = TRUE, file = "", use.sp.names = FALSE){

  if(length(fossils$sp) < 1) stop("Number of specimens must be > 0")

  if(!is.null(traits) && length(unique(fossils$sp)) != length(traits))
    stop("Number of trait values must equal the number of species")

  if(!is.null(traits)) traits = data.frame(sp = unique(fossils$sp), trait = round(traits, digits = 6))

  if(!is.null(cutoff) && !(cutoff > 0)) stop("Cutoff must be > 0")

  # tolerance for extant specimens age bounds = 0
  tol = 1e-8

  if(!is.null(min)){
    fossils$hmin[which(fossils$hmin != fossils$hmax && fossils$hmin < tol)] = min
    if(any(fossils$hmin > fossils$hmax)) stop("min value has generated intervals with hmax < hmin")
  }

  # deal with specimen age options
  if(!is.null(cutoff) && cutoff > 0){
    cat(paste0("Excluding ", length(which(fossils$hmax - fossils$hmin > cutoff)), " occurrences based on cutoff...\n"))
    fossils = fossils[which(fossils$hmax - fossils$hmin < cutoff),]
    if(length(fossils$sp) < 1) stop("Number of specimens after cutoff must be > 0")
    else if(exclude.extant.singletons && length(fossils[-which(fossils$hmin == fossils$hmax & fossils$hmin < tol),]$sp) < 1) stop("Number of specimens after cutoff must be > 0")
  }

  if(python){
    # use random ages
    if(random && any(fossils$hmax != fossils$hmin)){
      fossils$hmin = round(unlist(lapply(1:dim(fossils)[1], function(x) { runif(1, fossils$hmin[x], fossils$hmax[x]) })), digits = 6)
      fossils$hmax = fossils$hmin
    } else if (any(fossils$hmax != fossils$hmin)) { # use median ages
      fossils$hmin = round(unlist(lapply(1:dim(fossils)[1], function(x) { mean( c(fossils$hmin[x], fossils$hmax[x])) })), digits = 6)
      fossils$hmax = fossils$hmin
    } else {
      fossils$hmin = round(fossils$hmin, digits = 6)
      fossils$hmax = round(fossils$hmax, digits = 6)
    }
  }

  total = length(unique(fossils$sp))

  if(python){
    cat("#!/usr/bin/env python", "from numpy import * ", "",  sep = "\n", file = file, append = FALSE)

    cat("data_1=[", sep = "", file = file, append = TRUE)

    data_1 = c()
    trait1 = c()
    taxa_names = c()
    for(i in 1:total){
      sp = unique(fossils$sp)[i]
      occs = fossils[which(fossils$sp == sp),]
      times = c()
      values = c()
      for(j in 1:length(occs$sp)){
        times = c(times, occs$hmin[j])
      }
      if(!is.null(traits)) values = c(values, traits$trait[which(traits$sp == sp)])

      if(exclude.extant.singletons && length(times) == 1 && times[1] == 0) next

      if(length(times) > 0){
        data_1 = c(data_1, paste0("array([", paste(times, collapse = ","), "])"))
        if(use.sp.names)
          taxa_names = c(taxa_names, paste0("'", sp, "'"))
        else
          taxa_names = c(taxa_names, paste0("'taxa", sp, "'"))
        if(!is.null(traits)) trait1 = c(trait1, paste0("array([", paste(values, collapse = ","), "])"))
      }
    }

    # print fossil ages
    cat(paste(data_1, collapse = ",\n"), "]\n", sep = "\n", file = file, append = TRUE)

    cat("d=[data_1]", paste0("names=['", file ,"']"), "def get_data(i): return d[i]", "def get_out_name(i): return  names[i]",
        sep = "\n", file = file, append = TRUE)

    # print species names
    cat("taxa_names=[", paste(taxa_names, collapse = ","), "]\n", sep = "", file = file, append = TRUE)

    cat("def get_taxa_names(): return taxa_names", sep = "\n", file = file, append = TRUE)

    # print trait values
    if(!is.null(traits)){
      cat("\ntrait1=[", paste(trait1, collapse = ",\n"), "\n]\n\ntraits=[trait1]\ndef get_continuous(i): return traits[i]\n", sep = "", file = file, append = TRUE)
    }

  } else {

    if(is.null(traits))
      cat("Species\tStatus\tmin_age\tmax_age\n", sep = "", file = file, append = FALSE)
    else
      cat("Species\tStatus\tmin_age\tmax_age\ttraits\n", file = file, append = FALSE)

    for(i in 1:total){
      status = NULL
      sp = unique(fossils$sp)[i]
      occs = fossils[which(fossils$sp == sp),]
      # identify extant lineages
      if(any( (occs$hmin - occs$hmax < tol) & occs$hmin < tol))
        status = "extant"
      else status = "extinct"
      if(!is.null(traits)) tr = traits$trait[which(traits$sp == sp)]

      if(exclude.extant.singletons && length(occs$sp) == 1 && occs$hmax == 0) next

      for(j in 1:length(occs$sp)){
        # skip extant samples
        if(occs$hmin[j] == occs$hmax[j] && occs$hmin[j] < tol) next
        if(is.null(traits))
          cat(paste0("taxa",i), "\t", status, "\t", occs$hmin[j], "\t", occs$hmax[j], "\n", file = file, append = TRUE)
        else
          cat(paste0("taxa",i), "\t", status, "\t", occs$hmin[j], "\t", occs$hmax[j], "\t", tr, "\n", file = file, append = TRUE)
      }
    }
  }
}
