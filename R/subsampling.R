#' Remove fossil lineages from a tree
#'
#' @param tree an object of class "Phylo"
#'
#' @return an object of class "Phylo". If fossil lineages were found in the tree
#'   these will be pruned, if not then the original tree is returned.
#' @export
#'
prune.fossil.tips <- function(tree){
  if(ape::is.ultrametric(tree)){
    warning("Tree is ultrametric and no fossil tips can be pruned")
    return(tree)
  }

  # avoids geiger::is.extinct()
  tol <- min(tree$edge.length)/100
  tipHeight <- diag(ape::vcv.phylo(tree))

  foss <- names(which(tipHeight < max(tipHeight)-tol))
  newTree <- ape::drop.tip(tree, foss)
  return(newTree)
}

#' Obtain the tips that define each node in a tree
#'
#' @param tree an object of class "Phylo"
#'
#' @return A list of vectors, with one entry for each node consisting of the tip labels
#'   that define that node.
#' @export
#'
get.all.descs <- function(tree) {
  descs <-
    function(tree, node) {
      # function to get descendents from single node
      return(fetch.descendants(edge = node, tree = tree))
      # return(tree$tip.label[phangorn::Descendants(tree, node = node, "tips")[[1]]])
    }

  nodes <- unique(tree$edge[, 1])
  all_descs <- list()
  for (i in 1:length(nodes)) {
    all_descs[[i]] <- descs(tree, nodes[i])
  }

  names(all_descs) <- nodes
  return(all_descs)
}

#' Remove fossil samples that occur in the stem
#'
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences for "tree"
#' @param tree an object of class "Phylo"
#'
#' @return an object of class "fossils", containing only the fossil samples that
#'   occur in the crown.
#' @export
#'
remove.stem.fossils <- function(fossils, tree) {
  crown <- prune.fossil.tips(tree)
  crown <- crown$tip.label

  crownNode <- ape::getMRCA(tree, crown)
  crownTax <- fetch.descendants(crownNode, tree)
  stemTax <- setdiff(fetch.descendants(min(tree$edge[,1]), tree), crownTax)
  if(length(stemTax) == 0){
    #warning("No stem-group found in user supplied tree")
    return(fossils)
  }
  stem <- setdiff(fetch.descendants(tree$edge[,1], tree, TRUE), fetch.descendants(crownNode, tree, TRUE))

  remove <- which(fossils$sp %in% stem)
  if (length(remove > 0)) {
    fossils <- fossils[-remove, ]
    row.names(fossils) <- as.character(c(1:length(fossils$sp)))
  }

  return(fossils)
}

#' Remove stem lineages from a tree
#'
#' @param tree an object of class "Phylo"
#'
#' @return an object of class "Phylo", if stem lineages were found in the tree
#'   these will be pruned; if not then the original tree is returned.
#' @export
#'
remove.stem.lineages <- function(tree){
  crown <- prune.fossil.tips(tree)
  crown <- crown$tip.label
  crownNode <- ape::getMRCA(tree, crown)
  #crownTips <-
  #  tree$tip.label[phangorn::Descendants(tree, crownNode)[[1]]]
  crownTips <- fetch.descendants(crownNode, tree)
  if (length(crownTips) == length(tree$tip.label)) {
    warning("No stem lineages found, returning original tree")
    return(tree)
  }
  remove <- setdiff(tree$tip.label, crownTips)
  tree <- ape::drop.tip(tree, remove)
  return(tree)
}

#' Place fossil samples from one tree in another tree, or find the ancestral
#' node for each fossil sample in one tree
#'
#' If "ext.tree" is not supplied, this function will find the direct ancestral
#' node for each of the supplied fossil samples. If "ext.tree" is supplied, this
#' function will find the direct ancestral node for each fossil in "ext.tree".
#' This second behaviour is used for placing fossils simulated on a complete
#' Birth-Death tree in the extant-only counterpart tree. This results in fossil
#' samples being placed in the crown clades of the tree upon which they were
#' simulated. When "ext.tree" is supplied, any fossil samples appearing before
#' the MRCA of the crown group are discarded.
#'
#' @param tree an object of class "Phylo"
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences for the "tree" argument
#' @param ext.tree an object of class "Phylo" representing the extant
#'   countepart to "tree", this can be obtained with prune.fossil.tips(tree).
#' @return a vector of node numbers corresponding to the direct ancestor of each
#'   fossil sample in "fossils".
#' @export
#'
place.fossils <- function(tree, fossils, ext.tree) {

  if (any(fossils$sp == min(tree$edge[,1]))) {
    stop("Can't handle fossil samples on the root.edge")
  }

  # if placing in the extant tree is not required, then set it to be tree
  if (missing(ext.tree)) {
    ext.tree <- tree
  } else {
    # cant place stem group fossils in an extant tree
    fossils <- remove.stem.fossils(fossils, tree)

    if(!ape::is.ultrametric(ext.tree)){
      stop("User supplied extant tree is not ultrametric")
    }
  }

  if (!ape::is.binary.phylo(tree) |
      !ape::is.binary.phylo(ext.tree)) {
    stop("Both trees must be strictly bifurcating")
  }

  # Get the nodes that are suitable to fit fossils to
  d <- get.all.descs(ext.tree)
  nodes <- c()
  for (i in 1:length(d)) {
    nodes[i] <- ape::getMRCA(tree, d[[i]])
  }

  output_nodes <- c()

  # for each fossil, go backwards in the tree until we hit one of the suitable nodes
  for (i in 1:length(fossils$sp)) {
    #a <- phangorn::Ancestors(tree, node = fossils$sp[i], type = "all")
    a <- find.edges.inbetween(j = min(tree$edge[,1]), i = fossils$sp[i], tree = tree)[-1]
    if (length(a) == 1 && a[1] < min(nodes)) {
      # this error should not be met
      stop(paste0("fossil number ", i, " does not belong to the crown group"))
    }
    hit <- min(which(a %in% nodes))
    output_nodes[i] <- a[hit]
  }
  # Now find the comparable node in the second tree
  for (i in 1:length(output_nodes)) {
    #tmp <- tree$tip.label[phangorn::Descendants(tree,
    #        output_nodes[i])[[1]]][tree$tip.label[phangorn::Descendants(tree,
    #          output_nodes[i])[[1]]] %in% ext.tree$tip.label]
    tmp <-
      tree$tip.label[which(tree$tip.label %in% fetch.descendants(edge = output_nodes[i], tree = tree))][tree$tip.label[which(tree$tip.label %in% fetch.descendants(edge = output_nodes[i], tree = tree))] %in% ext.tree$tip.label]

    output_nodes[i] <- ape::getMRCA(ext.tree, tmp)
  }
  return(output_nodes)
}


#' Obtain a uniform random sample of fossil occurrences
#'
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences.
#' @param proportion the proportion of all fossil samples to return in the
#'   subsample.
#' @return an object of class "fossils" containing the subsampled fossil
#'   occurrences.
#' @export
#'
subsample.fossils.uniform <- function(fossils, proportion) {
  if (proportion > 1 | proportion < 0) {
    stop("proportion must be between 0 and 1")
  }
  smp <-
    sample(
      x = c(1:length(fossils$sp)),
      size = length(fossils$sp) * proportion,
      replace = FALSE
    )
  return(fossils[smp,])
}

#' Obtain a subsample of fossil occurrences containing the oldest fossil sample
#' in each node of the tree
#'
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences.
#' @param tree an object of class "Phylo", representing the tree upon which the
#'   fossil occurrences were simulated.
#' @param complete logical, if TRUE the oldest sample from each clade in the
#'   complete tree is returned, if FALSE the oldest sample from each clade in
#'   the extant only counterpart tree is returned.
#' @return an object of class "fossils" containing the subsampled fossil
#'   occurrences.
#' @export
#'
subsample.fossils.oldest <- function(fossils, tree, complete = TRUE){
  if(!complete) {
    ext <- prune.fossil.tips(tree)
    fossils <- remove.stem.fossils(fossils, tree)
    ancs <- place.fossils(tree, fossils, ext)
  } else {
    ancs <- place.fossils(tree, fossils)
  }

  smp <- c()
  for (i in 1:length(unique(ancs))) {
    x <- which(ancs == unique(ancs)[i])
    smp <- c(smp, which(fossils$hmax == max(fossils$hmax[x])))
  }
  out <- fossils[smp,]
  row.names(out) <- as.character(c(1:length(out$hmax)))
  return(out
  )
}

#' Obtain a subsample of fossil occurrences containing the youngest fossil
#' sample in each node of the tree
#'
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences.
#' @param tree an object of class "Phylo", representing the tree upon which the
#'   fossil occurrences were simulated.
#' @param complete logical, if TRUE the youngest sample from each clade in the
#'   complete tree is returned, if FALSE the youngest sample from each clade in
#'   the extant only counterpart tree is returned.
#' @return an object of class "fossils" containing the subsampled fossil
#'   occurrences.
#' @export
#'
subsample.fossils.youngest <- function(fossils, tree, complete = TRUE){
  if(!complete) {
    ext <- prune.fossil.tips(tree)
    fossils <- remove.stem.fossils(fossils, tree)
    ancs <- place.fossils(tree, fossils, ext)
  } else {
    ancs <- place.fossils(tree, fossils)
  }

  smp <- c()
  for (i in 1:length(unique(ancs))) {
    x <- which(ancs == unique(ancs)[i])
    smp <- c(smp, which(fossils$hmin == min(fossils$hmin[x])))
  }
  out <- fossils[smp,]
  row.names(out) <- as.character(c(1:length(out$hmin)))
  return(out
  )
}

#' Obtain a subsample of fossil occurrences containing the oldest and youngest
#' fossil sample for each clade of the tree
#'
#' @param fossils an object of class "fossils" that corresponds to fossil
#'   occurrences.
#' @param tree an object of class "Phylo", representing the tree upon which the
#'   fossil occurrences were simulated.
#' @param complete logical, if TRUE the oldest and youngest sample from each
#'   clade in the complete tree is returned, if FALSE the oldest and youngest
#'   sample from each clade in the extant only counterpart tree is returned.
#' @return an object of class "fossils" containing the subsampled fossil
#'   occurrences.
#' @export
subsample.fossils.oldest.and.youngest <- function(fossils, tree, complete = TRUE){
  if (!complete) {
    ext <- prune.fossil.tips(tree)
    fossils <- remove.stem.fossils(fossils, tree)
    ancs <- place.fossils(tree, fossils, ext)
  } else {
    ancs <- place.fossils(tree, fossils)
  }

  smp_1 <- c()
  smp_2 <- c()
  for (i in 1:length(unique(ancs))) {
    x <- which(ancs == unique(ancs)[i])
    smp_1 <- c(smp_1, which(fossils$hmax == max(fossils$hmax[x])))
    smp_2 <- c(smp_2, which(fossils$hmin == min(fossils$hmin[x])))
  }

  smp <- unique(c(smp_1, smp_2))
  out <- fossils[smp,]
  row.names(out) <- as.character(c(1:length(out$hmin)))
  return(out)
}

# mimics the performance of phangorn::Descendents(type="children")
get.dec.nodes <- function(tree, node){
  if(node <= length(tree$tip.label)){
    stop("node must be an internal node, not a tip")
  }

  return(tree$edge[tree$edge[, 1] == node, 2])
}

# Bind a new tip into an existing tree with a given label
# the new tip will appear as the sister taxon to the chosen tip
# "Where" is the node number of a tip
bind.to.tip <- function(tree, where, label = "Foss_1"){

  if(where > length(tree$tip.label)){
    stop("'where' must be the node number of a tip only")
  }

  tip <- ape::rtree(2)
  tip$tip.label <- c("=^%", label)
  tip <- ape::drop.tip(tip, "=^%")

  len <- which(tree$edge[, 2] == where)
  len <- tree$edge.length[len]/2

  x <- ape::bind.tree(tree, tip, where = where, position = len)

  return(x)
}



