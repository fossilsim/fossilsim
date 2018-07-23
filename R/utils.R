###########################################
# Tree functions
###########################################

#' Find the maximum age in a phylo object (root age or origin time)
#'
#' @description
#' Function returns the the root age or the origin time (if \code{root.edge = TRUE}).
#'
#' @param tree Phylo object.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return max age
#' @examples
#' t = ape::rtree(6)
#' tree.max(t, root.edge = FALSE)
#'
#' @export
tree.max = function(tree, root.edge = TRUE){
  node.ages<-n.ages(tree)
  if(root.edge && exists("root.edge",tree) )
    ba = max(node.ages) + tree$root.edge
  else
    ba = max(node.ages)

  return(ba)
}

# Function to calculate node ages of a non-ultrametric tree using the TreeSim function getx
n.ages <- function(tree){

  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))

  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)

  return(node.ages)
}

# find all egdes between two edges
# @param i the younger of the two edges
# @param j the older of the two edges
# @return a vector including the two edges and any edges in-between
find.edges.inbetween <- function(i,j,tree){
  if(i == j) return(i)
  d = fetch.descendants(j, tree, return.edge.labels = TRUE)
  if(!i %in% d) stop("i not a descendant of j")
  parent = ancestor(i,tree)
  edges = c(i)
  while(parent != j){
    edges = c(edges, parent)
    parent = ancestor(parent,tree)
  }
  edges = c(edges,j)
  return(edges)
}

# Identify parent nodes
ancestor <- function(edge,tree){
  parent = tree$edge[,1][which(tree$edge[,2]==edge)]
  return(parent)
}

# Identify tips
#
# @param taxa Edge label.
# @param tree Phylo object.
# @return Boolean (true/false).
# @examples
# t = ape::rtree(6)
# is.tip(t$edge[,2][6],t)
is.tip <- function(taxa,tree){
  return (length(which(tree$edge[,1]==taxa)) < 1)
}

# Identify extant tips
#
# @param taxa Edge label.
# @param tree Phylo object.
# @param tol Rounding error tolerance.
# @return Boolean (true/false).
# @examples
# t = ape::rtree(6)
# is.extant(t$edge[,2][6],t)
is.extant <- function(taxa,tree,tol=NULL){

  if(is.null(tol))
    tol = min((min(tree$edge.length)/100),1e-8)

  age = n.ages(tree)[taxa]

  return(abs(age) < tol)
}

# return tip.labels of extinct tips, given tolerance tol
# adapted from geiger
is.extinct <- function (phy, tol=NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("\"phy\" is not of class \"phylo\".");
    }
    if (is.null(phy$edge.length)) {
        stop("\"phy\" does not have branch lengths.");
    }
    if (is.null(tol)) {
        tol <- min(phy$edge.length)/100;
    }
    Ntip <- length(phy$tip.label)
    phy <- ape::reorder.phylo(phy);
    xx <- numeric(Ntip + phy$Nnode);
    for (i in 1:length(phy$edge[,1])) {
        xx[phy$edge[i,2]] <- xx[phy$edge[i,1]] + phy$edge.length[i];
    }
    aa <- max(xx[1:Ntip]) - xx[1:Ntip] > tol;
    if (any(aa)) {
        return(phy$tip.label[which(aa)]);
    } else {
        return(NULL);
    }
}

# Identify the root
root <- function(tree){
  return(length(tree$tip.label) + 1)
}

# Test is root
is.root <- function(edge,tree) {
  return(edge == root(tree))
}

# map a vector of node numbers from one topology to another
map_nodes<-function(x, t.old, t.new) {
  ret = x
  for(i in 1:length(ret)) {
    if(x[i] > length(t.old$tip.label)) {
      st = ape::extract.clade(t.old,x[i])$tip.label
      ret[i] = ape::getMRCA(t.new,st)
    }
    else {
      ret[i] = which(t.new$tip.label==t.old$tip.label[x[i]])
    }
  }
  ret
}

# Fetch descendant lineages in a symmetric tree
#
# @param edge Edge label.
# @param tree Phylo object.
# @param return.edge.labels If TRUE return all descendant edge labels instead of tips.
# @examples
# t = ape::rtree(6)
# fetch.descendants(7,t)
# fetch.descendants(7,t,return.edge.labels=TRUE)
# @return
# Vector of symmetric descendants
# @export
# required by find.edges.inbetween
fetch.descendants = function(edge, tree, return.edge.labels = F) {

  aux = function(node) {
    result = if(return.edge.labels) node else c()
    if(!return.edge.labels && is.tip(node,tree))  result = tree$tip.label[node]

    descendants = tree$edge[which(tree$edge[,1]==node),2]
    for(d in descendants) {
      result = c(result, aux(d))
    }
    result
  }

  result = c()
  descendants = tree$edge[which(tree$edge[,1]==edge),2]
  for(d in descendants) {
    result = c(result, aux(d))
  }

  result
}


###########################################
# Fossils functions
###########################################

#' Count the total number of fossils
#'
#' @param fossils Fossils object.
#' @return Number of extinct samples.
#'
#' @export
count.fossils = function(fossils){
  k = length(fossils$sp[which(fossils$h > 0)])
  return(k)
}

#' Count the total number of fossils per interval
#'
#' @param fossils Fossils object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#'
#' @return Vector of extinct samples corresponding to each interval. Note the last value corresponds to the number of samples > the maximum age of the oldest interval.
#'
#' @export
count.fossils.binned = function(fossils, interval.ages){
  intervals<-interval.ages

  k = rep(0, length(intervals))

  if(length(fossils$sp) == 0)
    return(k)

  for(i in 1:length(fossils$h)){
    if(fossils$h[i] != 0){
      j = assign.interval(intervals, fossils$h[i])
      k[j] = k[j] + 1
    }
  }
  return(k)
}


###########################################
# Taxonomy functions
###########################################

# find species start time in taxonomy obj
species.start = function(species, taxonomy){
  max(taxonomy$start[which(taxonomy$sp == species)])
}

# find species end time in taxonomy obj
species.end = function(species, taxonomy){
  min(taxonomy$end[which(taxonomy$sp == species)])
}

# find edge beginning species
edge.start = function(species, taxonomy){
  taxonomy$edge[which(taxonomy$sp == species
                      & taxonomy$start == species.start(species, taxonomy))]
}

# find edge ending species
edge.end = function(species, taxonomy){
  taxonomy$edge[which(taxonomy$sp == species
                      & taxonomy$end == species.end(species, taxonomy))]
}

# find which species is on branch at time according to taxonomy
find.species.in.taxonomy = function(taxonomy, branch, time = NULL) {
  possible = which(taxonomy$edge == branch)
  if(length(possible) == 1) return(taxonomy$sp[possible])
  if(is.null(time)) stop("Multiple species found on branch, please specify a time")

  for(x in possible) {
    if(taxonomy$start[x] > time && taxonomy$end[x] < time) return(taxonomy$sp[x])
  }
  stop("No species found, check that branch and time are compatible")
}

# get species record from taxonomy, i.e discard edge attributes
species.record.from.taxonomy = function(taxonomy) {
  spec = taxonomy
  spec$edge = spec$start = spec$end = NULL
  spec = unique(spec)
  spec$species.start = sapply(spec$sp, function(x) species.start(x, taxonomy))
  spec$species.end = sapply(spec$sp, function(x) species.end(x, taxonomy))
  spec
}
