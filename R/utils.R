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

# Function to calculate node ages of a non-ultrametric tree
n.ages <- function(tree){

  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))

  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)

  return(node.ages)
}

# find all nodes between two nodes
# @param i the younger of the two nodes
# @param j the older of the two nodes
# @return a vector including the two nodes and any nodes in-between
find.nodes.inbetween <- function(i,j,tree){
  if(i == j) return(i)
  d = fetch.descendants(j, tree, return.nodes = TRUE)
  if(!i %in% d) stop("i not a descendant of j")
  parent = ancestor(i,tree)
  nodes = c(i)
  while(parent != j){
    nodes = c(nodes, parent)
    parent = ancestor(parent,tree)
  }
  nodes = c(nodes,j)
  return(nodes)
}

# Identify parent node
ancestor <- function(node,tree){
  parent = tree$edge[,1][which(tree$edge[,2] == node)]
  return(parent)
}

# Identify tips
#
# @param taxa Node label.
# @param tree Phylo object.
# @return Boolean (true/false).
# @examples
# t = ape::rtree(6)
# is.tip(t$edge[,2][6],t)
is.tip <- function(taxa,tree){
  return (sum(tree$edge[,1] == taxa) < 1)
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
is.extant <- function(taxa, tree, tol = NULL){
  if(is.null(tol)) tol = max((min(tree$edge.length)/100), 1e-8)
  age = n.ages(tree)[taxa]
  return(abs(age) < tol)
}

# return tip labels of extinct tips, given tolerance tol
# adapted from geiger
get.extinct.tips <- function (phy, tol = NULL) {
    if (!"phylo" %in% class(phy)) {
        stop("\"phy\" is not of class \"phylo\".");
    }
    if (is.null(phy$edge.length)) {
        stop("\"phy\" does not have branch lengths.");
    }
    if (is.null(tol)) tol <- max((min(tree$edge.length)/100), 1e-8)

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
is.root <- function(node, tree) {
  return(node == root(tree))
}


# Fetch descendant lineages from specified node in a symmetric tree
#
# @param start.node Starting node
# @param tree Phylo object.
# @param return.nodes If FALSE (default) return all descendant tip labels, otherwise return all descendant nodes (including tips).
# @examples
# t = ape::rtree(6)
# fetch.descendants(7,t)
# fetch.descendants(7,t,return.nodes=TRUE)
# @return Vector of symmetric descendants
fetch.descendants = function(start.node, tree, return.nodes = F) {

  aux = function(node) {
    result = if(return.nodes) node else c()
    if(!return.nodes && is.tip(node,tree)) result = tree$tip.label[node]

    descendants = tree$edge[which(tree$edge[,1]==node),2]
    for(d in descendants) {
      result = c(result, aux(d))
    }
    result
  }

  result = c()
  descendants = tree$edge[which(tree$edge[,1]==start.node),2]
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
#' @param tol Tolerance for identifying extant tips.
#' @return Number of extinct samples.
#'
#' @export
count.fossils = function(fossils, tol = 1e-8){
  k = length(fossils$sp[which(fossils$hmax > tol)])
  return(k)
}

#' Count the total number of fossils per interval
#'
#' @param fossils Fossils object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param tol Tolerance for identifying extant tips.
#'
#' @return Vector of extinct samples corresponding to each interval. Note the last value corresponds to the number of samples > the maximum age of the oldest interval.
#'
#' @export
count.fossils.binned = function(fossils, interval.ages, tol = 1e-8){

  if(any(fossils$hmin != fossils$hmax)) stop("Function only works for fossils with exact ages i.e. hmin = hmax");

  k = rep(0, length(interval.ages))

  if(length(fossils$sp) == 0) return(k)

  for(i in 1:length(fossils$hmax)){
    if(fossils$hmax[i] > tol){
      j = assign.interval(interval.ages, fossils$hmax[i])
      k[j] = k[j] + 1
    }
  }
  return(k)
}


###########################################
# Taxonomy functions
###########################################

#' Find a species' start (i.e speciation) time from a taxonomy object
#'
#' @param species Species id (as written in \code{taxonomy$sp}).
#' @param taxonomy Taxonomy object.
#' @return Start time.
#'
#' @export
species.start = function(species, taxonomy){
  max(taxonomy$start[which(taxonomy$sp == species)])
}

#' Find a species' end (i.e extinction) time from a taxonomy object
#'
#' @param species Species id (as written in \code{taxonomy$sp}).
#' @param taxonomy Taxonomy object.
#' @return End time.
#'
#' @export
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

# function 'untangles' (or attempts to untangle) a tree with crossing branches
# adapted from phytools
untangle<-function(tree){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  obj<-attributes(tree)
  tree<-if(length(tree$tip.label)>1) ape::read.tree(text=ape::write.tree(tree)) else tree
  ii<-!names(obj)%in%names(attributes(tree))
  attributes(tree)<-c(attributes(tree),obj[ii])
  tree
}

# truncated normal distribution
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}
