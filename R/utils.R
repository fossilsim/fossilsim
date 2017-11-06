# Function to calculate node ages of a non-ultrametric tree using the TreeSim function getx
n.ages<-function(tree){

  node.ages <- TreeSim::getx(tree, sersampling = 1)[1:(tree$Nnode+length(tree$tip))]
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))

  return(node.ages)
}

# Identify parent nodes
ancestor<-function(edge,tree){

  parent<-tree$edge[,1][which(tree$edge[,2]==edge)]

  return(parent)
  #eof
}

# Identify tips
#
# @param taxa Edge label.
# @param tree Phylo object.
# @return Boolean (true/false).
# @examples
# t<-ape::rtree(6)
# is.tip(t$edge[,2][6],t)
is.tip<-function(taxa,phylo){

  tree<-phylo

  if (length(which(tree$edge[,1]==taxa)) < 1) {
    return(1)
  }
  else {
    return(0)
  }

  # EOF
}

# Identify the root
root<-function(tree){

  root = length(tree$tip.label) + 1

  return(root)
  #eof
}

# Test is root
is.root<-function(edge,tree){

  root=length(tree$tip.label)+1

  if(edge == root)
    return(TRUE)
  else
    return(FALSE)
  #eof
}

# fetch immediate descendants
descendants<-function(edge,tree){

  if(edge %in% tree$edge[,1])
    decs<-tree$edge[,2][which(tree$edge[,1]==edge)]
  else
    decs = NULL

  return(decs)
  #eof
}

# map a vector of node numbers from one topology to another
map_nodes<-function(x,t.old,t.new)
{
  ret = x
  for(i in 1:length(ret))
  {
    if(x[i] > length(t.old$tip.label))
    {
      st = ape::extract.clade(t.old,x[i])$tip.label
      ret[i] = phytools::findMRCA(t.new,st)
    }
    else
    {
      ret[i] = which(t.new$tip.label==t.old$tip.label[x[i]])
    }
  }
  ret
}

# find which species is on branch at time according to taxonomy
find.species.in.taxonomy = function(taxonomy, branch, time = NULL) {
  possible = which(taxonomy$edge == branch)
  if(length(possible) == 1) return(taxonomy$sp[possible])
  if(is.null(time)) stop("Multiple species found on branch, please specify a time")
  
  for(x in possible) {
    if(taxonomy$start[x] > time && taxonomy$end[x] < time) return(taxonomy$species[x])
  }
  stop("No species found, check that branch and time are compatible")
}

# get species record from taxonomy, i.e discard edge attributes
species.record.from.taxonomy = function(taxonomy) {
  taxonomy$edge = NULL
  unique(taxonomy)
}
