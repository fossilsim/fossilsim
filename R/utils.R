# Function to calculate node ages of a non-ultramteric tree using the TreeSim function getx
#' @importFrom TreeSim getx
n.ages<-function(tree){

  node.ages <- getx(tree, sersampling = 1)[1:tree$Nnode+length(tree$tip)]
  names(node.ages) <- 1:tree$Nnode+length(tree$tip)

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

  root=length(tree$tip.label)+1

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

# fetch immediate decendants
descendants<-function(edge,tree){

  if(edge %in% tree$edge[,1])
    decs<-tree$edge[,2][which(tree$edge[,1]==edge)]
  else
    decs = NULL

  return(decs)
  #eof
}
