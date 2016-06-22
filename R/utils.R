# Function to calculate node ages of a non-ultramteric tree written by Liam J. Revell
#' @importFrom phytools nodeHeights
n.ages<-function(tree) {

  node.ages<-max(nodeHeights(tree))-nodeHeights(tree)[match(1:tree$Nnode+length(tree$tip),tree$edge[,1]),1]
  names(node.ages)<-1:tree$Nnode+length(tree$tip)

  return(node.ages)

  # EOF
}

# Identify parent nodes
ancestor<-function(edge,tree){
  edge<-edge
  tree<-tree

  parent<-tree$edge[,1][which(tree$edge[,2]==edge)]

  return(parent)
  #eof
}

# Identify the root
root<-function(tree){
  tree<-tree

  root=length(tree$tip.label)+1

  return(root)
  #eof
}

# Test is root
is.root<-function(edge,tree){
  edge<-edge
  tree<-tree

  root=length(tree$tip.label)+1

  if(edge == root)
    return(TRUE)
  else
    return(FALSE)
  #eof
}



