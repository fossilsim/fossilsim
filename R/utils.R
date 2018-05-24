# Function to calculate node ages of a non-ultrametric tree using the TreeSim function getx
n.ages<-function(tree){

  depth = ape::node.depth.edgelength(tree)
  node.ages = max(depth) - depth
  names(node.ages) <- 1:(tree$Nnode+length(tree$tip))

  # adding possible offset if tree fully extinct
  if(!is.null(tree$root.time)) node.ages = node.ages + tree$root.time - max(node.ages)

  return(node.ages)
}

# Identify parent nodes
ancestor<-function(edge,tree){

  parent<-tree$edge[,1][which(tree$edge[,2]==edge)]

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
is.tip<-function(taxa,tree){
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
is.extant<-function(taxa,tree,tol=NULL){

  if(is.null(tol))
    tol = min((min(tree$edge.length)/100),1e-8)

  ages = n.ages(tree)

  end = ages[taxa]

  return(abs(end) < tol)
}


# Identify the root
root<-function(tree){

  root = length(tree$tip.label) + 1

  return(root)
}

# Test is root
is.root<-function(edge,tree){

  root=length(tree$tip.label)+1

  return(edge == root)
}

# fetch immediate descendants
descendants<-function(edge,tree){

  if(edge %in% tree$edge[,1])
    decs<-tree$edge[,2][which(tree$edge[,1]==edge)]
  else
    decs = NULL

  return(decs)
}

# map a vector of node numbers from one topology to another
map_nodes<-function(x, t.old, t.new)
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

# find all egdes between two edges
# @param i the younger of the two edges
# @param j the older of the two edges
# @return a vector including the two edges and any edges in-between
find.edges.inbetween<-function(i,j,tree){
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
fetch.descendants<-function(edge,tree,return.edge.labels=F){
  ancestor<-edge

  if(is.tip(edge, tree))
    return(NULL)

  # create vectors for nodes, tips & tracking descendents
  tips<-c()
  done<-c() # this vector contains descendants (nodes+tips)

  coi=ancestor # clade of interest
  process.complete=0
  # count=0 # debugging code

  if(is.tip(ancestor,tree)){
    tip.label=tree$tip[ancestor]
    tips<-c(tips,tip.label)
  }
  else{

    while(process.complete==0) {

      # fetch the two descendants
      row=which(tree$edge[,1]==ancestor)
      descendants=tree$edge[,2][row]
      d1<-descendants[1]
      d2<-descendants[2]

      if(!d1 %in% done) {
        if ((is.tip(d1,tree)) == 1) {
          done<-c(done, d1)
          tip.label=tree$tip[d1]
          tips<-c(tips,tip.label)
        }
        else {
          ancestor=d1
        }
      }
      else if (!d2 %in% done) {
        if ((is.tip(d2,tree)) == 1) {
          done<-c(done, d2)
          tip.label=tree$tip[d2]
          tips<-c(tips,tip.label)
        }
        else {
          ancestor=d2
        }
      }
      else {
        if(ancestor==coi){
          process.complete=1
        }
        else {
          done<-c(done, ancestor)
          row=which(tree$edge[,2]==ancestor)
          ancestor=tree$edge[,1][row]
        }
      }
      #	if (count==100) {
      #		process.complete=1
      #	}

      # count=count+1
    }
  }
  if(return.edge.labels) return(done)
  else return(tips)
  # eof
}

