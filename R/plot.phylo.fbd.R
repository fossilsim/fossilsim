#' plot.phylo.fbd: Plot object of class phylo.fbd
#'
#' @param tree FBD tree to plot.
#' @param complete Is this a complete FBD tree?
#' @keywords fossilized birth death
#' @export
plot.phylo.fbd <- function(tree, complete=FALSE){
  sa.labels = tree$tip.label[tree$edge[which(tree$edge.length == 0),2]]

  tree.sa = ape::drop.tip(tree, sa.labels, collapse.singles=F)
  tree.sa$root.edge = tree$root.edge

  tree = tree.sa

  node.ages = n.ages(tree)

  sa.nodes = as.numeric(names(which(table(tree$edge[,1])==1)))
  fossil.tips = as.numeric(which(node.ages[1:length(tree$tip.label)]>1e-7))
  extant.tips = as.numeric(which(node.ages[1:length(tree$tip.label)]<=1e-7))

  fossil.nodes = c(sa.nodes, fossil.tips)

  node_species = asymmetric.identities(tree)
  species = unique(node_species)

  root=length(tree$tip.label)+1
  origin=tree$Nnode+length(tree$tip.label)+1

  node.ages[origin] = max(node.ages)+tree$root.edge
  tree$edge = rbind(tree$edge, c(origin,root))

  # find bi, di
  bi = sapply(species, function(x) max(node.ages[tree$edge[which(tree$edge[,2] %in% which(node_species==x)),1]]))
  di = sapply(species, function(x) min(node.ages[which(node_species==x)]))

  tip.labels = asymmetric.identities(tree)[1:length(tree$tip.label)]
  node.labels = asymmetric.identities(tree)[(length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)]

  # find the index distance between
  # each birth time and its parent
  anc_branch <- function(x) {
  	a = which(bi > x)
  	i = which(bi == x)
  	a = a[a < i]
  	if( length(a) == 0 )
  		return(0)

  	i-max(a)
  }

  anc = sapply(bi, anc_branch )

  # get sampled nodes
  sampled.nodes = c(sa.nodes, extant.tips)
  if( complete == FALSE ){
    sampled.nodes = c(sampled.nodes, fossil.tips)
  }

  sampled.ages<-function(x) {
    ret = c()
    sampled_species_nodes = which(sampled.nodes %in% which(node_species==x))
    if( length(sampled_species_nodes) > 0 )
      ret = node.ages[sampled.nodes[sampled_species_nodes]]
    ret
  }

  # find oi, yi
  oi = c()
  yi = c()
  sp = c()

  sa.age = c()
  sa.sp = c()
  for( s in 1:length(species) ) {
    x = species[s]
    ages = sampled.ages(x)
    if( length( ages ) > 0 ) {
      min.age = min(ages)
      max.age = max(ages)

      oi = c(oi, max(ages))
      yi = c(yi, min(ages))
      sp = c(sp, s)

      sa.age = c(sa.age, min.age, max.age)
      sa.sp = c(sa.sp, s, s)
    }
  }

  # find last sampled age for each species
  pdi <- di

  is.sampled<-function(node) {
    ret = node %in% sampled.nodes

    desc = tree$edge[which(tree$edge[,1] == node), 2]
    for(d in desc) {
      s = is.sampled(d)
      if(s==FALSE){
        pdi[which(species==node_species[d])] <<- node.ages[node]
      }
      ret = ret || s
    }
    ret
  }

  is.sampled(origin)

  # actually do the plotting
  # make empty plot
  plot(x=NULL,y=NULL,type="n",xlim=c(max(bi),0),ylim=c(0,length(di)),axes=FALSE,xlab="",ylab="")
  par(lend=2)

  # plot species tree
  # get sampled bifurcations
  i.s = which(bi != pdi)
  bi.s = bi[i.s]
  anc.s = anc[i.s]
  y.s = (1:length(bi))[i.s]
  # get unsampled bifurcations
  i.u = which(bi == pdi)
  bi.u = bi[i.u]
  anc.u = anc[i.u]
  y.u = (1:length(bi))[i.u]
  # plot species bifurcations
  segments(bi.s,y.s,bi.s,y.s-anc.s)
  segments(bi.u,y.u,bi.u,y.u-anc.u, lty=3)
  # plot species branches
  segments(bi,1:length(bi),pdi,1:length(pdi))
  segments(pdi,1:length(pdi),di,1:length(di),lty=3)
  # plot stratigraphic ranges
  w = 0.1
  rect(oi, sp+w, yi, sp-w,col=rgb(0,0,1,0.2))
  points(sa.age, sa.sp, cex=1, pch=18)
}