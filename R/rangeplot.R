#' rangeplot: Make a stratigaphic range plot from an object of class phylo.fbd
#'
#' @param x phylo.fbd object to plot.
#' @param complete Plot unsampled species?
#' @export
rangeplot <- function(x, complete=FALSE){
  if(!("phylo.fbd" %in% class(x)) ){
    stop(paste('object "',class(x),'" is not of class "phylo.fbd"'))
  }

  sa.labels = x$tip.label[x$edge[which(x$edge.length == 0),2]]

  # collapse sampled ancestor tips into 2-degree nodes
  tree.sa = ape::drop.tip(x, sa.labels, collapse.singles=F)

  node.ages = n.ages(tree.sa)
  origin.age = x$root.edge + max(node.ages)

  node.species = asymmetric.identities(tree.sa)

  tips = 1:length(tree.sa$tip.label)
  sa.nodes = as.numeric(names(which(table(tree.sa$edge[,1])==1)))
  fossil.tips = as.numeric(which(node.ages[tips]>1e-7))

  # if we are plotting only sampled species
  # prune unsampled species tips
  if(x$complete && complete == FALSE){
    unsampled.species = node.species[fossil.tips][which(!(node.species[fossil.tips] %in% node.species[sa.nodes]))]
    unsampled.labels = tree.sa$tip.label[which(node.species[tips] %in% unsampled.species)]

    tree.sa = ape::drop.tip(x, unsampled.labels)
    tree.sa = ape::drop.tip(tree.sa, sa.labels, collapse.singles=F)

    node.ages = n.ages(tree.sa)
    node.species = asymmetric.identities(tree.sa)

    tips = 1:length(tree.sa$tip.label)
    sa.nodes = as.numeric(names(which(table(tree.sa$edge[,1])==1)))
    fossil.tips = as.numeric(which(node.ages[tips]>1e-7))
  }
  tree = tree.sa
  tree$root.edge = origin.age - max(node.ages)

  species = unique(node.species)
  num.species = length(species)
  extant.tips = as.numeric(which(node.ages[tips]<=1e-7))

  root=length(tips)+1
  origin=tree$Nnode+length(tips)+1

  node.ages[origin] = origin.age
  tree$edge = rbind(tree$edge, c(origin,root))

  # find bi, di
  bi = sapply(species, function(x) max(node.ages[tree$edge[which(tree$edge[,2] %in% which(node.species==x)),1]]))
  di = sapply(species, function(x) min(node.ages[which(node.species==x)]))

  # get sampled nodes
  sampled.nodes = c(sa.nodes, extant.tips)
  if( tree$complete == FALSE ) sampled.nodes = c(sampled.nodes, fossil.tips)

  sampled.ages<-function(x) {
    ret = c()
    sampled_species_nodes = which(sampled.nodes %in% which(node.species==x))
    if( length(sampled_species_nodes) > 0 )
      ret = node.ages[sampled.nodes[sampled_species_nodes]]
    ret
  }

  # find oi, yi
  oi = c()
  yi = c()
  ra.sp = c()

  sa.age = c()
  sa.sp = c()
  for( s in 1:num.species ) {
    x = species[s]
    ages = sampled.ages(x)
    if( length( ages ) > 0 ) {
      min.age = min(ages)
      max.age = max(ages)

      oi = c(oi, max(ages))
      yi = c(yi, min(ages))
      ra.sp = c(ra.sp, s)

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
        pdi[which(species==node.species[d])] <<- node.ages[node]
      }
      ret = ret || s
    }
    ret
  }

  is.sampled(origin)

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
  
  y = 1:num.species

  # get species coordinates 
  # lineages with sampled descendants
  i.s = which(bi != pdi)
  bi.s = bi[i.s]
  di.s = di[i.s]
  pdi.s = pdi[i.s]
  y.s = y[i.s]
  anc.s = anc[i.s]
  # lineages with no sampled descendants
  i.u = which(bi == pdi)
  bi.u = bi[i.u]
  di.u = di[i.u]
  y.u = y[i.u]
  anc.u = anc[i.u]

  # if we're only plotting lineages with sampled descendants
  # condense the coordinates
  if(complete == FALSE || length(i.u) == 0){
    bi = bi.s
    di = di.s
    anc.s = sapply(bi.s, anc_branch )

    y.s.new = order(y.s)

    y.map = y
    y.map[y.s] = y.s.new

    ra.sp = y.map[ra.sp]
    sa.sp = y.map[sa.sp]

    y.s=y.s.new
  }

  # actually do the plotting
  # make empty plot
  plot(x=NULL,y=NULL,type="n",xlim=c(max(bi),0),ylim=c(0,length(di)),axes=FALSE,xlab="",ylab="")
  par(lend=2)

  if( complete ){
    # plot bifurcations with no sampled descendants
    segments(bi.u,y.u,bi.u,y.u-anc.u, lty=3)
    # plot lineages with no sampled descendants
    segments(bi.u,y.u,di.u,y.u,lty=3)
  }
   # plot bifurcations with sampled descendants
  segments(bi.s,y.s,bi.s,y.s-anc.s)
  # plot lineages with sampled descendants
  segments(bi.s,y.s,pdi.s,y.s)
  segments(pdi.s,y.s,di.s,y.s,lty=3)

  # plot sampled ranges
  w = 0.1
  rect(oi, ra.sp+w, yi, ra.sp-w,col=rgb(0,0,1,0.2))
  # plot sampled points
  points(sa.age, sa.sp, cex=1, pch=18)
}