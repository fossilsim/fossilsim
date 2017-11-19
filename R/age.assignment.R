#' Reassign fossil ages
#'
#' Reassign fossil ages using the median age of a stratigraphic interval
#' or a random time point drawn from that interval (poisson = TRUE).
#' If a given species is only extant for a portion of the total interval duration
#' ages are reassigned using the median age of the branch duration
#' or a random time point drawn from that duration (poisson = TRUE).
#'
#' @param tree Phylo object.
#' @param fossils Fossils object.
#' @param basin.age Maximum age of the oldest horizon.
#' @param strata Number of stratigraphic intervals.
#' @param poisson If TRUE the function returns random ages.
#' @return Fossils object with reassigned fossil ages.
#' sp = edge labels. h = ages.
#' @examples
#' # simulate tree
#' t<-ape::rtree(8)
#' # assign a max age based on tree height
#' max = basin.age(t)
#' # simulate fossils
#' strata = 5
#' sampling = 1
#' f <- sim.fossils.intervals(t, basin.age = max, strata = strata, probabilities = rep(sampling,5))
#' # reassign ages
#' # TODO reassign.ages(t,f,max,strata)
#' @export
reassign.ages<-function(tree,fossils,basin.age,strata,poisson=FALSE,root.edge=TRUE) {

  interval.length = basin.age/strata
  interval.md = interval.length/2

  node.ages<-n.ages(tree)
  ld<-data.frame(lineage=numeric(),start=numeric(),end=numeric())

  speciation.mode = attr(fossils, "speciation")

  for (i in c(tree$edge[,2])){ # internal nodes + tips

    # work out the max age of the lineage (e.g. when that lineage became extant)
    # & get ancestor
    row=which(tree$edge[,2]==i)
    ancestor=tree$edge[,1][row]

    # get the age of the ancestor
    a=which(names(node.ages)==ancestor)
    lineage.start=node.ages[[a]]

    # work out the min age of the lineage (e.g. when that lineage became extinct)
    # & get the branch length
    b=tree$edge.length[row]
    lineage.end=lineage.start-b # branch length

    ld<-rbind(ld,data.frame(lineage=i,start=lineage.start,end=lineage.end))

  }

  if(root.edge && exists("root.edge",tree) ){
    root=root(tree)

    a=which(names(node.ages)==root)
    lineage.end=node.ages[[a]]

    b=tree$root.edge
    lineage.start=lineage.end+b

    ld<-rbind(ld,data.frame(lineage=root,start=lineage.start,end=lineage.end))
  }

  ld=round(ld,7)

  fossils.new=data.frame(h=numeric(),sp=numeric())

  for (f in 1:length(fossils[,1])) {

    f.id=fossils[,2][f]
    f.age=fossils[,1][f]

    l.start=ld$start[which(ld$lineage==f.id)]
    l.stop=ld$end[which(ld$lineage==f.id)]

    h.start=f.age # max
    h.stop=f.age-interval.length # min

    if ((l.start >= h.start) & (h.stop >= l.stop)){
      if(poisson)
        f.age=runif(1,min=h.stop,max=h.start)
      else
        f.age=h.start-interval.md
      fossils.new<-rbind(fossils.new,data.frame(h=f.age, sp=f.id))
    }
    else if (l.start < h.start) { # lineage originates during this interval
      if(poisson)
        f.age=runif(1,min=h.stop,max=l.start)
      else
        f.age=((l.start-h.stop)/2)+h.stop
      fossils.new<-rbind(fossils.new,data.frame(h=f.age, sp=f.id))
    }
    else if (l.stop > h.stop) { # lineage goes extinct during this interval
      if(poisson)
        f.age=runif(1,min=l.stop,max=h.start)
      else
        f.age=((h.start-l.stop)/2)+l.stop
      fossils.new<-rbind(fossils.new,data.frame(h=f.age, sp=f.id))
    }
    else {
      print("Why am I here? something went wrong")
    }
  }
  fossils.new<-fossils(fossils.new, age = "median", speciation.mode = speciation.mode)
  return(fossils.new)

  #eof
}


