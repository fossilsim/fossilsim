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
xreassign.ages<-function(tree,fossils,basin.age,strata,poisson=FALSE,root.edge=TRUE) {

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



#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param use.species.ages If TRUE reassigned fossil ages will respect the speciation times. Default = FALSE.
#'
#' @return An object of class fossils.

reassign.ages<-function(fossils, tree = NULL, species = NULL,
                        interval.ages = NULL, basin.age = NULL, strata = NULL, use.species.ages = FALSE){

  if(!is.null(fossils))
    stop("Specify fossils object")

  if(!is.null(fossils) && !"phylo" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(is.null(tree) && is.null(species))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  #TODO: we do not need species info if use.species.ages = FALSE
  #TODO: cross check species taxonomy and fossils object
  #TODO: check fossil ages are not already binned
  #TODO: throw an error if any fossils are older than the max

  if(is.null(species)){
    species = create.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  if(is.null(interval.ages) && (is.null(basin.age) || is.null(strata)))
    stop("Intervals need to be defined by specifying either interval.ages or basin.age and strata")
  if(!is.null(basin.age) && !is.null(strata)) {
    if(!is.null(interval.ages)) warning("Two interval definitions found, using interval.ages")
    else interval.ages <- seq(0, basin.age, length = strata + 1)
  }

  # for each fossil
  for(i in 1:length(fossils$hmin)){
    int = assign.interval(interval.ages, fossils$hmin[i])

    if(!use.species.ages){
      # assign hmin and hmax using interval ages
      fossils$hmin[i] = times[int]
      fossils$hmin[i] = times[int+1]
    } else { # { assign hmin and hmax that do not violate species ages }

    }

  }
  return(fossils)
}

