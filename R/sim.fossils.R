#' Simulate fossils under a Poisson sampling model
#'
#' @param tree Phylo object.
#' @param phi Sampling rate.
#' @return dataframe of sampled fossils.
#' sp = edge labels. h = ages.
#' @examples
#' t<-rtree(4)
#' sim.fossils.poisson(t,1)
#' @keywords uniform preseravtion
#' @export
sim.fossils.poisson<-function(tree,phi,root.edge=T){
  tree<-tree
  lambda<-phi
  root.edge<-root.edge

  node.ages<-n.ages(tree)

  fossils<-data.frame(h=numeric(),sp=numeric())

  for (i in tree$edge[,2]){ # internal nodes + tips

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

    # sample fossil numbers from the Poisson distribution
    rand=rpois(1,b*lambda)

    if(rand > 0){
      h=runif(rand,min=lineage.end,max=lineage.start)
      fossils<-rbind(fossils,data.frame(h=h,sp=i))
    }
  }

  if(root.edge && exists("root.edge",tree) ){

    root=length(tree$tip.label)+1
    a=which(names(node.ages)==root)
    lineage.end=node.ages[[a]]

    b=tree$root.edge
    lineage.start=lineage.end+b

    # sample fossil numbers from the Poisson distribution
    rand=rpois(1,b*lambda)

    if(rand > 0){
      h=runif(rand,min=lineage.end,max=lineage.start)
      fossils<-rbind(fossils,data.frame(h=h,sp=root))
    }

  }

  return(fossils) # in this data frame h=fossil age and sp=lineage
  # EOF
}

#' Simulate fossils under a uniform model of preservation
#'
#' @param tree Phylo object.
#' @param basin.age Maximum age fossils may be sampled.
#' @param strata Number of stratigraphic horizons.
#' @param sampling Probability of sampling/preservation.
#' @return dataframe of sampled fossils.
#' sp = edge labels. h = horizon labels (= max horizon age).
#' @examples
#' t<-rtree(4)
#' sim.fossils.unif(t,1,5,0.5)
#' @keywords uniform fossil preseravtion
#' @export
sim.fossils.unif<-function(tree,basin.age,strata,sampling){
  tree<-tree
  basin.age<-basin.age
  strata<-strata
  sampling<-sampling

  s1=basin.age/strata # horizon length (= max age of youngest horizon)
  horizons<-seq(s1, basin.age, length=strata)

  node.ages<-n.ages(tree)

  fossils<-data.frame(h=numeric(),sp=numeric())

  for(h in horizons){

    h.min<-h-s1
    h.max<-h

    for (i in tree$edge[,2]){ # internal nodes + tips

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

      # if the lineage is extant during this horizon
      if ( (lineage.start >= h.min) & (lineage.end <= h.max) ) {

        # 1. generate a random number from the uniform distribution
        random.number=runif(1)

        # 2. calculate the proportion of time during each horizon the lineage is extant
        # if lineage is extant for the entire duration of the horizon
        if ( (lineage.end <= h.min) & (lineage.start >= h.max) ) {
          pr=1
        }
        # if lineage goes extinct within the horizon
        else if (lineage.end >= h.min) {
          pr=h.max-lineage.end
          pr=pr/s1
        }
        # if lineage originates within the horizon
        else {
          pr=lineage.start-h.min
          pr=pr/s1
        }

        # 3. define the probabilty
        pr = pr * sampling

        # 4. if random.number < pr { record fossil as collected }
        if (random.number <= pr ) {
          fossils<-rbind(fossils,data.frame(h=h,sp=i))
        }
      }

    }
  }

  return(fossils) # in this data frame h=horizon and sp=lineage

  # EOF
}

#' Simulate water depth profile
#'
#' @param strata Number of stratigraphic horizons.
#' @return dataframe of sampled water depths.
#' @examples
#' wd<-sim.water.depth(100)
#' plot(wd, type="l")
#' @keywords non-uniform fossil preservation
#' @export
sim.water.depth<-function(strata,depth=2,cycles=2) {

  # define the x-axis values
  x=seq(0,2,length.out=strata)

  # define y-axis values
  # a - total depth excursion - amplitude
  # b - number of cycles - period
  # 1/c - defines the relative start time of each cycle - phase shift
  # y = a * sin (b * pi * (x-1/c))
  y=depth*sin(cycles*pi*(x-1/4))

  return(data.frame(x=c(1:strata),y=y))

  # EOF
}

# Part 2. simulate non-uniform occurrences
# note that this function returns an error if profile=profile (don't really know what to do about this)
# sim.non.unif<-function(tree,basin.age,strata,PA=.5,PD=.5,DT=.5,profile) {
#
#   s1=basin.age/strata # max age of youngest horizon/horizon length
#   horizons<-seq(s1, basin.age, length=strata)
#
#   # calculate node ages
#   node.ages<-n.ages(tree)
#
#   fossils<-data.frame(x=numeric(),y=numeric())
#
#   depth.counter=0
#
#   for(h in horizons){
#
#     # define the water depth
#     depth.counter=depth.counter+1
#     current.depth=profile$y[depth.counter]
#     #print(current.depth)
#
#     h.min<-h-s1
#     h.max<-h
#
#     #cat("Horizon max:",h.max,"\n") # delete
#     #cat("Horizon min:",h.min,"\n") # delete
#
#     for (i in tree$edge[,2]){ # internal nodes + tips
#
#       #print(i) # delete
#
#       # work out the max age of the lineage (e.g. when that lineage became extant)
#       # & get ancestor
#       row=which(tree$edge[,2]==i)
#       ancestor=tree$edge[,1][row]
#
#       # get the age of the ancestor
#       a=which(names(node.ages)==ancestor)
#       lineage.start=node.ages[[a]]
#
#       # work out the min age of the lineage (e.g. when that lineage became extinct)
#       # & get the branch length
#       b=tree$edge.length[row]
#       lineage.end=lineage.start-b # branch length
#
#       #cat("lineage.start:",lineage.start,"\n") # delete
#       #cat("lineage.end:",lineage.end,"\n") # delete
#
#       # if the lineage is extant during this horizon
#       if ( (lineage.start >= h.min) & (lineage.end <= h.max) ) {
#
#         # 1. generate a random number from the uniform distribution
#         random.number=runif(1)
#
#         # 2. calculate the proportion of time during each horizon the lineage is extant
#         # if lineage is extant for the entire duration of the horizon
#         if ( (lineage.end <= h.min) & (lineage.start >= h.max) ) {
#           pr=1
#         }
#         # if lineage goes extinct within the horizon
#         else if (lineage.end >= h.min) {
#           pr=h.max-lineage.end
#           pr=pr/s1
#         }
#         # if lineage originates within the horizon
#         else {
#           pr=lineage.start-h.min
#           pr=pr/s1
#         }
#
#         # 3. define the probabilty
#         pr = pr * PA * exp( (-(current.depth-PD)**2) / (2 * (DT ** 2)) )
#         # the equivalent command in perl
#         # calculate the probability of collection
#         # my $p_collection = $PA * exp( (-($current_waterdepth-$PD)**2) / (2 * ($DT ** 2)) );
#
#         # 4. if random.number < pr { record fossil as collected }
#         if (random.number <= pr ) {
#           fossils<-rbind(fossils,data.frame(x=h,y=i))
#         }
#       }
#
#     }
#   }
#
#   return(fossils) # in this data frame x=horizon and y=lineage
#
#   # EOF
# }
