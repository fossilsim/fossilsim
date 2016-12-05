#' Simulate fossils under a Poisson sampling model
#'
#' @param tree Phylo object.
#' @param phi Sampling rate.
#' @param root.edge Boolean indicating tree inclues a root edge.
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
#' @param basin.age Maximum age of the oldest horizon.
#' @param strata Number of stratigraphic horizons.
#' @param sampling Probability of sampling/preservation.
#' @return dataframe of sampled fossils.
#' sp = edge labels. h = horizon labels (= max horizon age).
#' @examples
#' t<-ape::rtree(4)
#' ba<-basin.age(t,root.edge=F)
#' sim.fossils.unif(t,ba,5,0.5)
#' @keywords uniform fossil preseravtion
#' @export
sim.fossils.unif<-function(tree,basin.age,strata,sampling,root.edge=T,convert.rate=FALSE){
  tree<-tree
  basin.age<-basin.age
  strata<-strata
  sampling<-sampling
  convert.rate<-convert.rate # convert prability to rate and generate k fossils

  s1=basin.age/strata # horizon length (= max age of youngest horizon)
  horizons<-seq(s1, basin.age, length=strata)

  # poisson rate under constant model
  rate = -log(1-sampling)/(basin.age/strata)

  node.ages<-n.ages(tree)
  root=length(tree$tip.label)+1

  fossils<-data.frame(h=numeric(),sp=numeric())

  brl = 0 # record total branch length for debugging

  for(h in horizons){

    h.min<-h-s1
    h.max<-h

    for(i in tree$edge[,2]){ # internal nodes + tips

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

        # calculate the proportion of time during each horizon the lineage is extant
        # lineage speciates and goes extinct in interval h
        if((lineage.end > h.min) && (lineage.start < h.max)){
          pr = (lineage.start-lineage.end)/s1
          f.max = lineage.start
          f.min = lineage.end
        }
        # lineage goes extinct in interval h
        else if(lineage.end > h.min){
          pr = (h.max-lineage.end)/s1
          f.max = h.max
          f.min = lineage.end
        }
        # lineage speciates in interval h
        else if(lineage.start < h.max){
          pr = (lineage.start-h.min)/s1
          f.max = lineage.start
          f.min = h.min
        }
        # lineage is extant the entire duration of interval h
        else{
          pr = 1
          f.max = h.max
          f.min = h.min
        }

        brl = brl + (s1*pr)

        if(convert.rate){
          # generate k fossils from a poisson distribution
          k = rpois(1,rate*s1*pr)
          if(k > 0){
            for(j in 1:k){
              age = runif(1,f.min,f.max)
              fossils<-rbind(fossils,data.frame(h=age,sp=i))
            }
          }
        } else {
          # define the probabilty
          pr = pr * sampling
          # if random.number < pr { record fossil as collected }
          if (runif(1) <= pr) {
            fossils<-rbind(fossils,data.frame(h=h,sp=i))
          }
        }
      }
    }

    if(root.edge){

      lineage.start = max(node.ages)+tree$root.edge
      lineage.end = max(node.ages)

      if ( (lineage.start >= h.min) & (lineage.end <= h.max) ) {

        # calculate the proportion of time during each horizon the lineage is extant
        # lineage speciates and goes extinct in interval h
        if((lineage.end > h.min) && (lineage.start < h.max)){
          pr = (lineage.start-lineage.end)/s1
          f.max = lineage.start
          f.min = lineage.end
        }
        # lineage goes extinct in interval h
        else if(lineage.end > h.min){
          pr = (h.max-lineage.end)/s1
          f.max = h.max
          f.min = lineage.end
        }
        # lineage speciates in interval h
        else if(lineage.start < h.max){
          pr = (lineage.start-h.min)/s1
          f.max = lineage.start
          f.min = h.min
        }
        # lineage is extant the entire duration of interval h
        else{
          pr = 1
          f.max = h.max
          f.min = h.min
        }

        brl = brl + (s1*pr)

        if(convert.rate){
          # generate k fossils from a poisson distribution
          k = rpois(1,rate*s1*pr)
          if(k > 0){
            for(j in 1:k){
              age = runif(1,f.min,f.max)
              fossils<-rbind(fossils,data.frame(h=age,sp=i))
            }
          }
        } else {
          # define the probabilty
          pr = pr * sampling
          # if random.number < pr { record fossil as collected }
          if (runif(1) <= pr) {
            fossils<-rbind(fossils,data.frame(h=h,sp=i))
          }
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
#' @param depth Maximum water depth.
#' @param cycles Number of cycles (transgressions and regressions)
#' @return dataframe of sampled water depths.
#' @examples
#' wd<-sim.water.depth(100)
#' plot(wd, type="l")
#' @keywords non-uniform fossil preservation
#' @export
sim.water.depth<-function(strata,depth=2,cycles=2){

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

#' Simulate fossils under a non-uniform model of preservation (Holland, 1995)
#'
#' @param tree Phylo object.
#' @param basin.age Maximum age of the oldest horizon.
#' @param strata Number of stratigraphic horizons.
#' @param PA Peak adbundance parameter.
#' @param PD Preferred depth parameter.
#' @param DT Depth tolerance parameter.
#' @param root.edge Boolean indicating tree inclues a root edge.
#' @return dataframe of sampled fossils.
#' sp = edge labels. h = horizon labels (= max horizon age).
#' @keywords uniform fossil preseravtion
#' @export
sim.fossils.non.unif<-function(tree,basin.age,strata,profile,PA=.5,PD=.5,DT=.5,root.edge=T){
  tree<-tree
  basin.age<-basin.age
  strata<-strata
  profile<-profile
  PA<-PA
  PD<-PD
  PA<-PA
  root.edge<-root.edge

  s1=basin.age/strata # horizon length (= max age of youngest horizon)
  horizons<-seq(s1, basin.age, length=strata)

  node.ages<-n.ages(tree)
  root=length(tree$tip.label)+1

  fossils<-data.frame(h=numeric(),sp=numeric())

  depth.counter=0

  for(h in horizons){

    depth.counter=depth.counter+1
    current.depth=profile$y[depth.counter]

    h.min<-h-s1
    h.max<-h

    for(i in tree$edge[,2]){ # internal nodes + tips

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
        pr = pr * PA * exp( (-(current.depth-PD)**2) / (2 * (DT ** 2)) )

        # 4. if random.number < pr { record fossil as collected }
        if (random.number <= pr ) {
          fossils<-rbind(fossils,data.frame(h=h,sp=i))
        }
      }

    }

    if(root.edge){

      lineage.start = max(node.ages)+tree$root.edge
      lineage.end = max(node.ages)

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
        pr = pr * PA * exp( (-(current.depth-PD)**2) / (2 * (DT ** 2)) )

        # 4. if random.number < pr { record fossil as collected }
        if (random.number <= pr ) {
          fossils<-rbind(fossils,data.frame(h=h,sp=root))
        }
      }

    }

  }

  return(fossils) # in this data frame h=horizon and sp=lineage

  # EOF

}

#' Select a sensible basin age based tree age
#'
#' @param tree Phylo object.
#' @param root.edge Boolean indicating tree inclues a root edge.
#' @return basin age
#' @export
basin.age<-function(tree,root.edge=T){
  node.ages<-n.ages(tree)
  if(root.edge)
    ba = max(node.ages)+tree$root.edge
  else
    ba = max(node.ages)

  ba = round(ba,1)+0.1
  return(ba)
}

