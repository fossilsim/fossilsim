#' Simulate fossils under a Poisson sampling model
#'
#' @param tree Phylo object.
#' @param sampling Poisson sampling rate.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return An object of class fossils.
#' sp = edge labels. h = ages.
#' @examples
#' # simulate tree
#' t<-ape::rtree(4)
#' # simulate fossils
#' sampling = 2
#' f<-sim.fossils.poisson(t, sampling)
#' plot(f, t)
#' @keywords uniform preseravtion
#' @export
sim.fossils.poisson<-function(tree,sampling,root.edge=TRUE){
  tree<-tree
  lambda<-sampling
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
  fossils<-fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  return(fossils) # in this data frame h=fossil age and sp=lineage
  # EOF
}

#' Simulate fossils under a uniform model of preservation for a set of equal length intervals
#'
#' @param tree Phylo object.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param sampling Probability of sampling/preservation.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @param convert.rate If TRUE convert per interval sampling probability into a per interval Poisson rate (default = FALSE).
#' @return An object of class fossils.
#' sp = edge labels. h = fossil or interval ages. If convert.rate = TRUE, h = specimen age, if convert.rate = FALSE, h = max interval age.
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # assign a max age based on tree height
#' max<-basin.age(t)
#' # simulate fossils
#' strata = 4
#' sampling = 0.7
#' f<-sim.fossils.unif(t, max, strata, sampling)
#' plot(f, t, binned = TRUE, strata = strata)
#' @keywords uniform fossil preseravtion
#' @export
sim.fossils.unif<-function(tree,basin.age,strata,sampling,root.edge=T,convert.rate=FALSE){
  tree<-tree
  basin.age<-basin.age
  strata<-strata
  sampling<-sampling
  convert.rate<-convert.rate # convert prability to rate and generate k fossils

  if(!((sampling >= 0) & (sampling <= 1)))
      stop("Sampling must be a probability between 0 and 1")

  if(sampling == 1)
    sampling = 0.9999

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

    if(root.edge && exists("root.edge",tree) ){

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
              fossils<-rbind(fossils,data.frame(h=age,sp=root))
            }
          }
        } else {
          # define the probabilty
          pr = pr * sampling
          # if random.number < pr { record fossil as collected }
          if (runif(1) <= pr) {
            fossils<-rbind(fossils,data.frame(h=h,sp=root))
          }
        }
      }
    }
  }

  if(convert.rate)
    fossils<-fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  else
    fossils<-fossils(fossils, age = "interval.max", speciation.mode = "symmetric")
  return(fossils) # in this data frame h=horizon and sp=lineage
  # EOF
}

#' Simulate fossils under a non-uniform model of preservation
#'
#' @param tree Phylo object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param sampling Vector of Poisson sampling rates. The first number corresponds to the youngest interval. The length of the vector should 1 less than the length of interval.ages.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return An object of class fossils.
#' sp = edge labels. h = fossil ages.
#' @keywords non-uniform fossil preseravtion
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # assign a max age based on tree height
#' max = basin.age(t)
#' # assign interval times & rates
#' times = seq(0, max, length.out = 4)
#' rates = c(5, 3, 1)
#' # simulate fossils
#' f<-sim.fossils.non.unif(t, times, rates)
#' plot(f, t)
#' @export
sim.fossils.non.unif<-function(tree, interval.ages, sampling, root.edge = TRUE){
  tree<-tree
  interval.ages<-interval.ages
  rate<-sampling
  root.edge<-root.edge

  if(length(rate) != (length(interval.ages) - 1 ))
    stop("something went wrong when you specified the inteval rate and times vectors")

  horizons.min = head(interval.ages, -1)
  horizons.max = interval.ages[-1]

  node.ages<-n.ages(tree)
  root = length(tree$tip.label)+1

  fossils<-data.frame(h=numeric(),sp=numeric())

  brl = 0 # record total branch length for debugging

  for(h in 1:length(horizons.min)){

    h.min<-horizons.min[h]
    h.max<-horizons.max[h]
    s1 = h.max - h.min # horizon length

    if(rate[h]==0)
      next

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

        brl = brl + (s1*pr) # debugging code

        # generate k fossils from a poisson distribution
        k = rpois(1,rate[h]*s1*pr)
        if(k > 0){
          for(j in 1:k){
            age = runif(1,f.min,f.max)
            fossils<-rbind(fossils,data.frame(h=age,sp=i))
          }
        }
      }
    } # end of lineage

    if(root.edge && exists("root.edge",tree) ){

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

        # generate k fossils from a poisson distribution
        k = rpois(1,rate[h]*s1*pr)
        if(k > 0){
          for(j in 1:k){
            age = runif(1,f.min,f.max)
            fossils<-rbind(fossils,data.frame(h=age,sp=root))
          }
        }
      }
    } # end of root edge

  } # end of horizon

  fossils<-fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  return(fossils)
}

#' Simulate fossils under a non-uniform model of preservation (Holland, 1995)
#'
#' @description
#' This function uses a three parameter Guassion model to simulate non-uniform fossil preservation along a specified phylogeny.
#' Preservation varies with respect to water depth, which is used as a proxy for changes in sedimentary environment.
#' The per interval probability of sampling is \deqn{P(collection) = PA e ^ (-(d - PD)^2 / 2*DT^2 ) }
#' where \emph{PA} is species peak abundance, \emph{PD} is preferred depth, \emph{DT} is depth tolerance and \emph{d} is current water depth.
#' \emph{PD} is the depth at which the species is most likely to be found and is equivalent to the mean of the distribution.
#' \emph{PA} is the probability of sampling an occurrence at this depth.
#' \emph{DT} is the potential of a species to be found at a range of depths and is equivalent to the standard deviation. \cr \cr
#' Non-uniform interval ages can be specified as a vector (\code{interval.ages}) or a uniform set of interval ages can be specified using
#' maximum interval age (\code{basin.age}) and the number of intervals (\code{strata}), where interval length \eqn{= basin.age/strata}.
#'
#' @param tree Phylo object.
#' @param profile Vector of relative water depth. The first number corresponds to the youngest interval. The length of the vector should 1 less than the length of interval.ages.
#' @param PA Peak adbundance parameter.
#' @param PD Preferred depth parameter.
#' @param DT Depth tolerance parameter.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @param convert.rate If TRUE convert per interval sampling probability into a per interval Poisson rate (default = FALSE).
#' @return An object of class fossils.
#' sp = edge labels. h = fossil or interval ages. If convert.rate = TRUE, h = specimen age, if convert.rate = FALSE, h = max horizon age.
#'
#' @references
#' Holland, S.M. 1995. The stratigraphic distribution of fossils. Paleobiology 21: 92–109.
#'
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#' # assign a max age based on tree height
#' max = basin.age(t)
#' # generate water depth profile
#' strata = 7
#' wd<-sim.water.depth(strata)
#' # simulate fossils
#' f<-sim.fossils.non.unif.depth(t, wd, PA = 1, PD = 0.5, DT = 1, basin.age = max, strata = strata, convert.rate = TRUE)
#' plot(f,t, show.proxy = T, proxy.data = wd, strata = strata, show.strata = T)
#' @keywords non-uniform fossil preseravtion
#' @export
sim.fossils.non.unif.depth<-function(tree, profile, PA=.5, PD=.5, DT=.5, interval.ages = NULL, basin.age = NULL, strata = NULL, root.edge = TRUE, convert.rate = FALSE){
  tree<-tree
  profile<-profile
  PA<-PA
  PD<-PD
  DT<-DT
  interval.ages<-interval.ages
  basin.age<-basin.age
  strata<-strata
  root.edge<-root.edge
  convert.rate<-convert.rate

  if( (is.null(interval.ages)) && (is.null(basin.age)) && (is.null(strata)) )
    stop("Specify interval.ages OR basin.age and number of strata")
  else if ( (is.null(interval.ages)) && (is.null(basin.age)) )
    stop("Specify interval.ages OR basin.age and number of strata")
  else if ( (is.null(interval.ages)) && (is.null(strata)) )
    stop("Specify interval.ages OR basin.age and number of strata")
  # add warning about conflicting info

  if( (is.null(interval.ages)) ){
    s1 = basin.age/strata # horizon length (= max age of youngest horizon)
    horizons.max = seq(s1, basin.age, length = strata)
    horizons.min = horizons.max - s1
  } else {
    horizons.min = head(interval.ages, -1)
    horizons.max = interval.ages[-1]
  }

  if(length(wd) < length(horizons.max))
    stop("Water depth values < the number of intervals!")

  node.ages<-n.ages(tree)
  root=length(tree$tip.label)+1

  fossils<-data.frame(h=numeric(),sp=numeric())

  brl = 0 # record total branch length for debugging

  depth.counter = 0

  for(h in 1:length(horizons.min)){

    h.min<-horizons.min[h]
    h.max<-horizons.max[h]
    s1 = h.max - h.min # horizon length

    current.depth = profile[h]

    # calculate the interval rate
    sampling = PA * exp( (-(current.depth-PD)**2) / (2 * (DT ** 2)) )
    if(sampling >= 1)
      sampling = 0.9999
    rate = -log(1-sampling)/s1

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
            fossils<-rbind(fossils,data.frame(h=h.max,sp=i))
          }
        }
      }
    }

    if(root.edge && exists("root.edge",tree) ){

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
              fossils<-rbind(fossils,data.frame(h=age,sp=root))
            }
          }
        } else {
          # define the probabilty
          pr = pr * sampling
          # if random.number < pr { record fossil as collected }
          if (runif(1) <= pr) {
            fossils<-rbind(fossils,data.frame(h=h,sp=root))
          }
        }
      }
    }
    #eol
  }

  if(convert.rate)
    fossils<-fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  else
    fossils<-fossils(fossils, age = "interval.max", speciation.mode = "symmetric")
  return(fossils) # in this data frame h=horizon and sp=lineage

  #eof
}

#' Simulate water depth profile
#'
#' @description
#' Function returns water depth profile using the sine wave function \eqn{y = depth*sin(cycles*pi*(x-1/4))}.
#'
#' @param strata Number of stratigraphic intervals
#' @param depth Maximum water depth.
#' @param cycles Number of cycles (transgressions and regressions).
#' @return dataframe of sampled water depths.
#' @examples
#' strata = 100
#' wd<-sim.water.depth(strata)
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

  #return(data.frame(x=c(1:strata),y=y))
  return(y)

  # EOF
}

#' Select a sensible basin age based on tree height
#'
#' @description
#' Function returns an age slightly older than the root.age or origin time.
#'
#' @param tree Phylo object.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return basin age
#' @examples
#' t<-ape::rtree(6)
#' basin.age(t, root.edge = FALSE)
#' @export
basin.age<-function(tree,root.edge=TRUE){
  node.ages<-n.ages(tree)
  if(root.edge && exists("root.edge",tree) )
    ba = max(node.ages)+tree$root.edge
  else
    ba = max(node.ages)

  ba = round(ba,1)+0.1
  return(ba)
}

#' Count the total number of fossils
#'
#' @param fossils Fossils object.
#' @return Number of extinct samples.
#'
#' @export
count.fossils<-function(fossils){
  k = length(fossils$sp[which(fossils$h > 0)])
  return(k)
}

#' Count the total number of fossils per interval
#'
#' @param fossils Fossils object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#'
#' @return Vector of extinct samples corresponding to each interval. Note the last value corresponds to the number of samples > the maximum age of the oldest interval.
#'
#' @export
count.fossils.binned<-function(fossils, interval.ages){
  intervals<-interval.ages
  k = c()
  for(i in 1:length(intervals)){
    k = c(k, 0)
  }
  for(i in 1:length(fossils$h)){
    if(fossils$h[i] != 0){
      j = assign.interval(intervals, fossils$h[i])
      k[j] = k[j] + 1
    }
  }
  return(k)
}

# assign any given age to one of a set of intervals
assign.interval<-function(intervals, t){
  i = -1
  for(j in 1:length(intervals)){
    if(t >= intervals[j])
      i = j
  }
  return(i)
}

