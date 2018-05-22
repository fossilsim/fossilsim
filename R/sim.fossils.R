#' Simulate fossils under a Poisson sampling model
#'
#' Simulate fossils for a phylo (\code{tree}) or taxonomy object (\code{species}).
#' If both are specified, the function uses taxonomy.
#' If no taxonomic information is provided, the function assumes all speciation is symmetric (i.e. budding, \code{beta = 1}).
#'
#' @param rate Poisson sampling rate.
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#'
#' @return An object of class fossils.
#'
#' @examples
#' # simulate tree
#' t <- ape::rtree(6)
#'
#' # simulate fossils using the tree
#' rate = 2
#' f <- sim.fossils.poisson(rate, tree = t)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s <- create.taxonomy(t, 0.5, 1, 0.5)
#' f <- sim.fossils.poisson(rate, species = s)
#' plot(f, t)
#'
#' @keywords uniform preservation
#' @export
#'
#' @importFrom stats rpois runif
sim.fossils.poisson<-function(rate, tree = NULL, species = NULL, root.edge = TRUE) {

  if(is.null(tree) && is.null(species))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
      stop("tree must be an object of class \"phylo\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  if(is.null(species) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(species) && !ape::is.rooted(tree))
    stop("tree must be rooted")

  if(is.null(species)){
    species = create.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  # If TRUE use exact sampling times.
  # If FALSE hmin and hmax will equal the start and end times of the corresponding edge.
  use.exact.times = TRUE

  fdf = fossils()

  lineages = unique(species$sp)

  for (sp in lineages){ # internal nodes + tips

    start = species$start[which(species$sp == sp)][1]
    end = species$end[which(species$sp == sp)][1]
    origin = species$origin[which(species$sp == sp)][1]
    edges = species[which(species$sp == sp), ]

    blength = start - end

    # sample fossil numbers from the Poisson distribution
    rand = rpois(1, blength*rate)

    if(rand > 0) {
      h = runif(rand, min = end, max = start)
      edge = sapply(h, function(x) edges$edge[which(edges$edge.start > x & edges$edge.end < x)])
      if(use.exact.times) {
        fdf <- rbind(fdf, data.frame(sp = sp, edge = edge, origin = origin, hmin = h, hmax = h, stringsAsFactors = F))
      } else {
        fdf <- rbind(fdf, data.frame(sp = sp, edge = edge, origin = origin, hmin = rep(end, rand), hmax = rep(start, rand), stringsAsFactors = F))
      }
    }
  }
  fdf <- as.fossils(fdf, from.taxonomy)
  return(fdf)
}

#' Simulate fossils under a non-uniform model of preservation for a given set of consecutive time intervals
#'
#' Intervals can be specified by specifying the interval boundaries using \code{interval.ages} or specifying both \code{basin.age} and \code{strata}.
#' In the second scenario all intervals will be of equal length.
#' Preservation can be specified using \code{rates}, which represent the rates of a Poisson process in each interval,
#' or \code{probabilities}, which represent the probabilities of sampling per interval.
#' When using \code{probabilities}, at most one fossil per species will be sampled per interval. \cr \cr
#' Simulate fossils for a phylo (\code{tree}) or taxonomy object (\code{species}).
#' If both are specified, the function uses taxonomy.
#' If no taxonomic information is provided, the function assumes all speciation is symmetric (i.e. budding, \code{beta = 1}).
#'
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param rates Poisson sampling rate for each interval. The number of rates should match the number of intervals.
#' @param probabilities Probability of sampling/preservation in each interval. The number of probabilities should match the number of intervals.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#' @param use.exact.times If TRUE use exact sampling times. If FALSE \code{hmin} and \code{hmax} will equal the start and end times of the corresponding interval. Default = TRUE.
#' @return An object of class fossils.
#'
#' @examples
#' # simulate tree
#' t <- ape::rtree(6)
#'
#' # assign a max age based on tree height
#' max.age <- basin.age(t)
#'
#' # simulate fossils using basin.age and strata & probabilities
#' strata = 4
#' probability = rep(0.7, 4)
#' f <- sim.fossils.intervals(t, basin.age = max.age, strata = strata, probabilities = probability)
#' plot(f, t, strata = strata, show.strata = TRUE)
#'
#' # simulate fossils using interval.ages & rates
#' times = seq(0, max.age, length.out = 4)
#' rates = c(5, 0, 5)
#' f <- sim.fossils.intervals(t, interval.ages = times, rates = rates)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s <- create.taxonomy(t, 0.1, 0.1, 1)
#' f <- sim.fossils.intervals(species = s, interval.ages = times, rates = rates)
#' plot(f, t)
#'
#' @keywords uniform fossil preservation
#' @keywords non-uniform fossil preservation
#' @export
sim.fossils.intervals<-function(tree = NULL, species = NULL,
                                interval.ages = NULL, basin.age = NULL, strata = NULL,
                                probabilities = NULL, rates = NULL,
                                root.edge = TRUE, use.exact.times = TRUE){

  if(is.null(tree) && is.null(species))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  if(is.null(species) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(species) && !ape::is.rooted(tree))
    stop("tree must be rooted")

  if(is.null(interval.ages) && (is.null(basin.age) || is.null(strata)))
    stop("Intervals need to be defined by specifying either interval.ages or basin.age and strata")
  if(!is.null(basin.age) && !is.null(strata)) {
    if(!is.null(interval.ages)) warning("Two interval definitions found, using interval.ages")
    else interval.ages <- seq(0, basin.age, length = strata + 1)
  }

  if(is.null(probabilities) && is.null(rates)) stop("Either rates or probabilities need to be specified")

  if(is.null(species)){
    species = create.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  use.rates = FALSE
  if(!is.null(probabilities) && !is.null(rates)) {
    rates = NULL
    warning("Both probabilities and rates found, using probabilities")
  }
  if(!is.null(rates)) {
    use.rates = TRUE
    if(length(rates) != (length(interval.ages) - 1 )) stop("Length mismatch between interval ages and sampling rates")
  } else {
    if(length(probabilities) != (length(interval.ages) - 1 )) stop("Length mismatch between interval ages and sampling probabilities")
    if(any(probabilities < 0) || any(probabilities > 1)) stop("Sampling probabilities must be between 0 and 1")
  }

  if(is.null(species))
    species = create.taxonomy(tree, beta = 1, root.edge = root.edge)

  fdf = fossils()

  lineages = unique(species$sp)

  for (sp in lineages) { # internal nodes + tips

    start = species$start[which(species$sp == sp)][1]
    end = species$end[which(species$sp == sp)][1]
    origin = species$origin[which(species$sp == sp)][1]
    edges = species[which(species$sp == sp), ]

    blength = start - end

    #possible intervals covered by species
    for (i in 1:(length(interval.ages) - 1)) {
      if(interval.ages[i+1] < end) next
      if(interval.ages[i] > start) break

      min.time = max(end, interval.ages[i])
      max.time = min(start, interval.ages[i+1])

      if(use.rates) {
        # generate k fossils from a poisson distribution
        k = rpois(1, rates[i]*(max.time - min.time))
        ages = runif(k, min.time, max.time)
        edge = sapply(ages, function(x) edges$edge[which(edges$edge.start > x & edges$edge.end < x)])
        if(k > 0) {
          if(use.exact.times) {
            fdf <- rbind(fdf, data.frame(sp = sp, edge = edge, origin = origin, hmin = ages, hmax = ages, stringsAsFactors = F))
          } else {
            min.time = rep(interval.ages[i], k)
            max.time = rep(interval.ages[i+1], k) # this is kind of redundant
            fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = min.time, hmax = max.time, stringsAsFactors = F))
          }
        }
      } else {
        # scale the probability
        pr = probabilities[i] * (max.time - min.time)/(interval.ages[i+1] - interval.ages[i])
        ages = runif(1, min.time, max.time)
        edge = sapply(ages, function(x) edges$edge[which(edges$edge.start > x & edges$edge.end < x)])
        # if random.number < pr { record fossil as collected during interval }
        if (runif(1) <= pr) {
          if(use.exact.times) {
            fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = ages, hmax = ages, stringsAsFactors = F))
          } else { # { use interval ages }
            min.time = interval.ages[i]
            max.time = interval.ages[i+1]
            fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = min.time, hmax = max.time, stringsAsFactors = F))
          }
        }
      }
    }
  }
  fdf <- as.fossils(fdf, from.taxonomy)
  return(fdf)
}

#' Simulate fossils under a non-uniform model of preservation (Holland, 1995)
#'
#' @description
#' This function uses a three parameter Guassian model to simulate non-uniform fossil preservation along a specified phylogeny.
#' Preservation varies with respect to water depth, which is used as a proxy for changes in sedimentary environment.
#' The per interval probability of sampling is \deqn{P(collection) = PA e ^ (-(d - PD)^2 / 2*DT^2 ) }
#' where \emph{PA} is species peak abundance, \emph{PD} is preferred depth, \emph{DT} is depth tolerance and \emph{d} is current water depth.
#' \emph{PD} is the depth at which the species is most likely to be found and is equivalent to the mean of the distribution.
#' \emph{PA} is the probability of sampling an occurrence at this depth.
#' \emph{DT} is the potential of a species to be found at a range of depths and is equivalent to the standard deviation. \cr \cr
#' Non-uniform interval ages can be specified as a vector (\code{interval.ages}) or a uniform set of interval ages can be specified using
#' maximum interval age (\code{basin.age}) and the number of intervals (\code{strata}), where interval length \eqn{= basin.age/strata}. \cr \cr
#' Simulate fossils for a phylo (\code{tree}) or taxonomy object (\code{species}).
#' If both are specified, the function uses taxonomy.
#' If no taxonomic information is provided, the function assumes all speciation is symmetric (i.e. budding, \code{beta = 1}).
#'
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param depth.profile Vector of relative water depth. The first number corresponds to the youngest interval. The length of the vector should be 1 less than the length of interval.ages.
#' @param PA Peak adbundance parameter.
#' @param PD Preferred depth parameter.
#' @param DT Depth tolerance parameter.
#' @param use.rates If TRUE convert per interval sampling probability into a per interval Poisson rate. Default = FALSE.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#' @param use.exact.times If TRUE use exact sampling times. If FALSE \code{hmin} and \code{hmax} will equal the start and end times of the corresponding interval. Default = TRUE.
#'
#' @return An object of class fossils.
#'
#' @references
#' Holland, S.M. 1995. The stratigraphic distribution of fossils. Paleobiology 21: 92-109.
#'
#' @examples
#' # simulate tree
#' t<-ape::rtree(6)
#'
#' # assign a max age based on tree height
#' max.age = basin.age(t)
#'
#' # generate water depth profile
#' strata = 7
#' wd<-sim.water.depth(strata)
#'
#' # simulate fossils using tree & basin.age and strata
#' f<-sim.fossils.non.unif.depth(t, basin.age = max.age, strata = strata,
#'  depth.profile = wd, PA = 1, PD = 0.5, DT = 1, use.rates = TRUE)
#' plot(f,t, show.proxy = TRUE, proxy.data = wd, strata = strata, show.strata = TRUE)
#'
#' # simulate fossils using taxonomy & interval.ages
#' s <- create.taxonomy(t, 0.1, 0.1, 1)
#' times = seq(0, max.age, length.out = strata + 1)
#' f <- sim.fossils.non.unif.depth(species = s, interval.ages = times,
#'      depth.profile = wd, PA = 1, PD = 0.5, DT = 1, use.rates = TRUE)
#' plot(f,t)
#'
#' @keywords non-uniform fossil preseravtion
#' @export
sim.fossils.non.unif.depth<-function(tree = NULL, species = NULL,
                                interval.ages = NULL, basin.age = NULL, strata = NULL,
                                depth.profile = NULL, PA = 0.5, PD = 0.5, DT = 0.5, use.rates = FALSE,
                                root.edge = TRUE, use.exact.times = TRUE){

  if(is.null(tree) && is.null(species))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(species))
    warning("tree and species both defined, using species taxonomy")

  if(is.null(species) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(species) && !ape::is.rooted(tree))
    stop("tree must be rooted")

  if(is.null(interval.ages) && (is.null(basin.age) || is.null(strata)))
    stop("Intervals need to be defined by specifying either interval.ages or basin.age and strata")
  if(!is.null(basin.age) && !is.null(strata)) {
    if(!is.null(interval.ages)) warning("Two interval definitions found, using interval.ages")
    else interval.ages <- seq(0, basin.age, length = strata + 1)
  }

  if(is.null(depth.profile)) stop("No water depth profile specified")
  if(length(depth.profile) != (length(interval.ages)-1))
    stop("Mismatch between the number of intervals and depth profile values")

  # calculate per interval probabilities
  probabilities = sapply(depth.profile, function(x) {PA * exp( (-(x-PD)**2) / (2 * (DT ** 2)) )})

  #TODO: still not sure if this is appropriate
  if(use.rates){
    s = sapply(1:length(interval.ages[-1]), function(x) { interval.ages[x+1] - interval.ages[x] })
    if(any(probabilities >= 1)){
      probabilities[which(probabilities >= 1)] = 0.99999
    }
    rates = -log(1-probabilities)/s
  }

  if(is.null(species)){
    species = create.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  fdf = fossils()

  lineages = unique(species$sp)

  for (sp in lineages) { # internal nodes + tips

    start = species$start[which(species$sp == sp)][1]
    end = species$end[which(species$sp == sp)][1]
    origin = species$origin[which(species$sp == sp)][1]
    edges = species[which(species$sp == sp), ]

    blength = start - end

    #possible intervals covered by species
    for (i in 1:(length(interval.ages) - 1)) {
      if(interval.ages[i+1] < end) next
      if(interval.ages[i] > start) break

      min.time = max(end, interval.ages[i])
      max.time = min(start, interval.ages[i+1])

      if(use.rates) {
        # generate k fossils from a poisson distribution
        k = rpois(1, rates[i]*(max.time - min.time))
        ages = runif(k, min.time, max.time)
        edge = sapply(ages, function(x) edges$edge[which(edges$edge.start > x & edges$edge.end < x)])
        if(k > 0) {
          if(use.exact.times) {
            fdf <- rbind(fdf, data.frame(sp = sp, edge = edge, origin = origin, hmin = ages, hmax = ages, stringsAsFactors = F))
          } else { # { use interval ages }
            min.time = rep(interval.ages[i], k)
            max.time = rep(interval.ages[i+1], k) # this is kind of redundant
            fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = min.time, hmax = max.time, stringsAsFactors = F))
          }
        }
      } else {
        # scale the probability
        pr = probabilities[i] * (max.time - min.time)/(interval.ages[i+1] - interval.ages[i])
        ages = runif(1, min.time, max.time)
        edge = sapply(ages, function(x) edges$edge[which(edges$edge.start > x & edges$edge.end < x)])
        # if random.number < pr { record fossil as collected during interval }
        if (runif(1) <= pr) {
          if(use.exact.times) {
            fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = ages, hmax = ages, stringsAsFactors = F))
          } else { # { use interval ages }
            min.time = interval.ages[i]
            max.time = interval.ages[i+1]
            fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = min.time, hmax = max.time, stringsAsFactors = F))
          }
        }
      }
    }
  }
  fdf <- as.fossils(fdf, from.taxonomy)
  return(fdf)
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
  x = seq(0,2,length.out=strata)

  # define y-axis values
  # a - total depth excursion - amplitude
  # b - number of cycles
  # 1/c - defines the relative start time of each cycle - phase shift
  # y = a * sin (b * pi * (x-1/c))
  y = depth*sin(cycles*pi*(x-1/4))

  #return(data.frame(x=c(1:strata),y=y))
  return(y)

  # EOF
}

#' Define a basin age based on tree height
#'
#' @description
#' Function returns an age slightly older than the root.age or origin time using the formula
#' \eqn{round(max,1) + 0.1}, where max is the root.age or origin time.
#'
#' @param tree Phylo object.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return basin age
#' @examples
#' t<-ape::rtree(6)
#' basin.age(t, root.edge = FALSE)
#'
#' @export
basin.age<-function(tree,root.edge=TRUE){
  node.ages<-n.ages(tree)
  if(root.edge && exists("root.edge",tree) )
    ba = max(node.ages) + tree$root.edge
  else
    ba = max(node.ages)

  ba = round(ba,1) + 0.1
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

  k = rep(0, length(intervals))

  if(length(fossils$sp) == 0)
    return(k)

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

  if(is.null(intervals) || is.null(t))
     stop("specify intervals and time t")

  if(any(intervals < intervals[1]))
    stop("specify intervals from youngest to oldest")

  i = -1
  for(j in 1:length(intervals)){
    if(t >= intervals[j])
      i = j
  }
  return(i)
}

#' Reconcile existing fossil and taxonomy objects
#'
#' This function uses edge identifiers (\code{edge}) and fossil sampling times (\code{hmin}) to reassign fossil species identifiers (\code{sp, origin}) using an existing taxonomy object.
#' It can only be used if exact fossil sampling times are known (i.e. \code{hmin = hmax}), otherwise edges containing multiple species may be indistinguishable.
#'
#' @param fossils Fossils object.
#' @param species Taxonomy object.
#'
#' @return An object of class fossils.
#' @examples
#' # simulate tree
#' t <- ape::rtree(6)
#'
#' # simulate fossils using the tree
#' rate = 2
#' f <- sim.fossils.poisson(rate, tree = t)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s <- create.taxonomy(t, 0.5, 1, 0.5)
#' f <- reconcile.fossils.taxonomy(f, s)
#' plot(f, t)
#'
#' @export
reconcile.fossils.taxonomy<-function(fossils, species){

  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(!is.null(species) && !"taxonomy" %in% class(species))
    stop("species must be an object of class \"taxonomy\"")

  # if fossils contain edges not in species
  if(!all(fossils$edge %in% species$edge))
    stop("incompatible fossils and taxonomy: not all fossil edges found in species")

  if(!identical(fossils$hmin, fossils$hmax))
    stop("exact fossil sampling times must be specified to use this function (i.e. hmin = hmax)")

  if(attr(fossils,"from.taxonomy"))
    warning("fossils already assigned based on taxonomy")

  # for each fossil identify the edge
  for(i in 1:length(fossils$edge)){
    edge = fossils$edge[i]
    # identify the edges in the corresponding taxonomy obj
    j = which(species$edge == edge)
    if(length(j) == 1){
      # reassign species
      fossils$sp[i] = species$sp[j]
      fossils$origin[i] = species$origin[j]
    } else { # {more than one species is associated with the edge }
      age = fossils$hmin[i]
      edges = species$edge[j]
      j = which(species$edge %in% edges & species$start > age & species$end < age)
      # reassign species
      fossils$sp[i] = species$sp[j]
      fossils$origin[i] = species$origin[j]
    }
  }
  fossils = as.fossils(fossils, from.taxonomy = TRUE)
  return(fossils)
}

# Simulate fossils under an exponential sampling model
#
# @param tree Phylo object.
# @param rate Exponential sampling rate.
# @param root.edge If TRUE include the root edge (default = TRUE).
# @return An object of class fossils.
# sp = node labels. h = ages.
# The label is for the node just below the sampled fossil.
# @examples
# # simulate tree
# t<-ape::rtree(4)
# # simulate fossils
# rate = 2
# f<-sim.fossils.exponential(t, rate)
# plot(f, t)
# @keywords uniform preservation
#
# @importFrom stats rexp
sim.fossils.exponential<-function(tree,rate,root.edge=TRUE){

  node.ages<-n.ages(tree)

  fossils<-data.frame(h=numeric(),sp=numeric())

  root = length(tree$tip.label) + 1

  if(root.edge && exists("root.edge",tree) ){

    lineages = c(tree$edge[,2], root)

  } else lineages = tree$edge[,2]

  for (i in lineages){ # internal nodes + tips

    if(i == root){

      # root age
      a=which(names(node.ages)==root)
      lineage.end=node.ages[[a]]

      # origin time
      b=tree$root.edge
      lineage.start=lineage.end+b

    } else {

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
    }

    t = 0
    while(TRUE){
      t = t + rexp(1, rate);
      if (t < b) { # make fossil
        fossils<-rbind(fossils, data.frame(h=(lineage.start-t),sp=i))
      }
      else break
    }
  }

  fossils <- as.fossils(fossils, from.taxonomy)
  return(fossils) # in this data frame h=fossil age and sp=lineage
  # EOF
}
