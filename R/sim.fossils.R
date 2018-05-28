#' Simulate fossils under a Poisson sampling model
#'
#' Simulate fossils for a phylo (\code{tree}) or taxonomy object (\code{taxonomy}).
#' If both are specified, the function uses taxonomy.
#' If no taxonomic information is provided, the function assumes all speciation is symmetric (i.e. bifurcating, \code{beta = 1}).
#' A vector of rates can be specified to allow for rate variation across lineages.
#'
#' @param rate A single Poisson sampling rate or a vector of rates. If a vector is provided, each entry will apply to each unique species in the order they appear in the corresponding taxonomy object.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#'
#' @return An object of class fossils.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils using the tree
#' rate = 2
#' f = sim.fossils.poisson(rate, tree = t)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s = sim.taxonomy(t, 0.5, 1, 0.5)
#' f = sim.fossils.poisson(rate, taxonomy = s)
#' plot(f, t)
#'
#' # simulate fossils with rate variation across lineages
#' rate = runif(length(s$sp), min = 0, max = 5)
#' f = sim.fossils.poisson(rate, species = s)
#' plot(f, t)
#'
#' @keywords uniform preservation
#' @seealso \code{\link{sim.fossils.intervals}}, \code{\link{sim.fossils.non.unif.depth}}
#' @export
#'
#' @importFrom stats rpois runif
sim.fossils.poisson = function(rate, tree = NULL, taxonomy = NULL, root.edge = TRUE) {

  if(is.null(tree) && is.null(taxonomy))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
      stop("tree must be an object of class \"phylo\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(taxonomy))
    warning("tree and taxonomy both defined, using taxonomy")

  if(is.null(taxonomy) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(taxonomy) && !ape::is.rooted(tree))
    stop("tree must be rooted")

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  if(length(rate) > 1 && length(rate) != length(unique(species$sp)))
    stop("vector of rates provided that doesn't correspond to the number of species")
  else if(length(rate) == 1)
    rate = rep(rate, length(unique(species$sp)))

  # If TRUE use exact sampling times.
  # If FALSE hmin and hmax will equal the start and end times of the corresponding edge.
  use.exact.times = TRUE

  fdf = fossils()

  lineages = unique(taxonomy$sp)

  for (i in 1:length(lineages)){

    sp = lineages[i]
    start = taxonomy$start[which(taxonomy$sp == sp)][1]
    end = taxonomy$end[which(taxonomy$sp == sp)][1]
    origin = taxonomy$origin[which(taxonomy$sp == sp)][1]
    edges = taxonomy[which(taxonomy$sp == sp), ]

    blength = start - end

    # sample fossil numbers from the Poisson distribution
    rand = rpois(1, blength*rate[i])

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
#' Simulate fossils for a phylo (\code{tree}) or taxonomy object (\code{taxonomy}).
#' If both are specified, the function uses taxonomy.
#' If no taxonomic information is provided, the function assumes all speciation is symmetric (i.e. bifurcating, \code{beta = 1}).
#'
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
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
#' t = ape::rtree(6)
#'
#' # assign a max age based on tree height
#' max.age = basin.age(t)
#'
#' # simulate fossils using basin.age and strata & probabilities
#' strata = 4
#' probability = rep(0.7, 4)
#' f = sim.fossils.intervals(t, basin.age = max.age, strata = strata, probabilities = probability)
#' plot(f, t, strata = strata, show.strata = TRUE)
#'
#' # simulate fossils using interval.ages & rates
#' times = seq(0, max.age, length.out = 4)
#' rates = c(5, 0, 5)
#' f = sim.fossils.intervals(t, interval.ages = times, rates = rates)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s = sim.taxonomy(t, 0.1, 0.1, 1)
#' f = sim.fossils.intervals(taxonomy = s, interval.ages = times, rates = rates)
#' plot(f, t)
#'
#' @keywords uniform fossil preservation
#' @keywords non-uniform fossil preservation
#' @seealso \code{\link{sim.fossils.poisson}}, \code{\link{sim.fossils.non.unif.depth}}
#' @export
sim.fossils.intervals = function(tree = NULL, taxonomy = NULL,
                                interval.ages = NULL, basin.age = NULL, strata = NULL,
                                probabilities = NULL, rates = NULL,
                                root.edge = TRUE, use.exact.times = TRUE){

  if(is.null(tree) && is.null(taxonomy))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(taxonomy))
    warning("tree and taxonomy both defined, using taxonomy")

  if(is.null(taxonomy) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(taxonomy) && !ape::is.rooted(tree))
    stop("tree must be rooted")

  if(is.null(interval.ages) && (is.null(basin.age) || is.null(strata)))
    stop("Intervals need to be defined by specifying either interval.ages or basin.age and strata")
  if(!is.null(basin.age) && !is.null(strata)) {
    if(!is.null(interval.ages)) warning("Two interval definitions found, using interval.ages")
    else interval.ages <- seq(0, basin.age, length = strata + 1)
  }

  if(is.null(probabilities) && is.null(rates)) stop("Either rates or probabilities need to be specified")

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)
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

  if(is.null(taxonomy))
    taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)

  fdf = fossils()

  lineages = unique(taxonomy$sp)

  for (sp in lineages) {

    start = taxonomy$start[which(taxonomy$sp == sp)][1]
    end = taxonomy$end[which(taxonomy$sp == sp)][1]
    origin = taxonomy$origin[which(taxonomy$sp == sp)][1]
    edges = taxonomy[which(taxonomy$sp == sp), ]

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
#' This function uses a three parameter Gaussian model to simulate non-uniform fossil preservation along a specified phylogeny.
#' Preservation varies with respect to water depth, which is used as a proxy for changes in sedimentary environment.
#' The per interval probability of sampling is \deqn{P(collection) = PA e ^ (-(d - PD)^2 / 2*DT^2 ) }
#' where \emph{PA} is species peak abundance, \emph{PD} is preferred depth, \emph{DT} is depth tolerance and \emph{d} is current water depth.
#' \emph{PD} is the depth at which the species is most likely to be found and is equivalent to the mean of the distribution.
#' \emph{PA} is the probability of sampling an occurrence at this depth.
#' \emph{DT} is the potential of a species to be found at a range of depths and is equivalent to the standard deviation. \cr \cr
#' Non-uniform interval ages can be specified as a vector (\code{interval.ages}) or a uniform set of interval ages can be specified using
#' maximum interval age (\code{basin.age}) and the number of intervals (\code{strata}), where interval length \eqn{= basin.age/strata}. \cr \cr
#' Simulate fossils for a phylo (\code{tree}) or taxonomy object (\code{taxonomy}).
#' If both are specified, the function uses taxonomy.
#' If no taxonomic information is provided, the function assumes all speciation is symmetric (i.e. bifurcating, \code{beta = 1}).
#'
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param depth.profile Vector of relative water depth. The first number corresponds to the youngest interval. The length of the vector should be 1 less than the length of interval.ages.
#' @param PA Peak abundance parameter.
#' @param PD Preferred depth parameter.
#' @param DT Depth tolerance parameter.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#'
#' @return An object of class fossils, where \code{hmin} and \code{hmax} will equal the start and end times of the corresponding interval.
#'
#' @references
#' Holland, S.M. 1995. The stratigraphic distribution of fossils. Paleobiology 21: 92-109.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # assign a max age based on tree height
#' max.age = basin.age(t)
#'
#' # generate water depth profile
#' strata = 7
#' wd = sim.water.depth(strata)
#'
#' # simulate fossils using tree & basin.age and strata
#' f = sim.fossils.non.unif.depth(t, basin.age = max.age, strata = strata,
#' depth.profile = wd, PA = 1, PD = 0.5, DT = 1, use.rates = TRUE)
#' plot(f, t, show.proxy = TRUE, proxy.data = wd, strata = strata, show.strata = TRUE)
#'
#' # simulate fossils using taxonomy & interval.ages
#' s = sim.taxonomy(t, 0.1, 0.1, 1)
#' times = seq(0, max.age, length.out = strata + 1)
#' f = sim.fossils.non.unif.depth(taxonomy = s, interval.ages = times,
#'      depth.profile = wd, PA = 1, PD = 0.5, DT = 1)
#' plot(f,t)
#'
#' @keywords non-uniform fossil preseravtion
#' @seealso \code{\link{sim.fossils.poisson}}, \code{\link{sim.fossils.intervals}}
#' @export
sim.fossils.non.unif.depth = function(tree = NULL, taxonomy = NULL,
                                      interval.ages = NULL, basin.age = NULL, strata = NULL,
                                      depth.profile = NULL, PA = 0.5, PD = 0.5, DT = 0.5,
                                      root.edge = TRUE){

  if(is.null(tree) && is.null(taxonomy))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(taxonomy))
    warning("tree and taxonomy both defined, using taxonomy")

  if(is.null(taxonomy) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(taxonomy) && !ape::is.rooted(tree))
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

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  fdf = fossils()

  lineages = unique(taxonomy$sp)

  for (sp in lineages) {

    start = taxonomy$start[which(taxonomy$sp == sp)][1]
    end = taxonomy$end[which(taxonomy$sp == sp)][1]
    origin = taxonomy$origin[which(taxonomy$sp == sp)][1]
    edges = taxonomy[which(taxonomy$sp == sp), ]

    blength = start - end

    #possible intervals covered by taxonomy
    for (i in 1:(length(interval.ages) - 1)) {
      if(interval.ages[i+1] < end) next
      if(interval.ages[i] > start) break

      min.time = max(end, interval.ages[i])
      max.time = min(start, interval.ages[i+1])

      # scale the probability
      pr = probabilities[i] * (max.time - min.time)/(interval.ages[i+1] - interval.ages[i])
      # assign fossils to edges
      ages = runif(1, min.time, max.time)
      edge = sapply(ages, function(x) edges$edge[which(edges$edge.start > x & edges$edge.end < x)])
      # if random.number < pr { record fossil as collected during interval }
      if (runif(1) <= pr) {
        # use interval ages
        min.time = interval.ages[i]
        max.time = interval.ages[i+1]
        fdf <- rbind(fdf,data.frame(sp = sp, edge = edge, origin = origin, hmin = min.time, hmax = max.time, stringsAsFactors = F))
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
#' wd = sim.water.depth(strata)
#' plot(wd, type="l")
#' @keywords non-uniform fossil preservation
#' @export
sim.water.depth = function(strata, depth = 2, cycles = 2){

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

#' Simulate fossils recovery rates with variation across lineages
#'
#' Available models include autocorrelated rates, independent rates and a "jump" model in which rate changes are coincident with speciation events.
#' Under the \code{autocorrelated} rates model, rates evolve along lineages according to a Brownian motion process, where the strength of the relationship between ancestor and descendant rates is determined by the parameter \eqn{\nu} (\code{v}).
#' If \eqn{\nu} is small rates will be more similar between ancestor and descendants, and if \eqn{\nu} is zero all rates will be equal.
#' For a given species \eqn{i} with ancestor \eqn{j}, a new rate \eqn{\kappa_i} is drawn from a lognornal distribution with
#' \deqn{\kappa_i = LN( \kappa_j - (\sigma^2/2), \sigma)}
#' where \eqn{\sigma = \nu * t_i} and \eqn{t_i} is the lineage duration of the species.
#' This fossil recovery model has been described in Heath et al. (2014) and is equivalent to the autocorrelated relaxed clock model described in Kishino et al. (2001).
#' Under the \code{independent} rates model a new rate is drawn for each species from any valid user-specified distribution (\code{dist}).
#' Under the \code{jump} model rates change at each speciation event with a given probability (\code{jump.pr}) and new rates are drawn from any valid user-specified distribution (\code{dist}).
#'
#' @param rate Initial rate at the origin or root of the phylo or taxonomy object. Default = 1.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#' @param model Model used to simulate rate variation across lineages. Options include "autocorrelated" (default), "independent" or "jump".
#' @param v Brownian motion parameter \eqn{v} used in the autocorrelated rates model. Default = 0.01.
#' @param dist Distribution of rates used to draw new rates under the "independent" and "jump" models. This parameter is ignored if \code{model = "autocorrealted"}. The default is a uniform distribution with \emph{U(0, 2)}. The distribution function must return a single positive value.
#' @param jump.pr Probability that fossil recovery rate changes at speciation events. Default = 0.01.
#' @return A vector of rates.
#' Rates are output for each species in the order they appear in the corresponding taxonomy object.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate taxonomy
#' s = sim.taxonomy(t, 0.5, 1, 0.5)
#'
#' # simualte rates under the autocorrelated rates model
#' rate = 1
#' rates = sim.species.rates(rate = rate, taxonomy = s, v = 1)
#' f = sim.fossils.poisson(rates, species = s)
#' plot(f, t)
#'
#' # simualte rates under the independent rates model
#' dist = function() { rlnorm(1, log(rate), 1) }
#' rates = sim.species.rates(rate = rate, taxonomy = s, model = "independent", dist = dist)
#' f = sim.fossils.poisson(rates, species = s)
#' plot(f, t)
#'
#' # simualte rates under the jump model
#' rates = sim.species.rates(rate = rate, taxonomy = s, model = "jump", dist = dist, jump.pr = 0.1)
#' f = sim.fossils.poisson(rates, species = s)
#' plot(f, t)
#'
#' @references
#' Heath et al. 2014. The fossilized birth-death process for coherent calibration of divergence-time estimates. PNAS 111:E2957-E2966.\cr
#' Kishino et al. 2001. Performance of a divergence time Estimation method under a probabilistic model of rate evolution MBE 18:352–361.
#'
#' @export
sim.species.rates = function(rate = 1, tree = NULL, taxonomy = NULL, root.edge = TRUE,
                             model = "autocorrelated", v = 0.01,
                             dist = function(){runif(1,0,2)}, jump.pr = 0.01, return.rates = TRUE){

  if(is.null(tree) && is.null(taxonomy))
    stop("Specify phylo or taxonomy object")

  if(!is.null(tree) && !"phylo" %in% class(tree))
    stop("tree must be an object of class \"phylo\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("species must be an object of class \"taxonomy\"")

  if(!is.null(tree) && !is.null(taxonomy))
    warning("tree and species both defined, using species taxonomy")

  if(is.null(taxonomy) && is.null(tree$edge.length))
    stop("tree must have edge lengths")

  if(is.null(taxonomy) && !ape::is.rooted(tree))
    stop("tree must be rooted")

  if(model != "autocorrelated" && model != "independent" && model != "jump")
    stop("specify a valid model option = 'autocorrelated', 'independent' or 'jump'")

  if(!(jump.pr >= 0 & jump.pr <= 1))
    stop("jump.pr must be a probability between 0 and 1")

  if((model == "independent" || model == "jump") & ( length(dist()) != 1 || !(is.numeric(dist()))))
    stop("specify a valid distribution function that returns a single +ve value")

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  aux = function(sp, t, r) {

    # simulate fossils
    start = t$start[which(t$sp == sp)][1]
    end = t$end[which(t$sp == sp)][1]
    origin = t$origin[which(t$sp == sp)][1]
    edges = t[which(t$sp == sp), ]

    blength = start - end

    # generate a new rate for sp
    if(model == "autocorrelated"){
      # this follows the molecular clock model of Kishino et al 2001
      # and the preservation model described in Heath et al 2014 (supplementary material)
      r = rlnorm(1, mean = log(r) - ((blength*v)/2), sdlog = sqrt(blength*v))
    } else if (model == "jump") {
      if(runif(1) < jump.pr)
        r = dist()
    } else { # independent rates
      r = dist()
    }

    if(r < 0) stop("specify a valid distribution function that returns a single +ve value")

    t[which(t$sp == sp),]$rate = r

    # fetch descendants
    descendants = unique(t$sp[which(t$parent == sp)])

    if(length(descendants) == 0) {
      return(t)
    }

    for(i in descendants){
      t = aux(i, t, r)
    }
    return(t)
  }

  taxonomy$rate = NA

  root = unique(taxonomy$sp[which(taxonomy$parent == 0)])

  taxonomy = aux(root, taxonomy, rate)

  # extract unique rates
  rates = unique(cbind(taxonomy["sp"],taxonomy["rate"]))$rate

  return(rates)
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
#' t = ape::rtree(6)
#' basin.age(t, root.edge = FALSE)
#'
#' @export
basin.age = function(tree,root.edge=TRUE){
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
count.fossils = function(fossils){
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
count.fossils.binned = function(fossils, interval.ages){
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
assign.interval = function(intervals, t){

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
#' @param taxonomy Taxonomy object.
#'
#' @return An object of class fossils.
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate fossils using the tree
#' rate = 2
#' f = sim.fossils.poisson(rate, tree = t)
#' plot(f, t)
#'
#' # simulate fossils using taxonomy
#' s = sim.taxonomy(t, 0.5, 1, 0.5)
#' f = reconcile.fossils.taxonomy(f, s)
#' plot(f, t)
#'
#' @export
reconcile.fossils.taxonomy = function(fossils, taxonomy){

  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
    stop("taxonomy must be an object of class \"taxonomy\"")

  # if fossils contain edges not in taxonomy
  if(!all(fossils$edge %in% taxonomy$edge))
    stop("incompatible fossils and taxonomy: not all fossil edges found in taxonomy")

  if(!identical(fossils$hmin, fossils$hmax))
    stop("exact fossil sampling times must be specified to use this function (i.e. hmin = hmax)")

  if(attr(fossils,"from.taxonomy"))
    warning("fossils already assigned based on taxonomy")

  # for each fossil identify the edge
  for(i in 1:length(fossils$edge)){
    edge = fossils$edge[i]
    # identify the edges in the corresponding taxonomy obj
    j = which(taxonomy$edge == edge)
    if(length(j) == 1){
      # reassign species
      fossils$sp[i] = taxonomy$sp[j]
      fossils$origin[i] = taxonomy$origin[j]
    } else { # {more than one species is associated with the edge }
      age = fossils$hmin[i]
      edges = taxonomy$edge[j]
      j = which(taxonomy$edge %in% edges & taxonomy$start > age & taxonomy$end < age)
      # reassign species
      fossils$sp[i] = taxonomy$sp[j]
      fossils$origin[i] = taxonomy$origin[j]
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
# t = ape::rtree(4)
# # simulate fossils
# rate = 2
# f = sim.fossils.exponential(t, rate)
# plot(f, t)
# @keywords uniform preservation
#
#' @importFrom stats rexp
sim.fossils.exponential = function(tree,rate,root.edge=TRUE){

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

  fossils <- as.fossils(fossils, FALSE)
  return(fossils) # in this data frame h=fossil age and sp=lineage
  # EOF
}
