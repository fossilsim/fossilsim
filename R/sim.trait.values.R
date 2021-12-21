#' Simulate trait values with variation across lineages
#'
#' @description
#' Fossil recovery rates or other parameter values can be simulated for a phylo (\code{tree}) or taxonomy (\code{taxonomy}) object.
#' Under the \code{autocorrelated} model, trait values evolve along lineages according to a Brownian motion process, where the strength of the relationship between ancestor and descendant values is determined by the parameter \eqn{\nu} (\code{v}).
#' If \eqn{\nu} is small values will be more similar between ancestor and descendants, and if \eqn{\nu} is zero all trait values will be equal.
#' For a given species \eqn{i} with ancestor \eqn{j}, a new trait value \eqn{\kappa_i} is drawn from a lognormal distribution with
#' \deqn{\kappa_i ~ LN( ln([\kappa_j] - (\sigma^2/2), \sigma)}
#' where \eqn{\sigma = \nu * t_i} and \eqn{t_i} is the lineage duration of the species.
#' This fossil recovery model is described in Heath et al. (2014) and is equivalent to the autocorrelated relaxed clock model described in Kishino et al. (2001).
#' Under the \code{BM} and \code{OU} models traits are simulated under a standard Brownian motion or Ornstein-Uhlenbeck process with rate parameter \eqn{\nu} (\code{v}).
#' The OU model has the additional parameter \code{alpha}, which determines the strength with which trait values are attracted to the mean. Note the \code{init} argument will specify both the value at the root and the mean of the process under the OU model.
#' Under the \code{independent} model a new trait value is drawn for each species from any valid user-specified distribution (\code{dist}).
#' \code{change.pr} is the probability that a trait value will change at each speciation event.
#' If \code{change.pr = 1} trait values will be updated at every speciation events.
#' Finally, traits can be simulated under the standard Lewis Mk model (\code{Mk}), with symmetric rates of change. The rate is specified using \code{v} and number of states using \code{k}.
#'
#' @param init Initial value at the origin or root of the phylo or taxonomy object. Default = 1.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#' @param model Model used to simulate rate variation across lineages. Options include "autocorrelated" (default), "BM" (Brownian motion), "OU" (Ornstein-Uhlenbeck), "independent" or the Lewis "Mk" model.
#' @param v Brownian motion parameter \eqn{v} used in the autocorrelated, BM and OU models. Or rate change under the Mk model. Default = 0.01.
#' @param alpha Ornstein-Uhlenbeck parameter \eqn{alpha}. Determines the strength with which trait values are pulled back towards the mean.
#' @param dist Distribution of trait values used to draw new values under the "independent" model. This parameter is ignored if \code{model = "autocorrealted"}. The default is a uniform distribution with \emph{U(0, 2)}. The distribution function must return a single value.
#' @param change.pr Probability that trait values change at speciation events. Default = 1.
#' @param k Number of states used for the Mk model. Default = 2.
#' @return A vector of parameter values.
#' Values are output for each species in the order in which they appear in the taxonomy object (if taxonomy was provided) or for each edge in the order in which they appear in the tree object.
#' If the tree object has a root edge (\code{root.edge}), the first entry in the vector will correspond to this edge.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # simulate taxonomy
#' s = sim.taxonomy(t, 0.5, 1, 0.5)
#'
#' # simulate rates under the autocorrelated trait values model
#' rate = 2
#' rates = sim.trait.values(rate, taxonomy = s, v = 1)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, t)
#'
#' # simulate rates under the independent trait values model
#' dist = function() { rlnorm(1, log(rate), 1) }
#' rates = sim.trait.values(rate, taxonomy = s, model = "independent", dist = dist)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, t)
#'
#' # simulate rates under the independent trait values model with infrequent changes
#' rates = sim.trait.values(rate, taxonomy = s, model = "independent",
#'                         dist = dist, change.pr = 0.1)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, t)
#'
#' # simulate traits under Brownian motion and convert into rates
#' traits = sim.trait.values(0, taxonomy = s, model = "BM", v = 2)
#' function for translating states into rates
#' translate.states = function(traits, low, high) sapply(traits, function(t) if(t < 0) low else high)
#' sampling rates
#' low = 0.1
#' high = 2
#' rates = translate.states(traits, low, high)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, tree = t)
#'
#' @references
#' Heath et al. 2014. The fossilized birth-death process for coherent calibration of divergence-time estimates. PNAS 111:E2957-E2966.\cr
#' Kishino et al. 2001. Performance of a divergence time estimation method under a probabilistic model of rate evolution MBE 18:352-361.
#'
#' @export
sim.trait.values = function(init = 1, tree = NULL, taxonomy = NULL, root.edge = TRUE,
                             model = "autocorrelated", v = 0.01, alpha = 0.1,
                             dist = function(){runif(1,0,2)}, change.pr = 1, k = 2){

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

  if(model != "autocorrelated" && model != "independent" && model != "BM" && model != "OU" && model != "Mk" )
    stop("specify a valid model option = 'autocorrelated', 'BM', 'independent' or 'Mk'")

  if(!(change.pr >= 0 & change.pr <= 1))
    stop("change.pr must be a probability between 0 and 1")

  if((model == "independent") & ( length(dist()) != 1 || !(is.numeric(dist()))))
    stop("specify a valid distribution function that returns a single value")

  if((model == "Mk") & ( init > k ))
    stop("Initial value incompatible with number of states k")

  if(is.null(taxonomy)){
    taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)
    from.taxonomy = FALSE
  } else
    from.taxonomy = TRUE

  aux = function(sp, t, r) {

    start = max(t$start[which(t$sp == sp)])
    end = min(t$end[which(t$sp == sp)])
    edges = t[which(t$sp == sp), ]

    blength = start - end

    # generate a new rate for sp
    if(model == "autocorrelated"){
      # this follows the molecular clock model of Kishino et al 2001
      # and the preservation model described in Heath et al 2014 (supplementary material)
      r = rlnorm(1, meanlog = log(r) - ((blength*v)/2), sdlog = sqrt(blength*v))
    } else if (model == "BM"){
      # regular Brownian motion
      r = rnorm(1, mean = r, sd = sqrt(blength * v))
    } else if (model == "OU") {
      mean = r * exp(-alpha)
      sd = sqrt(v/(2 * alpha) * (1 - exp(-2 * alpha)))
      r = rnorm(1, mean = mean, sd = sqrt(blength * v))
    } else if (model == "Mk") {
      # Pr of difference under symmetric Mk model
      if ( runif(1) < (( (k - 1) / k ) * (1 - exp(-( k / (k - 1) ) * v * blength) ) ) )
        r = sample(c(0:k-1)[-r], 1)
      # slow alternative for cross checking: x = rpois(1, blength * v); if(x > 0) for (i in 1:x) { r = sample(c(1:k)[-r], 1) }
    } else if (change.pr < 1) {
      if(runif(1) < change.pr)
        r = dist()
    } else { # independent rates
      r = dist()
    }

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

  taxonomy = aux(root, taxonomy, init)

  # extract unique rates
  rates = unique(cbind(taxonomy["sp"],taxonomy["rate"]))$rate

  if(!from.taxonomy) {
    rates = rates[order(as.numeric(taxonomy$sp))] # sort rates by node 1, node 2, etc
    if(!is.null(tree$root.edge)) rates = rates[c(root(tree), tree$edge[,2])] # sort rates according to tree
    else rates = rates[tree$edge[,2]]
  }

  return(rates)
}

#' Simulate an environmental gradient
#'
#' @description
#' Function returns a vector using the sine wave function \eqn{y = depth*sin(cycles*pi*(x-1/4))}
#' for a given set of intervals.
#' This vector can be used as a gradient to simulate fossils under an environment-dependent model of fossil recovery using the
#' function \code{sim.fossils.environment}.
#'
#' @param strata Number of stratigraphic intervals.
#' @param depth Maximum water depth.
#' @param cycles Number of cycles (e.g. transgressions and regressions).
#' @return vector of sampled water depths.
#' @examples
#' strata = 100
#' wd = sim.gradient(strata)
#' plot(wd, type="l")
#' @keywords non-uniform fossil preservation
#' @export
#' @seealso \code{\link{sim.fossils.environment}}
sim.gradient = function(strata, depth = 2, cycles = 2){

  # define the x-axis values
  x = seq(0,2,length.out=strata)

  # define y-axis values
  # a - total depth excursion - amplitude
  # b - number of cycles
  # 1/c - defines the relative start time of each cycle - phase shift
  # y = a * sin (b * pi * (x-1/c))
  y = depth*sin(cycles*pi*(x-1/4))

  return(y)
}
