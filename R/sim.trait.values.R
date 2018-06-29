#' Simulate trait values with variation across lineages
#'
#' @description
#' Fossil recovery rates or other parameter values can be simulated for a phylo (\code{tree}) or taxonomy (\code{taxonomy}) object.
#' Available models include "autocorrelated", "independent" and a "innovative" model in which changes in parameter values are coincident with speciation events.
#' Under the \code{autocorrelated} model, trait values evolve along lineages according to a Brownian motion process, where the strength of the relationship between ancestor and descendant values is determined by the parameter \eqn{\nu} (\code{v}).
#' If \eqn{\nu} is small values will be more similar between ancestor and descendants, and if \eqn{\nu} is zero all trait values will be equal.
#' For a given species \eqn{i} with ancestor \eqn{j}, a new trait value \eqn{\kappa_i} is drawn from a lognormal distribution with
#' \deqn{\kappa_i ~ LN( ln([\kappa_j] - (\sigma^2/2), \sigma)}
#' where \eqn{\sigma = \nu * t_i} and \eqn{t_i} is the lineage duration of the species.
#' This fossil recovery model is described in Heath et al. (2014) and is equivalent to the autocorrelated relaxed clock model described in Kishino et al. (2001).
#' Under the \code{independent} model a new trait value is drawn for each species from any valid user-specified distribution (\code{dist}).
#' Under the \code{innovative} model values change at each speciation event with a given probability (\code{innovative.pr}) and new values are drawn from any valid user-specified distribution (\code{dist}).
#'
#' @param rate Initial value at the origin or root of the phylo or taxonomy object. Default = 1.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param root.edge If TRUE include the root edge. Default = TRUE.
#' @param model Model used to simulate rate variation across lineages. Options include "autocorrelated" (default), "independent" or "innovative".
#' @param v Brownian motion parameter \eqn{v} used in the autocorrelated model. Default = 0.01.
#' @param dist Distribution of trait values used to draw new values under the "independent" and "innovative" models. This parameter is ignored if \code{model = "autocorrealted"}. The default is a uniform distribution with \emph{U(0, 2)}. The distribution function must return a single value.
#' @param change.pr Probability that trait values change at speciation events. Default = 0.01.
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
#' # simulate rates under the autocorrelated rates model
#' rate = 1
#' rates = sim.trait.values(rate = rate, taxonomy = s, v = 1)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, t)
#'
#' # simulate rates under the independent rates model
#' dist = function() { rlnorm(1, log(rate), 1) }
#' rates = sim.trait.values(rate = rate, taxonomy = s, model = "independent", dist = dist)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, t)
#'
#' # simualte rates under the innovative model
#' rates = sim.trait.values(rate = rate, taxonomy = s, model = "innovative",
#'                         dist = dist, change.pr = 0.1)
#' f = sim.fossils.poisson(rates, taxonomy = s)
#' plot(f, t)
#'
#' @references
#' Heath et al. 2014. The fossilized birth-death process for coherent calibration of divergence-time estimates. PNAS 111:E2957-E2966.\cr
#' Kishino et al. 2001. Performance of a divergence time estimation method under a probabilistic model of rate evolution MBE 18:352-361.
#'
#' @export
sim.trait.values = function(rate = 1, tree = NULL, taxonomy = NULL, root.edge = TRUE,
                             model = "autocorrelated", v = 0.01,
                             dist = function(){runif(1,0,2)}, change.pr = 0.01){

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

  if(model != "autocorrelated" && model != "independent" && model != "innovative")
    stop("specify a valid model option = 'autocorrelated', 'independent' or 'innovative'")

  if(!(change.pr >= 0 & change.pr <= 1))
    stop("change.pr must be a probability between 0 and 1")

  if((model == "independent" || model == "innovative") & ( length(dist()) != 1 || !(is.numeric(dist()))))
    stop("specify a valid distribution function that returns a single value")

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
    } else if (model == "innovative") {
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

  taxonomy = aux(root, taxonomy, rate)

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
