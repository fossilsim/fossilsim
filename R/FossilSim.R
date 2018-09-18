#' FossilSim: A package for simulating and plotting fossil and taxonomy data
#'
#' This package provides functions for simulating both taxonomy and fossil data from an existing phylogeny.
#' 
#' @section Simulating taxonomy:
#' Taxonomy can be simulated in FossilSim under a mixed model of speciation that can incorporate three modes of speciation -- 
#' budding (or asymmetric), bifurcating (or symmetric) and anagenetic -- in addition to cryptic speciation. 
#' A description of the resulting taxonomy objects and simulation functions can be found in the "Simulating taxonomy" vignette.
#' 
#' @section Simulating fossil data:
#' Fossils can be simulated from a phylogeny or a taxonomy under a model of constant fossilization rate or 
#' time-dependent, environment-dependent or species-dependent fossilization rates.
#' A description of the resulting fossil objects and simulation functions can be found in the "Simulating fossils" vignette.
#' 
#' @section Plotting functions:
#' Both taxonomy and fossil objects are provided with custom plotting functions that highlight important features of the simulated objects 
#' along the original phylogeny. More details about these functions can be found in the vignettes or by calling 
#' \code{?plot.taxonomy} and \code{?plot.fossils}.
#' 
#' @section Compatibility with other packages:
#' FossilSim is designed to use phylogenies in the ape format. It provides functions to convert to and from the fossilRecordSimulation 
#' format used by the package paleotree (see the vignette "Converting from and to paleotree format"), as well as functions to convert to the 
#' zero-edge format used for instance by BEAST2 and RevBayes (see the vignette "Exporting sampled ancestor trees").
#' 
#' @examples 
#' # simulate a tree using TreeSim conditioned on tip number
#' t = TreeSim::sim.bd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0.2)[[1]]
#' 
#' # simulate taxonomy under mixed speciation
#' s = sim.taxonomy(tree = t, beta = 0.5, lambda.a = 1, kappa = 0.1)
#' # plot the result
#' plot(s, tree = t, legend.position = "topleft")
#' 
#' # simulate fossils using the phylogeny and a constant fossilization rate
#' f = sim.fossils.poisson(rate = 3, tree = t)
#' # plot the result
#' plot(f, tree = t)
#' 
#' # simulate fossils using the taxonomy and a constant fossilization rate
#' f = sim.fossils.poisson(rate = 3, taxonomy = s)
#' # plot the result
#' plot(f, tree = t, taxonomy = s, show.taxonomy = TRUE)
#' 
#' # simulate fossils using time-dependent fossilization rates
#' f = sim.fossils.intervals(tree = t, rates = c(1, 0.1, 1, 0.1), max.age = tree.max(t), strata = 4)
#' # plot the result, with the time intervals
#' plot(f, tree = t, show.strata = TRUE, max.age = tree.max(t), strata = 4)
#'
#' @docType package
#' @name FossilSim
NULL