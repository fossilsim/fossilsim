#' Reassign fossil ages to user-specified stratigraphic intervals
#'
#' Reassign exact fossil ages using the minimum and maximum ages of a set of stratigraphic intervals.
#' If \code{use.species.ages = TRUE} the function will respect species durations and will not
#' return minimum and maximum bounds that may be younger or older than the species durations.
#' This requires supplying a phylo or taxonomy object.
#'
#' @param fossils Fossil object.
#' @param tree Phylo object.
#' @param taxonomy Taxonomy object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param max.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param use.species.ages If TRUE reassigned fossil ages will respect the speciation times. Default = FALSE.
#' @param root.edge If TRUE include root edge.
#' @param sim.extant If TRUE simulate age uncertainty for extant samples as well, default FALSE.
#'
#' @return An object of class fossils.
#'
#' @examples
#' # simulate tree
#' t = ape::rtree(8)
#'
#' # simulate fossils
#' rate = 2
#' f = sim.fossils.poisson(rate, t)
#' plot(f, t)
#'
#' # assign a max age based on tree height
#' max.age = tree.max(t)
#'
#' # define intervals
#' times = seq(0, max.age, length.out = 5)
#'
#' # reassign ages
#' f = sim.interval.ages(f, t, interval.ages = times)
#'
#' # plot output
#' plot(f, t, interval.ages = times)
#'
#' @export
sim.interval.ages = function(fossils, tree = NULL, taxonomy = NULL,
                        interval.ages = NULL, max.age = NULL, strata = NULL,
                        use.species.ages = FALSE, root.edge = TRUE, sim.extant = FALSE){

  if(is.null(fossils))
    stop("Specify fossils object")

  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(use.species.ages){

    if(is.null(tree) && is.null(taxonomy))
      stop("Specify phylo or taxonomy object to use species ages")

    if(!is.null(tree) && !"phylo" %in% class(tree))
      stop("tree must be an object of class \"phylo\"")

    if(!is.null(taxonomy) && !"taxonomy" %in% class(taxonomy))
      stop("taxonomy must be an object of class \"taxonomy\"")

    if(!is.null(tree) && !is.null(taxonomy))
      warning("tree and taxonomy both defined, using taxonomy")

    if(is.null(taxonomy)){
      taxonomy = sim.taxonomy(tree, beta = 1, root.edge = root.edge)
      from.taxonomy = FALSE
    } else from.taxonomy = TRUE

    if(any(!(fossils$sp %in% taxonomy$sp)))
      stop("incompatible fossil species and taxonomy")
  }

  # check fossil ages are not already binned
  if(!identical(fossils$hmin, fossils$hmax))
    stop("exact fossil sampling times must be specified to use this function (i.e. hmin = hmax)")

  if(is.null(interval.ages) && (is.null(max.age) || is.null(strata)))
    stop("Intervals need to be defined by specifying either interval.ages or max.age and strata")
  if(!is.null(max.age) && !is.null(strata)) {
    if(!is.null(interval.ages)) warning("Two interval definitions found, using interval.ages")
    else interval.ages <- seq(0, max.age, length = strata + 1)
  }

  if(any(fossils$hmin > max(interval.ages)))
    stop("some fossils are older than the oldest bin")
  if(any(fossils$hmin < min(interval.ages)))
    stop("some fossils are younger than the youngest bin")

  # for each fossil
  for(i in 1:length(fossils$hmin)){
    if(!sim.extant && fossils$hmin[i] < 1e-8) next

    int = assign.interval(interval.ages, fossils$hmin[i])

    if(use.species.ages){ # { assign hmin and hmax that do not violate species start and end times }

      sp = fossils$sp[i]
      start = max(taxonomy$start[which(taxonomy$sp == sp)])
      end = min(taxonomy$end[which(taxonomy$sp == sp)])
      start.int = assign.interval(interval.ages, start)
      end.int = assign.interval(interval.ages, end)

      # this is incorrect
      if(start.int == int || end.int == int){
        if(start.int == int & end.int == int){ # speciation and extinction occur within the bin }
          fossils$hmin[i] = end
          fossils$hmax[i] = start
        } else if (start.int == int) { # { speciation occurs within the bin }
          fossils$hmin[i] = interval.ages[int]
          fossils$hmax[i] = start
        } else if (end.int == int) { # { exinction occurs within the bin }
          fossils$hmin[i] = end
          fossils$hmax[i] = interval.ages[int+1]
        }
      } else {
        fossils$hmin[i] = interval.ages[int]
        fossils$hmax[i] = interval.ages[int+1]
      }

    } else { # { else ignore species ages }
      fossils$hmin[i] = interval.ages[int]
      fossils$hmax[i] = interval.ages[int+1]
    }
  }
  if(use.species.ages)
    fossils <- as.fossils(fossils, from.taxonomy)
  return(fossils)
}


# Calculates minimum and maximum ages of fossils, based on stratigraphic intervals.
# Allows a fraction of error on the correct stratigraphic interval.
#
# @param fossils Fossils object.
# @param intervals Stratigraphic intervals as a vector of interval boundaries in increasing order.
# @param error_fraction Probability of error in interval identification per fossil. Default 0.
# @param error_range Maximum deviation from the correct interval, expressed in number of intervals. Default 1.
# @param sd_intervals Standard deviations on interval boundaries. If null, all deviations are assumed 0.
# @param save_params Whether simulation parameters (error_fraction and error_range) should be saved as attributes in the returned object.
# @return Fossils object with min and max ages of fossils.
# @examples
# # simulate tree
# t = ape::rtree(6)
#
# # simulate fossils
# f = sim.fossils.poisson(rate = 2, tree = t)
#
# # intervals
# max.age = tree.max(t)
# times = seq(0, max.age, length.out = 4)
#
# # simulate uncertainty
# f2 = sim.uncertain.ages(f, times)
# plot(f2, t, binned = TRUE, interval.ages = times)
sim.uncertain.ages = function(fossils, intervals, error_fraction = 0, error_range = 1, sd_intervals = NULL, save_params = F) {
  if(!is.null(sd_intervals) && length(intervals) != length(sd_intervals)) {
    stop("Number of intervals do not match number of interval deviations")
  }
  if(error_fraction > 0 && error_range == 0) {
    warning("Error range is 0, no errors will be added")
  }

  fossils$h = (fossils$hmin + fossils$hmax)/2
  idx = findInterval(fossils$h,intervals)
  if(any(idx == 0) || any(idx == length(intervals))) {
    stop("All fossils must be inside of the specified intervals")
  }

  #integrating errors
  if(error_fraction >0 && error_range >0) {
    for(i in 1:length(idx)) {
      if(runif(1) < error_fraction) {
        repeat {
          #add a random error in range of error_range
          error_id = idx[i] + sample(c(1:error_range, -1:-error_range),1)
          #make sure the id+error is still in the intervals
          if(error_id >0 && error_id < length(intervals)) break
        }
        idx[i] = error_id
      }
    }
  }

  if(is.null(sd_intervals)) sd_intervals = rep(0, length(intervals))
  #hmin and hmax from intervals with maximum uncertainty
  fossils$hmin = intervals[idx] - sd_intervals[idx]
  fossils$hmax = intervals[idx + 1] + sd_intervals[idx + 1]

  #optional: save simulation parameters
  if(save_params) {
    attr(fossils, "error_fraction") = error_fraction
    attr(fossils, "error_range") = error_range
  }

  fossils
}

# assign any given age to one of a set of intervals
assign.interval = function(intervals, t) {

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
