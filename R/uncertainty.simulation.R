#' Calculates minimum and maximum ages of fossils, based on stratigraphic intervals.
#' Allows a fraction of error on the correct stratigraphic interval.
#'
#' @param fossils Fossils object.
#' @param intervals Stratigraphic intervals as a vector of interval boundaries in increasing order.
#' @param error_fraction Probability of error in interval identification per fossil. Default 0.
#' @param error_range Maximum deviation from the correct interval, expressed in number of intervals. Default 1.
#' @param sd_intervals Standard deviations on interval boundaries. If null, all deviations are assumed 0.
#' @param save_params Whether simulation parameters (error_fraction and error_range) should be saved as attributes in the returned object.
#' @return Fossils object with min and max ages of fossils.
#' @examples
#' # simulate tree
#' t <- ape::rtree(6)
#' # simulate fossils
#' f <- sim.fossils.poisson(rate = 2, tree = t)
#' # simulate uncertainty
#' f = sim.uncertain.ages(f, c(0,TreeSim::getx(t)))
#' @export
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
