#' Generate input file for analysis using the program pyrate
#'
#' @param fossils Fossils object.
#' @param python If TRUE the function outputs a file in python format.
#' If FALSE the function outputs a standard text file in the format used by pyrate (default).
#' @param traits Vector of trait values equal to the number of specimens, including extant samples if present in the fossils dataframe.
#' Each entry corresponds to the order in which they appear in fossils$sp.
#' @param cutoff Exclude occurrences with age uncertainty greater than this value. i.e. \code{hmax - hmin > cutoff}.
#' @param random If TRUE use a random number from within the interval U(hmin, hmax) for specimen ages,
#' otherwise use the midpoint of this interval (default). Applicable only when \code{python = TRUE} and for specimens with \code{hmin != hmax}.
#' @param min Value used to represent the minimum possible interval age of extinct specimens with \code{hmin = 0}. By default \code{min = NULL} and the function uses the sampling times in the fossils dataframe.
#' @param print.extant If TRUE also output sampling times and trait values for extant samples. Default = FALSE.
#' @param file Output file name.
#'
#' @examples
#'
#' set.seed(123)
#'
#' # simulate tree
#' t = ape::rtree(6)
#'
#' # assign a max age based on tree height
#' max.age = tree.max(t)
#'
#' # define a set of non-uniform length intervals
#' times = c(0, sort(runif(3, min = 0, max = max.age)), max.age)
#' rates = c(1,2,3,4)
#'
#' # simulate fossils reflect age uncertainty
#' f = sim.fossils.intervals(tree = t, interval.ages = times, rates = rates,
#'     use.exact.times = FALSE)
#'
#' # simulate extant samples
#' rho = 1
#' f = sim.extant.samples(f, t, rho = 1)
#'
#' plot(f, t)
#'
#' # generate input files for pyrate
#' pyrate.input(f)
#' pyrate.input(f, python = TRUE)
#'
#' # add trait values
#' traits = runif(length(f$sp))
#' pyrate.input(f, traits = traits)
#'
#' @export
pyrate.input = function(fossils, python = FALSE, traits = NULL, cutoff = NULL, random = FALSE, min = NULL, print.extant = FALSE, file = ""){

  if(length(fossils$sp) < 1) stop("Number of specimens must be > 0")

  if(!is.null(traits) && length(fossils$sp) != length(traits))
    stop("Number of trait values must equal the number of fossils")

  if(!is.null(traits)) fossils = round(cbind(fossils, traits), digits = 6)
  else fossils = round(fossils, digits = 6)

  if(!is.null(cutoff) && !(cutoff > 0)) stop("Cutoff must be > 0")

  # tolerance for extant specimens age bounds = 0
  tol = 1e-8

  if(!is.null(min)){
    fossils$hmin[which(fossils$hmin != fossils$hmax & fossils$hmin < tol)] = min
    if(any(fossils$hmin > fossils$hmax)) stop("min value has generated intervals with hmax < hmin")
  }

  # deal with specimen age options
  if(!is.null(cutoff) && cutoff > 0){
    cat(paste0("Excluding ", length(which(fossils$hmax - fossils$hmin > cutoff)), " occurrences based on cutoff...\n"))
    fossils = fossils[-which(fossils$hmax - fossils$hmin > cutoff),]
    # this logic is incorrect
    if(length(fossils$sp) < 1) stop("Number of specimens after cutoff must be > 0")
    else if(!print.extant && length(fossils[-which(fossils$hmax - fossils$hmin < tol),]$sp) < 1) stop("Number of specimens after cutoff must be > 0")
  }

  if(python){
    # use random ages
    if(random && any(fossils$hmax != fossils$hmin)){
      fossils$hmin = round(unlist(lapply(1:dim(fossils)[1], function(x) { runif(1, fossils$hmin[x], fossils$hmax[x]) })), digits = 6)
      fossils$hmax = fossils$hmin
    } else if (any(fossils$hmax != fossils$hmin)) { # use median ages
      fossils$hmin = round(unlist(lapply(1:dim(fossils)[1], function(x) { mean( c(fossils$hmin[x], fossils$hmax[x])) })), digits = 6)
      fossils$hmax = fossils$hmin
    }
  }

  total = length(unique(fossils$sp))

  if(python){
    cat("#!/usr/bin/env python", "from numpy import * ", "",  sep = "\n", file = file, append = FALSE)

    cat("data_1=[", sep = "", file = file, append = TRUE)

    data_1 = c()
    trait1 = c()
    taxa_names = c()
    for(i in 1:total){
      sp = unique(fossils$sp)[i]
      occs = fossils[which(fossils$sp == sp),]
      times = c()
      values = c()
      for(j in 1:length(occs$sp)){
        # skip extant samples
        if(occs$hmin[j] < tol && !print.extant) next
        else times = c(times, occs$hmin[j])
        if(!is.null(traits)) values = c(values, occs$traits[j])
      }

      if(length(times) > 0){
        data_1 = c(data_1, paste0("array([", paste(times, collapse = ","), "])"))
        taxa_names = c(taxa_names, paste0("'taxa", sp, "'"))
        if(!is.null(traits)) trait1 = c(trait1, paste0("array([", paste(values, collapse = ","), "])"))
      }
    }

    # print fossil ages
    cat(paste(data_1, collapse = ",\n"), "]\n", sep = "\n", file = file, append = TRUE)

    cat("d=[data_1]", "names=['example_1']", "def get_data(i): return d[i]", "def get_out_name(i): return  names[i]",
        sep = "\n", file = file, append = TRUE)

    # print species names
    cat("taxa_names=[", paste(taxa_names, collapse = ","), "]\n", sep = "", file = file, append = TRUE)

    cat("def get_taxa_names(): return taxa_names", sep = "\n", file = file, append = TRUE)

    # print trait values
    if(!is.null(traits)){
      cat("\ntrait1=[", paste(trait1, collapse = ",\n"), "\n]\n\ntraits=[trait1]\ndef get_continuous(i): return traits[i]\n", sep = "", file = file, append = TRUE)
    }

  } else {

    if(is.null(traits))
      cat("Species\tStatus\tmin_age\tmax_age\n", sep = "", file = file, append = FALSE)
    else
      cat("Species\tStatus\tmin_age\tmax_age\ttraits\n", file = file, append = FALSE)

    for(i in 1:total){
      status = NULL
      sp = unique(fossils$sp)[i]
      occs = fossils[which(fossils$sp == sp),]
      # identify extant lineages
      if(any(occs$hmin == occs$hmax & occs$hmin < tol))
        status = "extant"
      else status = "extinct"

      for(j in 1:length(occs$sp)){
        # skip extant samples
        if(occs$hmin[j] == occs$hmax[j] && occs$hmin[j] < tol && !print.extant) next
        if(is.null(traits))
          cat(paste0("taxa",i), "\t", status, "\t", occs$hmin[j], "\t", occs$hmax[j], "\n", file = file, append = TRUE)
        else
          cat(paste0("taxa",i), "\t", status, "\t", occs$hmin[j], "\t", occs$hmax[j], "\t", occs$traits[j], "\n", file = file, append = TRUE)
      }
    }
  }

}
