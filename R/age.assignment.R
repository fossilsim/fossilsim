#' Reassign fossil ages
#'
#' Reassign exact fossil ages using the minimum and maximum ages of a set of stratigraphic intervals.
#' If \code{use.species.ages = TRUE} the function will respect species durations and will not
#' return minimum and maximum bounds that may be younger or older than the species durations.
#' This requires supplying a phylo or taxonomy object.
#'
#' @param fossils Fossil object.
#' @param tree Phylo object.
#' @param species Taxonomy object.
#' @param interval.ages Vector of stratigraphic interval ages, starting with the minimum age of the youngest interval and ending with the maximum age of the oldest interval.
#' @param basin.age Maximum age of the oldest stratigraphic interval.
#' @param strata Number of stratigraphic intervals.
#' @param use.species.ages If TRUE reassigned fossil ages will respect the speciation times. Default = FALSE.
#' @param root.edge If TRUE include root edge.
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
#' max.age = basin.age(t)
#'
#' # define intervals
#' times = seq(0, max.age, length.out = 5)
#'
#' # reassign ages
#' f = reassign.ages(f, t, interval.ages = times)
#'
#' @export
reassign.ages<-function(fossils, tree = NULL, species = NULL,
                        interval.ages = NULL, basin.age = NULL, strata = NULL,
                        use.species.ages = FALSE, root.edge = TRUE){

  if(is.null(fossils))
    stop("Specify fossils object")

  if(!is.null(fossils) && !"fossils" %in% class(fossils))
    stop("fossils must be an object of class \"fossils\"")

  if(use.species.ages){

    if(is.null(tree) && is.null(species))
      stop("Specify phylo or taxonomy object to use species ages")

    if(!is.null(tree) && !"phylo" %in% class(tree))
      stop("tree must be an object of class \"phylo\"")

    if(!is.null(species) && !"taxonomy" %in% class(species))
      stop("species must be an object of class \"taxonomy\"")

    if(!is.null(tree) && !is.null(species))
      warning("tree and species both defined, using species taxonomy")

    #if(!is.null(tree) && is.null(species))
    #  warning("using tree, assumming all speciation is symmetric")

    if(is.null(species)){
      species = create.taxonomy(tree, beta = 1, root.edge = root.edge)
      from.taxonomy = FALSE
    } else from.taxonomy = TRUE

    if(any(!(fossils$sp %in% species$sp)))
      stop("incompatible fossil species and taxonomy")
  }

  # check fossil ages are not already binned
  if(!identical(fossils$hmin, fossils$hmax))
    stop("exact fossil sampling times must be specified to use this function (i.e. hmin = hmax)")

  if(is.null(interval.ages) && (is.null(basin.age) || is.null(strata)))
    stop("Intervals need to be defined by specifying either interval.ages or basin.age and strata")
  if(!is.null(basin.age) && !is.null(strata)) {
    if(!is.null(interval.ages)) warning("Two interval definitions found, using interval.ages")
    else interval.ages <- seq(0, basin.age, length = strata + 1)
  }

  if(any(fossils$hmin > max(interval.ages)))
    stop("some fossils are older than the oldest bin")
  if(any(fossils$hmin < min(interval.ages)))
    stop("some fossils are younger than the youngest bin")

  # for each fossil
  for(i in 1:length(fossils$hmin)){

    int = assign.interval(interval.ages, fossils$hmin[i])

    if(use.species.ages){ # { assign hmin and hmax that do not violate species start and end times }

      sp = fossils$sp[i]
      start = species$start[which(species$sp == sp)][1]
      end = species$end[which(species$sp == sp)][1]
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

