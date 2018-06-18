#' sim.fbd.age: Simulating fossilized birth-death trees of a fixed age.
#'
#' @param age Time since origin.
#' @param numbsim Number of trees to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction rate.
#' @param psi Fossil sampling rate.
#' @param complete Keep unsampled extinct lineages?
#' @return Array of 'numbsim' trees with the time since origin being 'age'.
#' If tree goes extinct or no tips are sampled, return value is '0'.
#' If only one extant tip is sampled, return value is '1'.
#' @examples
#' age = 1
#' lambda = 2.0
#' mu = 0.5
#' psi = 0.6
#' numbsim = 2
#' sim.fbd.age(age, numbsim, lambda, mu, psi)
#' @keywords fossilized birth death
#' @export
sim.fbd.age<-function(age, numbsim, lambda, mu, psi, complete=FALSE)
{
	trees = TreeSim::sim.bd.age(age,numbsim,lambda,mu,complete=T)
	for(i in 1:length(trees))
	{
		if(!is.numeric(trees[[i]]))
		{
			t = trees[[i]]
			f <- sim.fossils.poisson(tree = t, rate = psi)

			f <- f[order(f$hmin, decreasing = TRUE),] # replaced f$h

			num_fossils = length(f$sp) # replaced f[,2]
			h <- numeric(num_fossils)
			fl <- character(num_fossils)

			tree = t
			origin = max(n.ages(tree)) + tree$root.edge
			if(num_fossils > 0)
			{
				for(j in 1:num_fossils)
				{
					node.ages = n.ages(tree)
					a = which(names(node.ages) == f$sp[j]) # replaced f[j,2]
					lineage.end = node.ages[[a]]

					h = f$hmin[j] - lineage.end # replaced f[j,1]

					# replaced f[j,2]
					tmp = phytools::bind.tip(tree, paste("fossil",j), edge.length = 0.0, where = f$sp[j], position = h)

					f$sp = map_nodes(f$sp, tree, tmp) # replaced f[,2]

					tree = tmp
				}
			}
			node.ages = n.ages(tree)
			if( complete == FALSE )
			{
				fossil.tips = tree$tip.label[which(node.ages[1:length(tree$tip.label)]>1e-7)]
				unsampled.tips = fossil.tips[!grepl("fossil",fossil.tips)]

				tree = ape::drop.tip(tree, unsampled.tips)
				node.ages = n.ages(tree)
			}
			trees[[i]] = tree
			trees[[i]]$root.edge = origin - max(node.ages)

			trees[[i]]$complete = complete
		}
	}
	return(trees)
}

#' sim.fbd.rateshift.taxa: Simulating fossilized birth death trees incorporating rate shifts.
#'
#' @param n Number of extant sampled tips.
#' @param numbsim Number of trees to simulate.
#' @param lambda Vector of speciation rates, the rate in entry i
#' is the speciation rate prior (ancestral) to time times[i].
#' @param mu Vector of extinction rates, the rate in entry i
#' is the extinction rate prior (ancestral) to time times[i].
#' @param psi Vector of fossil sampling rates, the rate in entry i
#' is the fossil sampling rate prior (ancestral) to time times[i].
#' @param times Vector of mass extinction and rate shift times.
#' Time is 0 today and increasing going backwards in time. Specify the vector as times[i]
#' @param complete Keep unsampled extinct lineages?
#' @return List of numbsim simulated trees with n extant sampled tips.
#' @examples
#' n = 10
#' numbsim = 1
#' sim.fbd.rateshift.taxa(n,numbsim,c(2,1),c(0,0.3),c(1,0.1),c(0.3))
#' @keywords fossilized birth death
#' @export
sim.fbd.rateshift.taxa <- function(n, numbsim, lambda, mu, psi, times, complete=FALSE)
{
	if(length(lambda) != (length(times) + 1 ))
    	stop("Length mismatch between interval ages and birth rates")
    if(length(mu) != (length(times) + 1 ))
    	stop("Length mismatch between interval ages and death rates")
	if(length(psi) != (length(times) + 1 ))
    	stop("Length mismatch between interval ages and sampling rates")

	trees = TreeSim::sim.rateshift.taxa(n, numbsim, lambda, mu, 1, times, complete = TRUE)

	for(i in 1:length(trees))
	{
		t = trees[[i]]
		origin = max(n.ages(t)) + t$root.edge

		horizons = c(0, times, origin)

		f <- sim.fossils.intervals(tree = t, interval.ages = horizons, rates = psi) # reordered

		f <- f[order(f$hmin, decreasing = TRUE),] # replaced f$h

		num_fossils = length(f$sp) # # replaced f[,2]
		h <- numeric(num_fossils)
		fl <- character(num_fossils)

		tree = t
		if(num_fossils > 0)
		{
			for(j in 1:num_fossils)
			{
				node.ages = n.ages(tree)

				a = which(names(node.ages) == f$sp[j]) # replaced f[j,2]
				lineage.end = node.ages[[a]]

				h = f$hmin[j] - lineage.end # replaced f[j,1]

				# replaced f[j,2]
				tmp = phytools::bind.tip(tree, paste("fossil",j), edge.length = 0.0, where = f$sp[j], position = h)

				f$sp = map_nodes(f$sp,tree,tmp) # replaced f[,2]

				tree = tmp
			}
		}
		node.ages = n.ages(tree)
		if( complete == FALSE )
		{
			fossil.tips = tree$tip.label[which(node.ages[1:length(tree$tip.label)]>1e-7)]
			unsampled.tips = fossil.tips[!grepl("fossil",fossil.tips)]

			tree = ape::drop.tip(tree, unsampled.tips)
			node.ages = n.ages(tree)
		}
		trees[[i]] = tree
		trees[[i]]$root.edge = origin - max(node.ages)

		trees[[i]]$complete = complete
	}
	return(trees)
}

#' sim.fbd.taxa: Simulating fossilized birth-death trees on a fixed number of extant taxa.
#'
#' @param n Number of extant sampled tips.
#' @param numbsim Number of trees to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction rate.
#' @param psi Fossil sampling rate.
#' @param complete Keep unsampled extinct lineages?
#' @return List of numbsim simulated trees with n extant sampled tips.
#' @examples
#' n = 10
#' lambda = 2.0
#' mu = 0.5
#' psi = 0.6
#' numbsim = 2
#' sim.fbd.taxa(n, numbsim, lambda, mu, psi)
#' @keywords fossilized birth death
#' @export
sim.fbd.taxa <- function(n, numbsim, lambda, mu, psi, complete=FALSE)
{
	trees = TreeSim::sim.bd.taxa(n, numbsim, lambda, mu, complete = TRUE)
	for(i in 1:length(trees))
	{
		t = trees[[i]]
		f <- sim.fossils.poisson(tree = t, rate = psi)

		f <- f[order(f$hmin, decreasing = T),] # replaced f$h

		num_fossils = length(f$sp) # replaced f[,2]
		h <- numeric(num_fossils)
		fl <- character(num_fossils)

		tree = t
		origin = max(n.ages(tree)) + tree$root.edge
		if(num_fossils > 0)
		{
			for(j in 1:num_fossils)
			{
				node.ages = n.ages(tree)

				a = which(names(node.ages) == f$sp[j]) # replaced f[j,2]
				lineage.end = node.ages[[a]]

				h = f$hmin[j] - lineage.end # replaced f[j,1]

				# replaced f[j,2]
				tmp = phytools::bind.tip(tree, paste("fossil",j), edge.length = 0.0, where = f$sp[j], position = h)

				f$sp = map_nodes(f$sp,tree,tmp) # replaced f[,2]

				tree = tmp
			}
		}
		node.ages = n.ages(tree)

		if( complete == FALSE )
		{
			fossil.tips = tree$tip.label[which(node.ages[1:length(tree$tip.label)]>1e-7)]
			unsampled.tips = fossil.tips[!grepl("fossil",fossil.tips)]

			tree = ape::drop.tip(tree, unsampled.tips)
			node.ages = n.ages(tree)
		}
		trees[[i]] = tree
		trees[[i]]$root.edge = origin - max(node.ages)

		trees[[i]]$complete = complete
	}
	return(trees)
}
