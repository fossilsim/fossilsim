#' sim.fbd.age: Simulating fossilized birth-death trees of a fixed age.
#'
#' @param age Time since origin / most recent common ancestor.
#' @param numbsim Number of trees to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction rate.
#' @param psi Fossil sampling rate.
#' @param frac Sampling fraction: The actual number of tips is n/frac, but
#' only n tips are included (incomplete sampling).
#' @param mrca If mrca=FALSE: age is the time since origin. If mrca=TRUE:
#' age is the time since most recent common ancestor of the
#' extant tips.
#' @param complete If complete = TRUE, the tree with the extinct and non-sampled
#' lineages is returned. If complete = FALSE, the extinct and
#' non-sampled lineages are suppressed.
#' @param K If K=0, then lambda is constant. If K>0, density-dependent
#' speciation is assumed, with speciation rate = lambda(1-m/K)
#' when there are m extant species.
#' @return Array of 'numbsim' SAtrees with the time since origin / most
#' recent common ancestor being 'age'. If tree goes extinct or
#' no tips are sampled (only possible when mrca = FALSE), return
#' value is '0'. If only one extant and no extinct tips are
#' sampled, return value is '1'.
#' @examples
#' age = 1
#' lambda = 2.0
#' mu = 0.5
#' psi = 0.6
#' numbsim = 2
#' sim.fbd.age(age, numbsim, lambda, mu, psi)
#' @keywords fossilized birth death
#' @export
sim.fbd.age<-function(age, numbsim, lambda, mu, psi, frac = 1, mrca = FALSE, complete = FALSE, K = 0)
{
	trees = TreeSim::sim.bd.age(age,numbsim,lambda,mu,frac,mrca,complete=T,K)
	for(i in 1:length(trees))
	{
		if(!is.numeric(trees[[i]]))
		{
			t = trees[[i]]
			f <- sim.fossils.poisson(tree = t, rate = psi)
			
			tree = SAtree.from.fossils(t,f)

			node.ages = n.ages(tree)

			origin = max(n.ages(tree)) + tree$root.edge

			if( complete == FALSE )
			{
				fossil.tips = is.extinct(tree,tol=0.000001)
				extant.tips = tree$tip.label[!(tree$tip.label %in% fossil.tips)]

				sa.tips = tree$tip.label[tree$edge[,2][(tree$edge[,2] %in% 1:length(tree$tip.label)) & (tree$edge.length == 0.0)]]

				unsampled.tips = fossil.tips[!(fossil.tips %in% sa.tips)]

				if( frac < 1 )
				{
				  	n = round(length(extant.tips) * frac)
					unsampled.tips <- c( unsampled.tips, extant.tips[!(extant.tips %in% sample(extant.tips, n))] )
				}

				tree = ape::drop.tip(tree, unsampled.tips)
				node.ages = n.ages(tree)
			}
			trees[[i]] = tree
			if( mrca == FALSE )
			{
				trees[[i]]$root.edge = origin - max(node.ages)
			}

			trees[[i]] = SAtree(trees[[i]],complete)
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
#' @param complete If TRUE, the tree including the extinct lineages and
#' non-sampled lineages is returned (so the tree has
#' round(n/frac[1]) extant tips). If FALSE, the extinct lineages
#' and non-sampled lineages are suppressed.
#' @return List of numbsim simulated SAtrees with n extant sampled tips.
#' @examples
#' n = 10
#' numbsim = 1
#' sim.fbd.rateshift.taxa(n,numbsim,c(2,1),c(0,0.3),c(1,0.1),c(0.3))
#' @keywords fossilized birth death
#' @export
sim.fbd.rateshift.taxa <- function(n, numbsim, lambda, mu, psi, times, complete = FALSE)
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

		tree = SAtree.from.fossils(t,f)

		node.ages = n.ages(tree)
		
		if( complete == FALSE )
		{
			fossil.tips = is.extinct(tree,tol=0.000001)

			sa.tips = tree$tip.label[tree$edge[,2][(tree$edge[,2] %in% 1:length(tree$tip.label)) & (tree$edge.length == 0.0)]]

			unsampled.tips = fossil.tips[!(fossil.tips %in% sa.tips)]

			tree = ape::drop.tip(tree, unsampled.tips)
			node.ages = n.ages(tree)
		}
		trees[[i]] = tree
		trees[[i]]$root.edge = origin - max(node.ages)

		trees[[i]] = SAtree(trees[[i]],complete)
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
#' @param frac When complete = FALSE: The actual
#' number of extant tips is n/frac, but only n tips are included
#' (incomplete sampling). When complete = TRUE: all extinct and
#' non-sampled lineages are included, i.e. the tree has n/frac
#' extant tips.
#' @param complete If TRUE, the tree with the extinct and non-sampled lineages
#' is returned. If FALSE, the extinct lineages are suppressed.
#' @return List of numbsim simulated SAtrees with n extant sampled tips.
#' @examples
#' n = 10
#' lambda = 2.0
#' mu = 0.5
#' psi = 0.6
#' numbsim = 2
#' sim.fbd.taxa(n, numbsim, lambda, mu, psi)
#' @keywords fossilized birth death
#' @export
sim.fbd.taxa <- function(n, numbsim, lambda, mu, psi, frac = 1, complete = FALSE)
{
	trees = TreeSim::sim.bd.taxa(n, numbsim, lambda, mu, frac, complete = TRUE)
	for(i in 1:length(trees))
	{
		t = trees[[i]]
		f <- sim.fossils.poisson(tree = t, rate = psi)

		tree = SAtree.from.fossils(t,f)

		node.ages = n.ages(tree)

		origin = max(node.ages) + tree$root.edge

		if( complete == FALSE )
		{
			fossil.tips = is.extinct(tree,tol=0.000001)
			extant.tips = tree$tip.label[!(tree$tip.label %in% fossil.tips)]

			sa.tips = tree$tip.label[tree$edge[,2][(tree$edge[,2] %in% 1:length(tree$tip.label)) & (tree$edge.length == 0.0)]]

			unsampled.tips = fossil.tips[!(fossil.tips %in% sa.tips)]

			if( frac < 1 )
			{
				unsampled.tips <- c( unsampled.tips, extant.tips[!(extant.tips %in% sample(extant.tips, n))] )
			}

			tree = ape::drop.tip(tree, unsampled.tips)
			node.ages = n.ages(tree)
		}
		trees[[i]] = tree
		trees[[i]]$root.edge = origin - max(node.ages)

		trees[[i]] = SAtree(trees[[i]],complete)
	}
	return(trees)
}
