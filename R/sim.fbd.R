#' sim.cdfbd.taxa: Simulating Fossilized Birth Death Trees With Discrete Lineage-Specific Rate Regimes
#'
#' @param n Number of extant sampled tips.
#' @param numbsim Number of trees to simulate.
#' @param lambda Vector of lineage-specific speciation rates.
#' @param mu Vector of lineage-specific extinction rates.
#' @param psi Vector of lineage-specific fossil sampling rates.
#' @param rate Rate of switching between rate regimes.
#' @param pi Root state frequencies.
#' @param complete Keep unsampled extinct lineages?
#' @return List of numbsim simulated trees with n extant sampled tips.
#' @examples
#' n<-10
#' numbsim<-1
#' sim.cdfbd.taxa(n,numbsim,c(2,1),c(0,0.3),c(1,0.1),1)
#' @keywords fossilized birth death
#' @export
sim.cdfbd.taxa <- function(n,numbsim,lambda,mu,psi,rate,pi,complete=FALSE)
{
	k = length(psi)

	trees = TreeSim::sim.bd.taxa(n,numbsim,lambda,mu,complete=T)
	for(i in 1:length(trees))
	{
		t = trees[[i]]

		q = matrix(rep(rate,k*k),k,k,dimnames=list(LETTERS[1:k],LETTERS[1:k]))
		diag(q) <- rep(-rate,k)
		x = ape::rTraitDisc(t,"ER",k=k,rate=rate,root.value=sample.int(k,size=1,prob=pi))
		m = phytools::make.simmap(t,x,Q=q,message=FALSE)
		return(m)
		#simulate fossils
		f<-data.frame(h=numeric(),sp=numeric())

		root = length(t$tip.label) + 1
		lineages = c(t$edge[,2], root)

		node.ages = n.ages(t)

		for (j in lineages){ # internal nodes + tips
			if(j == root){

		      	# root age
		      	#a=which(names(node.ages)==root)
		      	#lineage.end=node.ages[[a]]

		      	# origin time
		      	#b=t$root.edge
		      	#lineage.start=lineage.end+b

		    } else {

		    	# work out the max age of the lineage (e.g. when that lineage became extant)
		      	# & get ancestor
		      	row=which(t$edge[,2]==j)
		      	ancestor=t$edge[,1][row]

		      	# get the age of the ancestor
		      	a=which(names(node.ages)==ancestor)
		      	lineage.start=node.ages[[a]]

		    	map = m$maps[[row]]

		    	for(kk in length(map))
		    	{
		    		# work out the min age of the lineage (e.g. when that lineage became extinct)
			      	# & get the branch length
			      	b=map[kk]
			      	lineage.end=lineage.start-b # branch length

			      	# sample fossil numbers from the Poisson distribution
			      	rate = psi[ which(LETTERS[1:k] == names(map)[kk]) ]
				    rand=rpois(1,b*rate)

				    if(rand > 0){
				      	h=runif(rand,min=lineage.end,max=lineage.start)
				      	f<-rbind(f,data.frame(h=h,sp=j))
				    }

				    lineage.start = lineage.end
				}
		    }	
		}


		f <- f[order(f$h,decreasing=T),]

		num_fossils = length(f[,2])
		h <- numeric(num_fossils)
		fl <- character(num_fossils)

		tree = t
		#origin = max(n.ages(tree)) + tree$root.edge
		if(num_fossils > 0)
		{
			for(j in 1:num_fossils)
			{
				node.ages = n.ages(tree)

				a=which(names(node.ages)==f[j,2])
				lineage.end = node.ages[[a]]

				h = f[j,1] - lineage.end
				
				tmp = phytools::bind.tip(tree,paste("fossil",j),edge.length=0.0,where=f[j,2],position=h)
				
				f[,2] = map_nodes(f[,2],tree,tmp)
				
				tree = tmp
			}
		}
		node.ages = n.ages(tree)
		if( complete == FALSE )
		{
			fossil.tips = tree$tip.label[which(node.ages[1:length(tree$tip.label)]>1e-7)]
			unsampled.tips = fossil.tips[!grepl("fossil",fossil.tips)]

			tree = ape::drop.tip(tree, unsampled.tips)
			#node.ages = n.ages(tree)
		}
		trees[[i]] = tree
		#trees[[i]]$root.edge = origin - max(node.ages)

		trees[[i]]$complete = complete
		class(trees[[i]]) <- c("phylo.fbd", class(trees[[i]]))
	}

	trees
}

#' sim.fbd.age: Simulating Fossilized Birth-Death Trees Of A Fixed Age.
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
#' age<-1
#' lambda <- 2.0
#' mu <- 0.5
#' psi <-0.6
#' numbsim<-2
#' sim.fbd.age(age,numbsim,lambda,mu,psi)
#' @keywords fossilized birth death
#' @export
sim.fbd.age<-function(age,numbsim,lambda,mu,psi,complete=FALSE)
{
	trees = TreeSim::sim.bd.age(age,numbsim,lambda,mu,complete=T)
	for(i in 1:length(trees))
	{
		if(!is.numeric(trees[[i]]))
		{
			t = trees[[i]]
			f <- sim.fossils.poisson(t, psi)
			f <- f[order(f$h,decreasing=T),]

			num_fossils = length(f[,2])
			h <- numeric(num_fossils)
			fl <- character(num_fossils)

			tree = t
			origin = max(n.ages(tree)) + tree$root.edge
			if(num_fossils > 0)
			{
				for(j in 1:num_fossils)
				{
					node.ages = n.ages(tree)

					a=which(names(node.ages)==f[j,2])
					lineage.end = node.ages[[a]]

					h = f[j,1] - lineage.end
					
					tmp = phytools::bind.tip(tree,paste("fossil",j),edge.length=0.0,where=f[j,2],position=h)
					
					f[,2] = map_nodes(f[,2],tree,tmp)
					
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
			class(trees[[i]]) <- c("phylo.fbd", class(trees[[i]]))
		}
	}

	trees
}

#' sim.fbd.rateshift.taxa: Simulating Fossilized Birth Death Trees Incorporating Rate Shifts.
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
#' n<-10
#' numbsim<-1
#' sim.rateshift.taxa(n,numbsim,c(2,1),c(0,0.3),c(1,0.1),c(0,0.3))
#' @keywords fossilized birth death
#' @export
sim.fbd.rateshift.taxa <- function(n,numbsim,lambda,mu,psi,times,complete=FALSE)
{
	if(length(psi) != (length(times) + 1 ))
    	stop("Length mismatch between interval ages and sampling rates")

	trees = TreeSim::sim.rateshift.taxa(n,numbsim,lambda,mu,1,times,complete=T)

	for(i in 1:length(trees))
	{
		t = trees[[i]]
		origin = max(n.ages(t)) + t$root.edge

		horizons = c(0, times, origin)

		f <- sim.fossils.non.unif(t, horizons, psi)
		f <- f[order(f$h,decreasing=T),]

		num_fossils = length(f[,2])
		h <- numeric(num_fossils)
		fl <- character(num_fossils)

		tree = t
		if(num_fossils > 0)
		{
			for(j in 1:num_fossils)
			{
				node.ages = n.ages(tree)

				a=which(names(node.ages)==f[j,2])
				lineage.end = node.ages[[a]]

				h = f[j,1] - lineage.end
				
				tmp = phytools::bind.tip(tree,paste("fossil",j),edge.length=0.0,where=f[j,2],position=h)
				
				f[,2] = map_nodes(f[,2],tree,tmp)
				
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
		class(trees[[i]]) <- c("phylo.fbd", class(trees[[i]]))
	}

	trees
}

#' sim.fbd.taxa: Simulating Fossilized Birth-Death Trees On A Fixed Number Of Extant Taxa.
#'
#' @param n Number of extant sampled tips.
#' @param numbsim Number of trees to simulate.
#' @param lambda Speciation rate.
#' @param mu Extinction rate.
#' @param psi Fossil sampling rate.
#' @param complete Keep unsampled extinct lineages?
#' @return List of numbsim simulated trees with n extant sampled tips.
#' @examples
#' n<-10
#' lambda <- 2.0
#' mu <- 0.5
#' psi <-0.6
#' numbsim<-2
#' sim.fbd.taxa(n,numbsim,lambda,mu,psi)
#' @keywords fossilized birth death
#' @export
sim.fbd.taxa <- function(n,numbsim,lambda,mu,psi,complete=FALSE)
{
	trees = TreeSim::sim.bd.taxa(n,numbsim,lambda,mu,complete=T)
	for(i in 1:length(trees))
	{
		t = trees[[i]]
		f <- sim.fossils.poisson(t, psi)
		f <- f[order(f$h,decreasing=T),]

		num_fossils = length(f[,2])
		h <- numeric(num_fossils)
		fl <- character(num_fossils)

		tree = t
		origin = max(n.ages(tree)) + tree$root.edge
		if(num_fossils > 0)
		{
			for(j in 1:num_fossils)
			{
				node.ages = n.ages(tree)

				a=which(names(node.ages)==f[j,2])
				lineage.end = node.ages[[a]]

				h = f[j,1] - lineage.end
				
				tmp = phytools::bind.tip(tree,paste("fossil",j),edge.length=0.0,where=f[j,2],position=h)
				
				f[,2] = map_nodes(f[,2],tree,tmp)
				
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
		class(trees[[i]]) <- c("phylo.fbd", class(trees[[i]]))
	}

	trees
}