library(devtools)
install_github("rachelwarnock/fossilsim")
# its a good idea to restart R after installing/updating packages via github or bitbucket

set.seed(123)

## simulate birth death tree
t <- TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]

## simulate fossils under a Poisson sampling process
# Poisson sampling rate
sampling = 3
f <- FossilSim::sim.fossils.poisson(t, sampling)
plot(f, t)

## simulate fossils under a unifom sampling model for a set of equal length intervals
# pick a sensible age for the max interval age
max<-FossilSim::basin.age(t)
# number of intervals
strata = 7
# sampling probability
sampling = 0.7
f <- FossilSim::sim.fossils.unif(t, max, strata, sampling)
plot(f, t, binned = TRUE, strata = strata, show.strata = TRUE)

## simulate fossils under a user-specified non-uniform preservation model & non-uniform intervals
# interval boundaries
times = c(0, 0.3, 1, 1.5, max)
# interval Poisson rates
rates = c(4, 1, 0, 1)
f <- FossilSim::sim.fossils.non.unif(t, times, rates)
plot(f, t, show.strata = TRUE, interval.ages = times)
# also plot proxy data
plot(f, t, show.strata = TRUE, interval.ages = times, show.proxy = TRUE, proxy.data = rates)

## simulate fossils under Steve Holland's non-uniform preservation model
# simulate water depth
wd <- FossilSim::sim.water.depth(strata)
f <- FossilSim::sim.fossils.non.unif.depth(t, wd, basin.age = max, strata = strata, PA = 1, PD = 0.5, DT = 1, convert.rate = TRUE)
plot(f, t, show.strata = TRUE, strata = strata, show.proxy = T, proxy.data = wd)
# plot ranges
plot(f, t, show.strata = TRUE, strata = strata, show.proxy = T, proxy.data = wd, show.ranges = T, show.fossils = F)

## count the number of fossils per intervals
FossilSim::count.fossils.binned(f, times)

#END
