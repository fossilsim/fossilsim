library(devtools)
install_github("rachelwarnock/fossilsim")

set.seed(123)

## simulate birth death tree
t<-TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]

## simulate fossils under a Poisson sampling process
# Poisson sampling rate
sampling = 3
f<-sim.fossils.poisson(t, sampling)
plot(f, t)

## simulate fossils under a unifom sampling model with equal length intervals
# pick a sensible age for the max interval age
max<-basin.age(t)
# number of intervals
strata = 4
# sampling probability
sampling = 0.7
f <- sim.fossils.unif(t, max, strata, sampling)
plot(f, t, binned = TRUE, strata = strata)

## simulate fossils under a user-specified non-uniform preservation model
# interval boundaries
times = c(0, 0.3, 1, 1.5, max)
# interval Poisson rates
rates = c(4, 1, 0, 1)
f<-sim.fossils.non.unif(t, times, rates)
plot(f, t, show.strata = TRUE, interval.ages = times)
# also plot proxy data
plot(f, t, show.strata = TRUE, interval.ages = times, show.proxy = T, proxy.data = rates)

## simulate fossils under Steve Holland's non-uniform preservation model
# simulate water depth
wd<-sim.water.depth(strata)

f <- sim.fossils.non.unif.depth(t, wd, times, PA = 1, PD = 0.5, DT = 1)
plot(f, t, show.strata = TRUE, interval.ages = times, show.proxy = T, proxy.data = wd, binned = TRUE)
