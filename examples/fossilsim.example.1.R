library(devtools)
install_github("rachelwarnock/fossilsim")
# its a good idea to restart R after installing/updating packages via github or bitbucket

set.seed(1234)

#pdf("fossilSim.example.pdf",height=6, width=6)

## simulate birth death tree
t <- TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]

## simulate fossils under a Poisson sampling process
# Poisson sampling rate
sampling = 3
# simulate fossils
f <- FossilSim::sim.fossils.poisson(t, sampling)
p1<-plot(f, t) # plot 1

## simulate fossils under Steve Holland's non-uniform sampling model
# number of intervals (user specified intervals also possible)
strata = 8 
# pick a sensible age for the max interval age
max<-FossilSim::basin.age(t)
# simulate water depth
wd <- FossilSim::sim.water.depth(strata)
# simulate fossils
f <- FossilSim::sim.fossils.non.unif.depth(t, wd, basin.age = max, strata = strata, PA = 1, PD = 0.5, DT = 1, convert.rate = TRUE)

# plot tree + fossils
plot(f, t, show.strata = TRUE, strata = strata) # plot 2
# add proxy data
plot(f, t, show.strata = TRUE, strata = strata, show.proxy = T, proxy.data = wd) # plot 3
# plot stratigraphic ranges
plot(f, t, show.strata = TRUE, strata = strata, show.proxy = T, proxy.data = wd, show.ranges = T, show.fossils = F) # plot 4

## simulate fossils under a user-specified non-uniform preservation model & non-uniform intervals
# interval boundaries
times = c(0, 0.3, 1, 1.5, max)
# interval Poisson rates
rates = c(4, 1, 0, 1)
# simulate fossils
f <- FossilSim::sim.fossils.non.unif(t, times, rates)
# also plot proxy data
plot(f, t, show.strata = TRUE, interval.ages = times, show.proxy = TRUE, proxy.data = rates) # plot 5

#dev.off()
#END