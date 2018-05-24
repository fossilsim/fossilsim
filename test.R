## simulate tree
t = TreeSim::sim.bd.taxa(8, 1, 1, 0.3)[[1]]

## simulate fossils under a Poisson sampling process
f = sim.fossils.poisson(rate = 3, tree = t)


# possible options
show.fossils = c(TRUE, FALSE)
show.tree = c(TRUE, FALSE)
show.ranges = c(TRUE, FALSE)
# age info
show.strata = c(TRUE, FALSE)
binned = c(TRUE, FALSE)
show.axis = c(TRUE, FALSE)
strata = 5
max = 4.8
show.proxy = c(TRUE, FALSE)

#proxy.data = NULL, show.preferred.environ = FALSE,
#preferred.environ = NULL,

#interval.ages = NULL

# proxy stuff

# taxonomy
#show.taxonomy = FALSE, species = NULL,
# tree appearance
#root.edge = TRUE, hide.edge = FALSE, edge.width = 1,
#show.tip.label = FALSE, align.tip.label = FALSE, fcex = 1.2,
# fossil appearance
#fcol = "darkorange", ecol = NULL

opts = expand.grid(show.fossils = show.fossils, show.tree = show.tree, show.ranges = show.ranges,
                   # age stuff
                   show.strata = show.strata, binned = binned, show.axis = show.axis, strata = strata, max = max,
                   # proxy stuff
                   show.proxy = show.proxy)

#strata = strata, max = max, interval.ages = interval.ages,

for(i in 1:length(opts[,1])){
  plot(f, t, show.fossils = opts["show.fossils"][i,],
       show.tree = opts["show.tree"][i,],
       show.ranges = opts["show.ranges"][i,],
       show.strata = opts["show.strata"][i,],
       binned = opts["binned"][i,],
       show.axis = opts["show.axis"][i,],
       strata = opts["strata"][i,],
       #max = opts["max"][i,],
       show.proxy = opts["show.proxy"][i,],
       proxy.data = c(1,2,3,1,2)
       )
}


#what's with the tiny fossils?
#can max override yy coordinates?




