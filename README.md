## FossilSim

R package for simulating fossil data on phylogenetic trees under mechanistic models of preservation and sampling

The latest version can be installed in R using the package devtools:

    library(devtools)
    install_github("rachelwarnock/fossilsim")

### Quick start

Simulating data using `FossilSim` can be as simple as the following code snippets.

```{r}
# simulate a tree using ape
tips = 8
t = ape::rtree(tips)
# simulate fossils using fossilsim
rate = 2
f = sim.fossils.poisson(rate, t)  
# plot the output
plot(f, t)
```

```{r}
# simulate taxonomy using fossilsim
beta = 0.5 # probability of symmetric speciation
lambda.a = 0.1  # rate of anagenesis
s = sim.taxonomy(t, beta, lambda.a)  
# plot the output
plot(s, t, legend.position = "bottomright")
```

### Package vignettes

The following vignettes are available via CRAN and provide detailed examples:

* [General introduction to the package](https://cran.r-project.org/web/packages/FossilSim/vignettes/intro.html)
* [Simulating taxonomy](https://cran.r-project.org/web/packages/FossilSim/vignettes/taxonomy.html)
* [Simulating fossils](https://cran.r-project.org/web/packages/FossilSim/vignettes/fossils.html)
* [Simulating FBD trees](https://cran.r-project.org/web/packages/FossilSim/vignettes/simfbd.html)
* [Exporting sampled ancestor trees](https://cran.r-project.org/web/packages/FossilSim/vignettes/SAtree.html)
* [Converting to paleotree format](https://cran.r-project.org/web/packages/FossilSim/vignettes/paleotree.html)

### For further information and examples see the package documentation

### Contributors
Joëlle Barido-Sottani  
Walker Pett  
Rachel Warnock


