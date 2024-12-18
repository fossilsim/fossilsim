## FossilSim

R package for simulating fossil data on phylogenetic trees under mechanistic models of preservation and sampling

The latest version can be installed in R using the package devtools:

    install.packages("BiocManager")
    BiocManager::install("ggtree")
    library(devtools)
    install_github("fossilsim/fossilsim")

### Quick start

Simulating data using `FossilSim` can be as simple as the following code snippets.

```{r}
# simulate a tree using ape
tips = 8
t = ape::rtree(tips)
# simulate fossils using fossilsim
rate = 2
f = sim.fossils.poisson(rate, t)  
# plot the complete output
plot(f, t)
# plot the reconstructed output
plot(f, t, reconstructed = TRUE)
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

* [General introduction to the package](https://CRAN.R-project.org/package=FossilSim/vignettes/intro.html)
* [Simulating taxonomy](https://CRAN.R-project.org/package=FossilSim/vignettes/taxonomy.html)
* [Simulating fossils](https://CRAN.R-project.org/package=FossilSim/vignettes/fossils.html)
* [Simulating FBD trees](https://CRAN.R-project.org/package=FossilSim/vignettes/simfbd.html)
* [Exporting sampled ancestor trees](https://CRAN.R-project.org/package=FossilSim/vignettes/SAtree.html)
* [Converting to paleotree format](https://CRAN.R-project.org/package=FossilSim/vignettes/paleotree.html)

### For further information and examples see the package documentation

### Contributors
* Joëlle Barido-Sottani  
* Walker Pett
* Joseph O'Reilly
* Ugnė Stolz
* Rachel Warnock

### Citation

Joëlle Barido-Sottani et al. 2019. *Methods Evolution & Ecology*. [FossilSim: an R package for simulating fossil occurrence data under mechanistic models of preservation and recovery](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13170)

