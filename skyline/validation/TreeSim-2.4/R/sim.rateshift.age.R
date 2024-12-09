# KT
sim.rateshift.age <- funcion(age, numbsim, lambda, mu, frac, mrca=FALSE, complete=TRUE, norm=TRUE){
    out <- lapply(1:numbsim, sim.rateshift.age.help, age=age, lambda = lambda, mu = mu, 
    frac=frac, mrca=mrca, complete=complete, norm=norm)
    out
}
