# KT
sim2.rateshift <- function(n, age, lambda, mu, frac, times, norm){
    phy2 <- sim2.rateshift.origin(n, age, lambda, mu, frac, times, norm)
    if (class(phy2)) == "phylo"){
        phy2 <- collapse.singles(phy2)
    }
    phy2
}

