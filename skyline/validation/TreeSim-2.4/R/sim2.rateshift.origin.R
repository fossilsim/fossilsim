# KT
sim2.rateshift.origin <- function(n, age, lambdavec, muvec, frac, times, norm = TRUE) {
    check <- 0
    sumratiolam <- 0
    for (j in 1:length(lambdavec)) {
        if (lambdavec[j] > 0) {
            sumratiolam <- sumratiolam + 1 / lambdavec[j]
        }
    }
    while (check == 0) {
        lambda <- lambdavec[1]
        lambda0 <- lambda
        edge <- c(-1, -2) # matrix of edges
        leaves <- c(-2) # list of extant leaves
        timecreation <- c(0, 0) # time when species -2, -3, ... -n was created after origin
        extinct <- vector() # list of extinct leaves
        time <- 0 # time after origin
        interval <- 2
        doneME <- 0
        maxspecies <- -2 # smallest species
        edge.length <- c(0) # edge length. if 0: leaf which didn't speciate /extinct yet
        extincttree <- 0
        stop <- 0

        if (K > 0) {
            lambda <- max(lambda0 * (1 - length(leaves) / K), 0)
        }
        mu <- muvec[1]
        if (length(frac) == 1) {
            doneME <- 1
            nextextinction <- 0
        } else {
            nextextinction <- times[interval]
        }
        while (stop == 0) {
            if (length(leaves) == 0) {
                # phy2 = 0
                if (age > 0) {
                    phy2 <- 0
                }
                extincttree <- 1
                stop <- 1
            } else {
                if ((lambda + mu) == 0) {
                    time <- nextextinction
                } else {
                    timestep <- rexp(1, (length(leaves) * (lambda + mu))) # time since last event
                }
                if ((age > 0 && (time + timestep) < age) || age == 0) {
                    # if (K > 0) {
                    #     lambda <- lambda0 * (1 - length(leaves) / K)
                    #     if (lambda < 0) {
                    #         lambda <- 0
                    #     }
                    # }
                    time <- time + timestep # time after origin
                    species <- sample(leaves, 1) # the leaf undergoing the next event
                    del <- which(leaves == species)
                    specevent <- runif(1, 0, 1) # event speciation or extinction
                    edgespecevent <- which(edge == species) - length(edge.length)
                    if ((lambda / (lambda + mu)) > specevent) {
                        edge.length[edgespecevent] <- time - timecreation[-species]
                        edge <- rbind(edge, c(species, maxspecies - 1))
                        edge <- rbind(edge, c(species, maxspecies - 2))
                        edge.length <- c(edge.length, 0, 0)
                        leaves <- c(leaves, maxspecies - 1, maxspecies - 2)
                        maxspecies <- maxspecies - 2
                        leaves <- leaves[-del]
                        # if (K > 0){
                        #     lambda <- max(lambda0*(1-length(leaves)/K),0)
                        # }
                        timecreation <- c(timecreation, time, time)
                        if (length(leaves) == n) {
                            stop <- 1
                        }
                    } else {
                        extinct <- c(extinct, leaves[del])
                        leaves <- leaves[-del]
                        edge.length[edgespecevent] <- time - timecreation[-species]
                    }
                } else {
                    stop <- 1
                }
            }
        }

        ## TODO: Insert rateshift changes here / check placement and contents vs. sim2.bd.rateshift.single
        nadd <- round(length(leaves) / frac[interval])
        if (nadd > length(leaves)) {
            for (j in 1:(nadd - length(leaves))) {
                timecreation <- c(timecreation, time)
                extinct <- extinct + 1
                leaves <- c(leaves, (maxleaf + extinct))
                nodes <- c(nodes, (maxleaf + extinct))
            }
        }
        lambda <- lambdavec[interval]
        lambda0 <- lambda
        mu <- muvec[interval]
        if (interval == length(frac)) {
            doneME <- 1
        } else {
            interval <- interval + 1
            nextextinction <- times[interval]
        }
        if (norm == TRUE) {
            sampletree <- runif(1, 0, 1)
            if (sampletree < (1 / lambda / (sumratiolam))) {
                check <- 1
            }
        } else {
            check <- 1
        }


        if (extincttree == 0) { # length pendant edge
            if (age == 0) {
                timestep <- rexp(1, (length(leaves) * (lambda + mu))) # time since last event
                time <- time + timestep # time after origin
            } else {
                time <- age
            }
        }


        if (extincttree == 0 || age == 0) {
            # assign pendant edge length
            for (j in (1:length(leaves))) {
                k <- which(edge == leaves[j]) - length(edge.length)
                edge.length[k] <- time - timecreation[-leaves[j]]
            }

            nodes <- (length(leaves) + length(extinct)) * 2
            leaf <- 1
            interior <- length(leaves) + length(extinct) + 1
            edgetemp <- edge

            if (nodes == 2) {
                # edge = c(2,1)
                # leaves = c(1,2)
                phy2 <- 1
            } else {
                for (j in (1:nodes)) {
                    if (sum(match(leaves, -j, 0)) + sum(match(extinct, -j, 0)) == 0) {
                        # replace all -j values in edge by interior
                        posvalue <- interior
                        interior <- interior + 1
                    } else {
                        posvalue <- leaf
                        leaf <- leaf + 1
                    }
                    replacel <- which(edge == -j)
                    if (length(replacel) > 0) {
                        for (k in 1:length(replacel)) {
                            if ((replacel[k] - 1) < length(edge.length)) {
                                edge[replacel[k], 1] <- posvalue
                            } else {
                                edge[(replacel[k] - length(edge.length)), 2] <- posvalue
                            }
                        }
                    }
                }


                phy <- list(edge = edge)
                phy$tip.label <- paste("t", sample(length(leaves) + length(extinct)), sep = "")
                phy$edge.length <- edge.length
                phy$Nnode <- length(leaves) + length(extinct)
                class(phy) <- "phylo"
                # phy2<-collapse.singles(phy)
                phy2 <- phy
            }
        }

        phy2
    }
}
