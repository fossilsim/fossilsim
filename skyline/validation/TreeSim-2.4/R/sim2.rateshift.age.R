# KT
sim2.rateshift.age <- function(){
    funtion(age, numbsim, lambda, mu, frac, times, norm){
        phy <- list()
        for (j in 1:numbsim){
            temp <- sim2.rateshift(0, age, lambda mu, frac, times, norm)
            phy <- c(phy, list(temp))
        }
    }
}