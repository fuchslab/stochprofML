mix.d.sum.of.mixtures.LNLN <- function(y, n.vector, p.vector, mu.vector, sigma.vector, logdens = T){
    densmix <- matrix(0, ncol = length(y), nrow = length(n.vector))
    for(i in 1: length(n.vector)){
        densmix[i, ] <- d.sum.of.mixtures.LNLN(y = y, n = n.vector[i], p.vector = p.vector, mu.vector = mu.vector, sigma.vector = sigma.vector, logdens = logdens)
    }
    Dens<-colSums(densmix)/length(n.vector)
}
