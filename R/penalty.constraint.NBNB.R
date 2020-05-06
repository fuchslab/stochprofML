penalty.constraint.NBNB <-
function(dataset,parameter,smoothingpar=10^5) {
   # this function should only be called if TY>1
   m <- ncol(dataset)
   TY <- (length(parameter)+1)/(2*m+1)
   size <- parameter[TY:((m+1)*TY-1)]
   mu <- parameter[(m+1)*TY:length(parameter)]

   pen <- 0
   # build a matrix such that the g.th column contains the size values for gene g
   size <- matrix(size,byrow=T,ncol=m)
   # same for mu
   mu <- matrix(mu,byrow=T,ncol=m)
   # for all genes 1,....m and all types 1,...,TY-1:
   for (g in 1:m) {
      max.datapoint <- max(dataset[,g])
      for (i in 1:(TY-1)) {
         # mode of lognormal i
         mode.ig <- max(0,floor((mu[i,g]*(size[i,g]-1))/size[i,g]))
         if (mode.ig<max.datapoint) {
            # sequence of values where to compare the densities
            x <- seq(mode.ig,max.datapoint,(max.datapoint-mode.ig)/500)
            # add mode of population i+1 to this sequence if it is on the right of the mode of population i
            mode.igplus1 <- max(0,floor((mu[i+1,g]*(size[i+1,g]-1))/size[i+1,g]))
            if (mode.igplus1>mode.ig) {
               x <- c(x,mode.igplus1)
               x <- x[order(x)]
            }
            # density of nbs i at these values
            d.NB1 <- d_snb(y = x, size = size[i, g], mu = mu[i, g])
            # density of nbs i+1 at these values
            d.NB2 <- d_snb(y=x, size = size[i+1, g], mu = mu[i+1, g])
            # penalize if d.NB2>d.NB1
            pen <- pen + sum(pmax(0,d.NB2-d.NB1)^2)
            # now check whether population i+1 is peaked; in that case, the function d.sum.of.lognormals
            # most probably smoothed the density too much
            if (sum(d.NB2)==0) {
               if (mode.igplus1>mode.ig) {
                  pen <- pen + max(0,dnbinom(mode.igplus1,size = size[i+1,g],mu = mu[i+1,g])-dnbinom(mode.igplus1,size= size[i,g], mu= mu[i,g]))^2
               }
            }
         }
      }
   }
   return(smoothingpar*pen)
}
