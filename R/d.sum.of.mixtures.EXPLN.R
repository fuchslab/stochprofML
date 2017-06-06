d.sum.of.mixtures.EXPLN <-
function(y,n,p.vector,mu.vector,sigma.vector,lambda,logdens=T) {
#
# Density of a sum of independent random variables, where each random variable is
# from the following mixture distribution. This is a mixture of zero, one or more
# lognormal distributions and one exponential distribution.
# More specifically, with probability p_i (i=1,...,T-1),
# a summand is of type i. In that case, it is lognormally
# distributed with log-mean mu_i and log-standard deviation sigma_i.
# With probability p_T, it is exponentially distributed with rate lambda.
#
# Parameters:
#
# - the density is evaluated at y, which can be multi-dimensional
# - n is the number of cells taken from each tissue sample. This can also be a vector stating how many
#   cells are in each sample seperatly.
# - p.vector is a vector containing the probabilities for each type of cell.
#   Its elements have to sum up to one.
# - mu.vector is the vector with mu-values for each lognormal type.
# - sigma.vector is the vector with sigma-values for each lognormal type.
# - lambda is the rate for the exponential type.
# - if logdens==T, the log of this density is returned

# The lengths mu.vector and sigma.vector have to be identical.
# p.vector has to have one component more.
# lambda has to be a scalar.
# These lengths automatically determine the number of different types.


   # check for some obvious errors
   if (round(sum(p.vector),4)!=1) {
      stop("d.sum.of.mixtures: Sum of p's does not equal one.")
   }
   if (sum(p.vector>1)+sum(p.vector<0)>0) {
      stop("d.sum.of.mixtures: There are p's which are not in [0,1]")
   }
   if ((length(p.vector)!=length(mu.vector)+1) || (length(p.vector)!=length(sigma.vector)+1)) {
      stop("d.sum.of.mixtures: p and mu and/or sigma are of contradicting lengths.")
   }
   if (length(lambda)!=1) {
      stop("d.sum.of.mixtures: lambda is not a scalar.")
   }


    if(length(n) != 1 && length(n) != dim(y)[1]){
        stop("d.sum.of.mixtures: n has to be either a finite natural number or vector, having the same length as the sample size.")
    }

    if(is.null(dim(y))) {y <- matrix(y, ncol =1)}

    TY <- length(p.vector)

   this.sum <- 0
   for(k in sort(unique(n))){
       index <- which(n == k)


   # all possible combinations for how many of the n random variables are of which type
       j.combis <- comb.summands(k,TY)


       for (i in 1:nrow(j.combis))  {
            this.j <- j.combis[i,]
            weight <- dmultinom(x=this.j,prob=p.vector,log=F)
            mixture.density <- lognormal.exp.convolution(y[index,],this.j,mu.vector,sigma.vector,lambda,logdens=F)
            this.sum <- this.sum + weight * mixture.density
       }
   }

   # sum up over all these possibilities
   if (logdens) { log(this.sum) } else { this.sum }
}
