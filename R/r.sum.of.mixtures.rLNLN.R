r.sum.of.mixtures.rLNLN <-
function(k,n,p.vector,mu.vector,sigma.vector,N.matrix = NULL) {
# Draws k i.i.d. random variables from the following distribution:
# Each random variable is the sum of another n independent random variables.
# These are from a mixture of lognormal distributions. More specifically,
# with probability p_i, a summand is of type i. In that case, it is lognormally
# distributed with log-mean mu_i and log-standard deviation sigma_i.
#
# Parameters:
#
# - k is the number of i.i.d. random variables returned by this function
#   (in the considered application: number of tissue samples)
# - n is the number of summands entering each of the k random variables. This can also be a vector stating how many
#   summands enter each random variable seperatly. (in the considered application: number of cells per tissue sample)
# - p.vector=(p1,p2,..,pT) is the probability for each type
# - mu.vector=(mu1,mu2,...,muT) is the log-mean for each type
# - sigma.vector=(sigma1,...,sigmaT) is the log-standard deviation for each type
# - N.matrix is the possible input of the decomposition of random variables coming from the different types
#
# The lengths of p.vector, mu.vector and sigma.vector have to be identical.
# Their lengths automatically determine the number of different types.

   # check for some obvious errors
   if (round(sum(p.vector),4)!=1) {
      stop("r.sum.of.mixtures: Sum of p's does not equal one.")
   }
   if (sum(p.vector>1)+sum(p.vector<0)>0) {
      stop("r.sum.of.mixtures: There are p's which are not in [0,1]")
   }
   if ((length(p.vector)!=length(mu.vector)) || (length(p.vector)!=length(sigma.vector))) {
      stop("r.sum.of.mixtures: p and mu and/or sigma are of different lengths.")
   }
   if(length(n) != 1 && length(n) != k) {
      stop("r.sum.of.mixtures: n has to be either a finite natural number or vector, having the same length as the sample size.")
      }
   if(!is.null(N.matrix)){
      if(!(dim(N.matrix) == c(k,length(p.vector)) && rowSums(N.matrix) == n )){
         stop("N.matrix has to be a matrix, where each column describes the decomposition of one observation, i.e.
           the number of cells of each contained population. Populations are shown in the rows ")
      }}

   # draw the number of summands of type i (i.th column) in each sample j (j.th row)
   # and check if N.matrix is given as an argument!
   if(is.null(N.matrix)){
      if(length(n) == 1){
         N.matrix <- t(rmultinom(n=k, size=n, prob=p.vector))
         n<-rep(n, k)
      } else {
         N.matrix <- matrix(0, nrow=k, ncol=length(p.vector))
         for(j in 1:k){
            N.matrix[j,] <- t(rmultinom(n=1, size=n[j], prob=p.vector))
         }
      }
   }
   if(length(n) == 1){
     n<-rep(n, k)
     }

   # draw the summands and sum up
   random <- rep(NA,k)
   for (j in 1:k) {
      # j.th sample
      N.vector <- N.matrix[j,]
      # expand mu.vector and sigma.vector
      full.mu.vector <- rep(0,n[j])
      full.sigma.vector <- full.mu.vector
      index <- 0
      for (i in 1:length(p.vector)) {
         if (index<n[j]) {
            full.mu.vector[index+(1:N.vector[i])] <- mu.vector[i]
            full.sigma.vector[index+(1:N.vector[i])] <- sigma.vector[i]
            index <- index + N.vector[i]
         }
      }
      # vector of summands
      Y <- rlnorm(n[j],full.mu.vector,full.sigma.vector)
      # sum up
      random[j] <- sum(Y)
   }

   return(random)
}
