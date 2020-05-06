transform.par.NBNB <-
function(this.par,m) {
# Transforms the parameter from its original scale to all real numbers, i.e.
# INPUT:
#    this.par=(p1,p2,...,p_(T-1),size,mu)
#
# OUTPUT:
#    this.theta=(w1,w2,...w_(T-1),log(size),log(mu))
# with
#    w_i = logit((p_1+...+p_i)/(p_1+...+p_(i+1))) for i=1,...,T-1
# and
#    size = (size_1_gene_1,size_2_gene_1,...,size_T_gene_1,...,size_1_gene_2,size_2_gene_2,...,size_T_gene_2)
# and
#    mu = (mu_1_gene_1,mu_2_gene_1,...,mu_T_gene_1,...,mu_1_gene_2,mu_2_gene_2,...,mu_T_gene_2)
#
#
# - this.par is the parameter on the original scale.
# - m is the number of genes analyzed

#



   # determine number of types of cells
   TY <- (length(this.par)+1)/(2*m+1)

   # (p_1,...,p_(T-1)) --> w
   if (TY>1) {
      this.p <- this.par[1:(TY-1)]
      this.w <- rep(NA,length(this.p))
      if (TY>2) {
         for (i in 1:(length(this.w)-1)) {
            num <- sum(this.p[1:i])
            denom <- sum(this.p[1:(i+1)])
            if (denom!=0) {
               this.w[i] <- stochprof.logit(num/denom)
            }
            else {
               this.w[i] <- stochprof.logit(0)
            }
         }
      }
      this.w[length(this.w)] <- stochprof.logit(sum(this.p))
   }
   else {
      this.w <- NULL
   }



   # size
   this.size <- this.par[(TY):((m+1)*TY-1)]


   # mu
   this.mu <- this.par[((m+1)*TY):((2*m+1)*TY-1)]

   # concatenate, transform
   this.theta <- c(this.w,log(this.size),log(this.mu))

   return(this.theta)
}
