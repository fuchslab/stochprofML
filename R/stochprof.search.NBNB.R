stochprof.search.NBNB <-
function(dataset,n,TY,method="grid",M=10,par.range=NULL,prev.result=NULL,fix.mu=F,fixed.mu,genenames=NULL,print.output=F,use.constraints=F) {
# Calculates the log likelihood function of all model parameters for a given dataset
# at certain parameter values. The so-obtained values are returned in a matrix with
# the following entries: Each row corresponds to one parameter combination. All columns
# but the last one contain the parameter values at which the log likelihood function has
# been computed. The column names are the parameter names. The last column ("target") is the
# negative log likelihood function computed at the respective parameter vector. For numerical
# reasons, this target value is set to the minimum of 10^7 and the actual value.
#
# The values at which the target function is calculated are randomly drawn from some range
# specified by "par.range". If method=="grid", the target function is simply evaluated
# at such a randomly drawn parameter vector. If method=="optim", this randomly drawn vector is
# passed to the Nelder-Mead algorithm as a starting value in order to search for a local
# maximum around it.
#
# Parameters:
#
# - dataset is a matrix which contains the cumulated expression data over all cells in a tissue sample.
#   Columns represent different genes, rows represent different tissue samples.
# - n is the number of cells taken from each tissue sample. This can also be a vector stating how many
#   cells are in each sample seperatly.
# - TY is the number of types of cells that is assumed in the stochastic model.
# - method (default="grid") determines whether a grid search or the Nelder-Mead algorithm should be applied:
#   If method=="grid", the log likelihood function is simply evaluated at certain parameter values that are
#   randomly drawn.
#   If method=="optim", a Nelder-Mead search starts at a randomly drawn set of parameter values in order to
#   find a local maximum. The resulting locally optimal parameter is stored in the results matrix as one row.
# - M (default=10) is the number of randomly drawn parameter combinations.
# - par.range (default=NULL) is the range from which the parameter values should be randomly drawn. This is
#   a matrix with the number of rows being equal to the number of model parameters. The first column contains
#   the lower bound, the second column the upper bound. If par.range==NULL, some rather large range is defined.
# - prev.result (default=NULL) can contain results from former calls of this function.
# - genenames (default=NULL) are the names of the genes in the dataset.
#   For genenames==NULL, the genes will simply be enumerated according to the column numbers in the dataset.
# - If print.output==T (default=F), interim results of the grid search and numerical optimization are printed
#   into the console throughout the estimation procedure.
# - If use.constraints==T, constraints on the densities of the populations will be applied.


   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   d.sum.of.mixtures <- NULL
   backtransform.par <- NULL
   penalty.constraint <- NULL
   draw.parameters <- NULL
   transform.par <- NULL
   rm(d.sum.of.mixtures)
   rm(backtransform.par)
   rm(penalty.constraint)
   rm(draw.parameters)
   rm(transform.par)


   if (M==0) return(NULL)

   ####################
   # general settings #
   ####################

   # number of genes
   m <- ncol(dataset)
   # gene names
   if (is.null(genenames)) {
      genenames <- 1:m
   }

   # names of variables
   if (TY==1) {
      varnames <- c(paste("size_",genenames,sep=""),
                    paste("mu_",genenames,sep="")
                    )
   }
   else {
      varnames <- paste("p_",1:(TY-1),sep="")
      for (i in 1:TY) {
         varnames <- c(varnames,paste("size",i,"gene",genenames,sep="_"))
      }
      for (i in 1:TY) {
          varnames <- c(varnames,paste("mu",i,"gene",genenames,sep="_"))
      }
   }

   ###########################################
   ## ML estimation: define target function ##
   ###########################################

   # loglikelihood for one gene
   loglikeli <- function(y,p,size,mu) {
      # p, size and mu are of length TY
      max(-10^7,sum(d.sum.of.mixtures(y, n, p, size, mu, logdens=T)))
   }

   # this function will be minimized (the function "to.minimize" below just
   # changes the parameterisation)
   target.function <- function(p,size,mu) {
      # p is of length TY
      # size is of length TY*m
      # mu is of length TY*m

      # build size.matrix such that the g.th row contains the values for gene g
      size.matrix <- matrix(size,byrow=F,ncol=TY)

      # build mu.matrix such that the g.th row contains the values for gene g
      mu.matrix <- matrix(mu,byrow=F,ncol=TY)

      # consider negative log likelihood because the target function will be minimized
      this.sum <- 0
      for (g in 1:m) {
         this.sum <- this.sum - loglikeli(dataset[,g, drop=F],p,size.matrix[g,,drop=T],mu.matrix[g,,drop=T])
      }
      return(this.sum)
   }

   # For identifiability purposes, we require the mu-values of the first gene to be
   # in descending order. Otherwise, the parameter combinations (p,size1_gene1,size2_gene1,mu1_gene1,mu2_gene1)
   # and (1-p,size2_gene1,size1_gene1,mu2_gene1,mu1_gene1) would yield identical values of the log-likelihood
   # function.
   # The randomly drawn parameters will always fulfil the above requirement. However, during the
   # Nelder-Mead optimization procedure it might happen that it is violated. In that case, the
   # optim algorithm might jump back and forth between two equivalent states.
   # In order to avoid this, a penalty term is introduced which is added to the actual target
   # function. This penalty is positive if mu_1^1 >= mu_2^1 >= ... >= mu_TY^1 is not fulfilled,
   # and 0 otherwise.
   penalty.mu <- function(mu,m,lambda=100) {
      # build a matrix such that the g.th column contains the mu values for gene g
      mu <- matrix(mu,byrow=T,ncol=m)
      # mu for gene 1
      mu.g1 <- mu[,1]
      # penalty
      differences <- mu.g1[-1]-mu.g1[-length(mu.g1)]
      pen <- pmax(0,differences)
      return(lambda*sum(pen^2))
   }

   # this function should be minimized
   to.minimize <- function(theta) {
      # theta=(w_1,...,w_{T-1},log(size), log(mu)).
      # Everything but size and mu is one-dimensional.
      # size and mu are both TY*m-dim.

      # backtransformation
      # afterwards, back.theta is full-dim, incl. mu
      back.theta <- backtransform.par(this.par=theta,m=m,fix.mu=fix.mu,fixed.mu=fixed.mu)

      if (TY>1) {
         p <- back.theta[1:(TY-1)]
         p <- c(p,1-sum(p))
      }
      else {
         p <- 1
      }
      size <- back.theta[TY:((m+1)*TY-1)]
      mu <- back.theta[((m+1)*TY):length(back.theta)]

      # penalty
      pen.mu <- 0
      if (TY>1) {
         pen.mu <- penalty.mu(mu,m)
      }
      pen.constr <- 0
      if ((TY>1) && (use.constraints)) {
         pen.constr <- penalty.constraint(dataset,parameter=c(p[-TY],size,mu))
      }

      # target
      a <- target.function(p,size,mu) + pen.mu + pen.constr
      a <- min(10^7,a)
      return(a)
   }

   ####################
   # previous results #
   ####################

   # are there previous results already?
   if (is.null(prev.result)) {
      # no, there aren't
      all.results <- matrix(nrow=0,ncol=length(varnames)+1)
      colnames(all.results) <- c(varnames,"target")
   }
   else {
      # yes, there are
      all.results <- prev.result
   }


   ###################################
   # parameter ranges to be searched #
   ###################################

   if (is.null(par.range)) {
      # determine some range
      ranges <- matrix(NA,ncol=2,nrow=length(varnames))
      # p
      if (TY>1) {
         ranges[1:(TY-1),1] <- 0
         ranges[1:(TY-1),2] <- 1
      }
      # mu
      ranges[TY:((m+1)*TY-1),1] <- 0.01
      ranges[TY:((m+1)*TY-1),2] <- 1000
      # sigma
      ranges[((m+1)*TY):nrow(ranges),1] <- 0.01
      ranges[((m+1)*TY):nrow(ranges),2] <- 1000
   }
   else {
      ranges <- par.range
   }


   ###################################
   ## estimation: optimization step ##
   ###################################

   #--------------#
   # optim method #
   #--------------#
   if (method=="optim") {
      for (i in 1:M) {
         # draw starting value
         par0 <- draw.parameters(ranges,m) # full-dim.

         if (print.output) {  # putput before optim starts so that I cans ee if something wents wrong in optim
            cat("---\n")
            cat("Start optim at:\n")
            cat(par0,"\n")
         }

         theta0 <- transform.par(this.par=par0,m=m)

         theta0[theta0==-Inf] <- -10^7
         theta0[theta0==Inf] <- 10^7

         # numerically optimize
         if (length(theta0)==1) {
            result <- optim(theta0,fn=to.minimize,control=list(maxit=10^5),method="Brent",hessian=F, lower=-10^7, upper=10^7)
         }
         else {
            result <- optim(theta0,fn=to.minimize,control=list(maxit=10^5),method="Nelder-Mead",hessian=F)
         }

         # result
         this.theta <- result$par
         this.par <- backtransform.par(this.par=this.theta,m=m)
         this.value <- result$value

         # attach new result to all former ones
         all.results <- rbind(all.results,c(this.par,this.value))

         if (print.output) {
            cat("Arrived at:\n")
            cat(this.par,"\n")
         }
      }
   }
   #-------------#
   # grid search #
   #-------------#
   else if (method=="grid") {
      for (i in 1:M) {
         # randomly draw parameter
         this.par <- draw.parameters(ranges,m) # full-dim
         this.theta <- transform.par(this.par=this.par,m=m)

         if (print.output) {
            cat("---\n")
            cat("Compute grid at:\n")
            cat(this.par,"\n")
         }

         # target function
         this.value <- to.minimize(this.theta)
         # attach new result to all former ones
         all.results <- rbind(all.results,c(this.par,this.value))
      }
   }

   return(all.results)
}
