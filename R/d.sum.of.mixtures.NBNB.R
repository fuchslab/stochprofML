d.sum.of.mixtures.NBNB <-
    function(y, n, p.vector, size.vector, mu.vector, logdens = TRUE) {
        # Density of a sum of independent random variables, where each random variable is
        # from the following mixture distribution: With probability p_i, it is of type i.
        # In that case, it is lognormally distributed with log-mean mu_i and log-standard
        # deviation sigma_i.
        #
        # Parameters:
        #
        # - the density is evaluated at y
        # - n is the number of cells taken from each tissue sample. This can also be a vector stating how many
        #   cells are in each sample seperatly.
        # - p.vector is a vector containing the probabilities for each type of cell.
        #   Its elements have to sum up to one.
        # - size.vector is the vector with size-values for each type.
        # - mu.vector is the vector with mu-values for each type.
        # - if logdens==T, the log of this density is returned

        # definition of variables (necessary for CMD R check)
        # (these variables will be initialized later, but they are not visible as global functions/data)
        d.sum.of.types <- NULL
        rm(d.sum.of.types)


        # check for some obvious errors
        if (round(sum(p.vector),4)!=1) {
            stop("d.sum.of.mixtures: Sum of p's does not equal one.")
        }
        if (sum(p.vector>1)+sum(p.vector<0)>0) {
            stop("d.sum.of.mixtures: There are p's which are not in [0,1]")
        }
        if ((length(p.vector)!=length(size.vector)) || (length(p.vector)!=length(mu.vector))) {
            stop("d.sum.of.mixtures: p and size and/or mu are of different lengths.")
        }


        if(length(n) != 1 && length(n) != length(y)){
            stop("d.sum.of.mixtures: n has to be either a finite natural number or vector, having the same length as the sample size.")
        }

        #if n is the same in all samples calculate as in previous stochasticML versions
        if(sum(n == n[1]) == length(n)){
            n <- n[1]
            j.combis <- comb.summands(n,length(p.vector))

            # sum up over all these possibilities
            this.sum <- 0
            for (i in 1:nrow(j.combis))  {
                this.j <- j.combis[i,]
                weight <- dmultinom(x=this.j,prob=p.vector,log=F)
                mixture.density <- d.sum.of.types(y,this.j,size.vector,mu.vector,logdens=F)
                this.sum <- this.sum + weight * mixture.density
            }
        } else{
            this.sum <- rep(0, length(y))
            for(k in sort(unique(n))){
                index <- which(n == k)


                # all possible combinations for how many of the n random variables are of which type
                j.combis <- comb.summands(k,length(p.vector))


                # sum up over all these possibilities

                for (i in 1:nrow(j.combis))  {
                    this.j <- j.combis[i,]
                    weight <- dmultinom(x=this.j,prob=p.vector,log=F)
                    mixture.density <- d.sum.of.types(y[index],this.j,size.vector,mu.vector,logdens=F)
                    this.sum[index] <- this.sum[index] + weight * mixture.density
                }
            }

        }
        # all possible combinations for how many of the n random variables are of which type





        if (logdens) { log(this.sum) } else { this.sum }
    }
