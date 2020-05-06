d.sum.of.types.NBNB <-
    function(y,j.vector,size.vector,mu.vector,logdens=T) {

        # Density for a sum of independent lognormally distributed variables,
        # using the approximation method by Fenton [1].
        #
        # - density is evaluated at y
        #
        # - j.vector=(j1,...,jT) is a vector indicating how many of the summands are of
        #   which type:
        #         j1 are of type I, ..., jT are of type T.
        #   The sum n=j1+...+jT implies how many summands are entering the sum.
        #
        # - size.vector=(size1,...,sizeT) is of the same length as j.vector. size_i is the
        #   size for type i.
        #
        # - mu.vector is defined analogously as size.vector. It contains the standard
        #   deviations (not variances!) of the according normal distributions.
        #
        # - if logdens==T, the log of this density is returned.
        #
        #



        # check for some obvious errors
        if (sum(j.vector)==0) {
            stop("d.sum.of.types: Sum of j's equals zero.")
        }
        if ((length(j.vector)!=length(size.vector)) || (length(j.vector)!=length(mu.vector))) {
            stop("d.sum.of.types: j and size and/or mu are of different lengths.")
        }

        # these are the "extended" mu and sigma vector, containing
        # j1 times mu1 (or sigma1, resp), j2 times mu2 (or sigma2, resp), ...
        # change to NB setting, reduce vector, as sum of j the same NB is again NB with size = j*size and mu = j*mu
        # if NB is defined via prob the prob would stay the same


        full.size.vector <- size.vector * j.vector
        full.mu.vector <- mu.vector *  j.vector

        if (sum(j.vector == 0) > 0){
            del <- which(j.vector == 0)
            full.size.vector <- full.size.vector[,-del, drop = FALSE]
            full.mu.vector <- full.mu.vector[,-del, drop = FALSE]
        }


        # pass these extended mu and sigma vectors to the actual approximation
        # of the density of sums of lognormals
        d <- d_snb(y, size_param = full.size.vector, mu_param = full.mu.vector)
        if (logdens) { log(d) } else { d }
    }
