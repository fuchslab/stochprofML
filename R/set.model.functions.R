set.model.functions <-
function(model) {
# Defines some global model-dependent functions.

   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   d.sum.of.mixtures <- NULL
   backtransform.par <- NULL
   stochprof.results <- NULL
   transform.par <- NULL
   calculate.ci <- NULL
   d.sum.of.types <- NULL
   draw.parameters <- NULL
   get.range <- NULL
   penalty.constraint <- NULL
   r.sum.of.mixtures <- NULL
   stochprof.search <- NULL
   mix.d.sum.of.mixtures <- NULL
   rm(d.sum.of.mixtures)
   rm(backtransform.par)
   rm(stochprof.results)
   rm(transform.par)
   rm(calculate.ci)
   rm(d.sum.of.types)
   rm(draw.parameters)
   rm(get.range)
   rm(penalty.constraint)
   rm(r.sum.of.mixtures)
   rm(stochprof.search)
   rm(mix.d.sum.of.mixtures)


   if (model=="LN-LN") {
      backtransform.par <<- backtransform.par.LNLN
      calculate.ci <<- calculate.ci.LNLN
      d.sum.of.mixtures <<- d.sum.of.mixtures.LNLN
      mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.LNLN
      d.sum.of.types <<- d.sum.of.types.LNLN
      draw.parameters <<- draw.parameters.LNLN
      get.range <<- get.range.LNLN
      penalty.constraint <<- penalty.constraint.LNLN
      r.sum.of.mixtures <<- r.sum.of.mixtures.LNLN
      stochprof.results <<- stochprof.results.LNLN
      stochprof.search <<- stochprof.search.LNLN
      transform.par <<- transform.par.LNLN
      mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.LNLN
   }
   else if (model=="rLN-LN") {
      backtransform.par <<- backtransform.par.rLNLN
      calculate.ci <<- calculate.ci.rLNLN
      d.sum.of.mixtures <<- d.sum.of.mixtures.rLNLN
      mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.rLNLN
      d.sum.of.types <<- d.sum.of.types.rLNLN
      draw.parameters <<- draw.parameters.rLNLN
      get.range <<- get.range.rLNLN
      penalty.constraint <<- penalty.constraint.rLNLN
      r.sum.of.mixtures <<- r.sum.of.mixtures.rLNLN
      stochprof.results <<- stochprof.results.rLNLN
      stochprof.search <<- stochprof.search.rLNLN
      transform.par <<- transform.par.rLNLN
      mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.rLNLN
   }
   else if (model=="EXP-LN") {
      backtransform.par <<- backtransform.par.EXPLN
      calculate.ci <<- calculate.ci.EXPLN
      d.sum.of.mixtures <<- d.sum.of.mixtures.EXPLN
      mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.EXPLN
      draw.parameters <<- draw.parameters.EXPLN
      get.range <<- get.range.EXPLN
      penalty.constraint <<- penalty.constraint.EXPLN
      r.sum.of.mixtures <<- r.sum.of.mixtures.EXPLN
      stochprof.results <<- stochprof.results.EXPLN
      stochprof.search <<- stochprof.search.EXPLN
      transform.par <<- transform.par.EXPLN
      mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.EXPLN
   }
   else if (model=="NB-NB") {
       backtransform.par <<- backtransform.par.NBNB
       calculate.ci <<- calculate.ci.NBNB
       d.sum.of.mixtures <<- d.sum.of.mixtures.NBNB
       mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.NBNB
       d.sum.of.types <<- d.sum.of.types.NBNB
       draw.parameters <<- draw.parameters.NBNB
       get.range <<- get.range.NBNB
       penalty.constraint <<- penalty.constraint.NBNB
       r.sum.of.mixtures <<- r.sum.of.mixtures.NBNB
       stochprof.results <<- stochprof.results.NBNB
       stochprof.search <<- stochprof.search.NBNB
       transform.par <<- transform.par.NBNB
       mix.d.sum.of.mixtures <<- mix.d.sum.of.mixtures.NBNB
   }
}
