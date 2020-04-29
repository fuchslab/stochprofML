# stochprofML 2.0.2
* d.sum.of.lognormals: If it is not a real sum but only one summand use dlnorm directly
* in d.sum.of.types (all models) and correspondingly d.sum.of.lognormal.types: bug fix, as mu.vector and sigma.vector  were wrongly filled for TY > 2.
* small foramting changes on Help pages



# stochprofML 2.0.1
* Deleted the argument "logdens" in mix.d.sum.of.mixtures because of a bug if set to TRUE.

# stochprofML 2.0.0

* Added a `NEWS.md` file to track changes to the package.
* Created a github repository: https://github.com/fuchslab/stochprofML
* Extended n the number of cells that was constant over all samples to be flexible to be different for each sample.
* Added possibility to put in the predefined number of cells when using the random number generators r.sum.of.mixtures by specifying N.matrix
* Add the mixed density calculator when calculating the density of a distribution with a specific n.vector as input.
* Changed all functions that use the fixed n to be able to use a vector n
* Added doi of Tirier et al. (2019) paper where this was first implemented and used
* Minor fixes of textings
