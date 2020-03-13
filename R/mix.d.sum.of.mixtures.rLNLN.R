mix.d.sum.of.mixtures.rLNLN <-
    function(y, n.vector, p.vector, mu.vector, sigma.vector){
# This function calculates a "mixed" density for a given n-vector and given data. So it
# calculates the density of a given random-variable vector y, with given cell numbers n.
# If each each value in y is a sum of n random variables from the following mixture
# distribution: With probability p_i, it is of type i.
# In that case, it is lognormally distributed with log-mean mu_i and log-standard
# deviation sigma_i.
# So in fact we calculate the density for y for each value given in the n-vector
# independently, assuming the value originates from a sum of this n-random variables
# from the mixture distribution. Then we take the weighted sum of these calculated densities.
#
# Parameters:
#
# - the density is evaluated at y
# - n.vector is the number of cells taken from each tissue sample. This can also be a vector stating how many
#   cells are in each sample seperatly.
# - p.vector is a vector containing the probabilities for each type of cell.
#   Its elements have to sum up to one.
# - mu.vector is the vector with mu-values for each type.
# - sigma.vector is the vector with sigma-values for each type.


densmix <- matrix(0, ncol = length(y), nrow = length(n.vector))
    for(i in 1: length(n.vector)){
        densmix[i, ] <- d.sum.of.mixtures.rLNLN(y = y, n = n.vector[i], p.vector = p.vector, mu.vector = mu.vector, sigma.vector = sigma.vector, logdens = FALSE)
    }
    Dens<-colSums(densmix)/length(n.vector)
}
