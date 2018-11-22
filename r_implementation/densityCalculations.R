
## densityCalculations.R
# implementation of the densities used in the gaussian mixture problem



# Mixture density using approximate predictive gaussian distribution
mixture_pdf_gaussian = function(model, data) {
    mixture = vector(mode = "numeric", length = NROW(data))
    for (k in 1:length(model$nu)) {
        tau_k   = model$W[,,k] * model$nu[k] # TODO: Is this right?
        mu_k    = model$m[, k]
        mixture = mixture + model$pi_k[k] * 
            dmvnorm(cbind(data$x,data$y), mean = mu_k, sigma = solve(tau_k))
    }
    return(mixture)
} # end of mixture_pdf_gaussian() function

# Mixture density using predictive t-distribution
mixture_pdf_t = function(model, data){
    mixture = vector(mode = "numeric", length = NROW(data))
    for (k in 1:length(model$nu)) {
        L_k = solve((((model$nu[k] + 1 - NROW(model$m)) * model$beta[k]) / 
                          (1 + model$beta[k])) * model$W[,,k])
        mixture = mixture + (model$alpha[k]/sum(model$alpha)) * 
            dmvt(x = cbind(data$x,data$y), delta = model$m[, k], 
                 sigma = L_k, df = model$nu[k] + 1 - NROW(model$m), 
                 log = FALSE, type = "shifted")
    }
    return(mixture)
} # end of mixture_pdf_t() function


# end of densityCalculations.R
