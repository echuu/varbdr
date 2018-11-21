
## misc.R

## commonly used calculations in Bayesian computation, including some tricks
## to stabilize calculations

# log_sum_exp():
# calculates expressions of the form log(sum(exp(x)))
log_sum_exp = function(x) { 
    offset = max(x)                         # scale by max to prevent overflow
    s = log(sum(exp(x - offset))) + offset
    i = which(!is.finite(s))                # check for any overflow
    if (length(i) > 0) {                    # replace inf values with max
        s[i] = offset 
    }
    
    return(s)
} # end of log_sum_exp()


# Wishart related calculations

## logB() function
# bishop appendix B.79
# log(B(W, nu)), B() is the normalizing constant of the Wishart distribution
logB = function(W, nu) {
    D = ncol(W)
    det_term = - 0.5 * nu * log(det(W))
    power_term = - (0.5 * nu * D * log(2) + 0.25 * D *(D - 1) * log(pi))
    prod_term = - sum(lgamma(0.5 * (nu + 1 - 1:D)))
    return(det_term + power_term + prod_term)
}







