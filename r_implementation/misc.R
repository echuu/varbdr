
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








