

## various checks on the ELBO
## input:
#          VERBOSE      : logical, if TRUE, then progress printed each iter
#          i            : current iteration
#          max_iter     : max # of iterations
#          epsilon_conv : tolerance used to assess convergence
checkELBO = function(theta, prior) {
    
    CONVERGE = FALSE
    
    i = theta$curr
    
    # display iteration, ELBO, change in ELBO
    if (theta$VERBOSE) {
        cat("It:\t",i,"\tLB:\t",L[i], "\tLB_diff:\t",L[i] - L[i - 1],"\n")
    }
    
    # Check if lower bound decreases
    if (theta$L[i] < theta$L[i - 1]) { 
        message("Warning: Lower bound decreases!\n")
    }
    # Check for convergence
    if (abs(theta$L[i] - theta$L[i - 1]) < theta$epsilon_conv) { 
        CONVERGE = TRUE 
    }
    
    # Check if VB converged in the given maximum iterations
    if (i == theta$max_iter) {
        warning("VB did not converge!\n")
    }
    
    return(CONVERGE)
} # end of checkELBO() function


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


# printVector():
# print a K-dim vector in the form: [x1, x2, ..., xK]
vecToStr = function(x, K) {
    
    l     = '['
    r     = ']'
    x_str = paste(round(x, 4), collapse = ', ')
    
    return(paste(l, paste(x_str, collapse = ', '), r, sep = ' '))
} # end printVector()


# quadMult(): quadratic form multiplication, e.g., x'Ax
# input: 
#         x  : (D x 1) vector
#         A  : (D x D) matrix
# output: resulting product of x'Ax
quadMult = function(x, A) {
    return(t(x) %*% A %*% x)
}


# lambda_xi(): evaluate: 1 / (4 * theta$xi) * tanh(0.5 * theta$xi)
# input: 
#         xi : variational paramter
# output: 
#         1 / (4 * theta$xi) * tanh(0.5 * theta$xi)
lambda_xi = function(xi) {
    return(1 / (4 * xi) * tanh(0.5 * xi))
}





