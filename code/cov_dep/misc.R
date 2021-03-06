
## misc.R -- helper functions

## included functions:
##     checkELBO       : various checks on the ELBO (convergence, etc.)
##     plotELBO()      : plot the ELBO over the iterations of CAVI
##     log_sum_exp()   : log(sum(exp(x)))
##     vecToStr()      : print a K-dim vector in the form: [x1, x2, ..., xK]
##     quadMult()      : quadratic form multiplication, e.g., x'Ax
##     lambda_xi()     : evaluate: 1 / (4 * theta$xi) * tanh(0.5 * theta$xi)


library(ggplot2)


## checkELBO() : various checks on the ELBO
## input:
#          VERBOSE      : logical, if TRUE, then progress printed each iter
#          i            : current iteration
#          max_iter     : max # of iterations
#          epsilon_conv : tolerance used to assess convergence
checkELBO = function(theta, prior) {
    
    CONVERGE = FALSE
    
    i = theta$curr
    
    # display iteration, ELBO, change in ELBO
    if ((prior$VERBOSE) && ((i - 1) %% 20 == 0)) {
        cat("It:\t",        i - 1,
            "\tLB:\t",      theta$L[i], 
            "\tLB_diff:\t", theta$L[i] - theta$L[i - 1],
            "\n")
    }
    
    # Check if lower bound decreases
    # if (theta$L[i] < theta$L[i - 1]) { 
    #     message("Warning: Lower bound decreases!\n")
    # }
    # Check for convergence
    if (abs(theta$L[i] - theta$L[i - 1]) < prior$tol) { 
        CONVERGE = TRUE 
    }
    
    # Check if VB converged in the given maximum iterations
    # if (i == prior$max_iter) {
    #     warning("VB did not converge!\n")
    # }
    
    return(CONVERGE)
} # end of checkELBO() function


## plotELBO() : plot the ELBO over the iterations
## input:
#           theta      : object containing variational parameters
## output: 
#           ggplot plot/object with the ELBO value at each iter until converge
plotELBO = function(theta) {
    
    elbo_df = data.frame(iter = 2:theta$curr, elbo = theta$L[2:theta$curr])
    
    elbo_plot = ggplot(elbo_df, aes(x = iter, y = elbo)) + 
        geom_point(size = 0.9)
    
    return(elbo_plot)
}



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
    return(t(x) %*% (A %*% x))
    # return(tcrossprod(x, A %*% x))
}


# lambda_xi(): evaluate: 1 / (4 * theta$xi) * tanh(0.5 * theta$xi)
# input: 
#         xi : variational paramter
# output: 
#         1 / (4 * theta$xi) * tanh(0.5 * theta$xi)
lambda_xi = function(xi) {
    return(1 / (4 * xi) * tanh(0.5 * xi))
}





