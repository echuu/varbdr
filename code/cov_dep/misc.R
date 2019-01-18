

## various checks on the ELBO
## input:
# VERBOSE      : logical, if TRUE, then progress printed each iter
# i            : current iteration
# max_iter     : max # of iterations
# epsilon_conv : tolerance used to assess convergence
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
