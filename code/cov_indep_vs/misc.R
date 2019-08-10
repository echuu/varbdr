
## misc.R

## misc calculations performed in CAVI for covariate-INDEPENDENT case


library(ggplot2)

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
    if (prior$VERBOSE) {
        cat("It:\t",        i,
            "\tLB:\t",      theta$L[i], 
            "\tLB_diff:\t", theta$L[i] - theta$L[i - 1],
            "\n")
    }
    
    # Check if lower bound decreases
    if (theta$L[i] < theta$L[i - 1]) { 
        message("Warning: Lower bound decreases!\n")
    }
    # Check for convergence
    if (abs(theta$L[i] - theta$L[i - 1]) < prior$tol) { 
        CONVERGE = TRUE 
    }
    
    # Check if VB converged in the given maximum iterations
    if (i == prior$max_iter) {
        warning("VB did not converge!\n")
    }
    
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



# computes the log of the normalizing constant of the Dirichlet density
# log(1/B(alpha)) = log(gamma(sum(alpha_k)) / prod(gamma(alpha_k)))
# intput:
#         alpha_k : (K x 1) or (1 x 1) if symm. Dirichlet; concentration paramt
#         K       : (1 x 1) dimension of dirichlet random variable
# output: 
#         log of the normalizing constant of Dirichlet density
log_dir_const = function(alpha_k, K) {
  
    log_const = 0
    
    if (K != length(alpha_k)) {      # symmetric dirichlet, alpha_k is scalar
        alpha_k = rep(alpha_k, K)    # in this case, turn it into k-dim vec
    }
    
    log_const = lgamma(sum(alpha_k)) - sum(lgamma(alpha_k))
      
    return(log_const)
    
} # end log_dir_const() function


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


# sigmoid(x) returns the sigmoid of the elements of x. The sigmoid
# function is also known as the logistic link function. It is the
# inverse of logit(x).
sigmoid <- function (x) {
    1/(1 + exp(-x))
}

# ----------------------------------------------------------------------
# logit(x) returns the logit of the elements of X. It is the inverse of
# sigmoid(x).
logit <- function (x) {
    log((x)/((1 - x)))
}





