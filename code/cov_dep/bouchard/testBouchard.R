


set.seed(1)



for (i in 1:100) {
    if (testBound()) {
        print("bound holds: ")
    } else {
        cat("Iter: ", i, "\t upper bound below the ground truth")
    }
}


testBound = function() {
    K = 1
    d = 1
    maxIter = 200
    tol     = 1e-4
    bd      = rep(Inf, maxIter) # store the upper bound for each iter
    VERBOSE = FALSE
    
    
    ## initialize values for alpha, xi_1,...,xi_K
    alpha   = 0         # alpha can be any real number
    xi      = rep(0, K) # xi_k in [0, inf)
    xSigmax = rep(0, K)
    
    x = runif(1, min = 0, max = 1)      # x ~ unif(0, 1)
    
    mu = rnorm(1, mean = 0, sd = 1)     # mu ~ N(0, 1)
    sigma_sq = rgamma(1, 1)             # sigma_sq ~ Ga(1, 1)
    
    true_val = x * mu
    
    # precompute some quantities that are used in every iteration of the 
    # variational loop
    mu_x = t(x) %*% mu             # K x 1 vector, k-th element is mu_k'x
    for (k in 1:K) {
        xSigmax[k] = x^2 * sigma_sq
    }
    
    ## coordinate ascent to minimize the upper bound
    ## start at 2 so the first upper bound is guaranteed to be less than Inf
    for (i in 2:1000) {
        
        # update xi_1,...,xi_K
        xi = sqrt(xSigmax + mu_x^2 + alpha^2 - 2 * alpha * mu_x)   # K x 1 vector
        lambda_xi = 1 / (4 * xi) * tanh(0.5 * xi)                  # K x 1 vector 
        
        # update alpha
        alpha = (0.5 * (0.5 * K - 1) + sum(lambda_xi * mu_x)) / sum(lambda_xi)
        
        
        # compute bound, store in bd[i]
        bd[i] = computeBound(mu_x, xSigmax, xi, lambda_xi, alpha, K)
        
        if (VERBOSE) {
            cat("Iter:  ", i, 
                "\talpha:\t", alpha, 
                "\txi:\t", xi, 
                "\tUB:  ", bd[i],
                "\tLB_diff:  ", bd[i] - bd[i - 1],
                "\n")
        }
        
        # check convergence:
        if (checkConvergence(bd, i, maxIter, tol)) {
            break
        }
    } # end of coordinate ascent
    
    return(bd[i] >= true_val)
}

