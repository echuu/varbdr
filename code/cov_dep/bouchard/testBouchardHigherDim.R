
# testBouchardHigherDim.R

# testing bouchard bound for higher dimensions, i.e., K, d > 1
# K = number of clusters
# d = dimension of covariates

library(cubature)
options(digits = 6)

# testing Bouchard bound for K = 2, d = 1
testBound_2 = function() {
    K = 5
    d = 1
    maxIter = 200
    tol     = 1e-4
    bd      = rep(Inf, maxIter) # store the upper bound for each iter
    VERBOSE = FALSE
    
    ## initialize values for alpha, xi_1,...,xi_K
    alpha   = 0         # alpha can be any real number
    xi      = rep(0, K) # xi_k in [0, inf)
    xSigmax = rep(0, K)
    
    x = runif(1, min = 0, max = 1)                     # x ~ unif(0, 1)
    
    mu = rnorm(K, mean = 0, sd = 1)                    # mu ~ N(0, 1)
    sigma_sq = rgamma(n = K, shape = 1, rate = 1)      # sigma_sq ~ Ga(1, 1)
    
    f = function(beta) {
        q_beta = numeric(K)
        for (k in 1:K) {
            q_beta[k] =  dnorm(beta[k], mean = mu[k], sd = sqrt(sigma_sq[k]))
        }
        log_sum_exp(x * beta) * prod(q_beta)
    }
    
    # numerical integration value
    numer_val = adaptIntegrate(f, lowerLimit = rep(-Inf, K), 
                               upperLimit = rep(Inf, K))
    
    # mc approximation:
    R = 1e5
    # mc_approx = mcBound(x, K, R , mu, sigma_sq)
    # mc_approx 
    
    # precompute some quantities that are used in every iteration of the 
    # variational loop
    mu_x = t(x) %*% mu             # K x 1 vector, k-th element is mu_k'x
    for (k in 1:K) {
        xSigmax[k] = t(x) %*% sigma_sq[k] %*% x    # x' Sigma_k x
    }
    
    ## coordinate ascent to minimize the upper bound
    ## start at 2 so the first upper bound is guaranteed to be less than Inf
    for (i in 2:1000) {
        
        # update xi_1,...,xi_K
        xi = sqrt(xSigmax + mu_x^2 + alpha^2 - 2 * alpha * mu_x)  # K x 1 vector
        lambda_xi = 1 / (4 * xi) * tanh(0.5 * xi)                 # K x 1 vector 
        
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
            cat("Iter:  ", i, 
                #"\talpha:\t", alpha, 
                #"\txi:\t", xi, 
                "\tUB:  ", round(bd[i], 5),
                "\t Numerical:  ", numer_val$integral,
                "\t MC:  ", mc_approx,
                "\t\t Delta: ", bd[i] - numer_val$integral,
                "\n")
            break
        }
    } # end of coordinate ascent
    
    return(bd[i] >= numer_val$integral)
}

set.seed(1)
for (i in 1:10) {
    if (testBound_2()) {
        # print("bound holds: ")
    } else {
        cat("Iter: ", i, "\t upper bound below the ground truth")
    }
}




