
## bouchard.R

library(MASS)

## testing the bound for E[logSumExp] proposed in Bouchard (2008)


D     = 2          # dimension of covariate, 2 <= d <= 100
K     = 2          # number of clusters
N     = 1          # number of samples
I_D   = diag(1, D) # D X D  identity matrix
df    = 100        # degrees of freedom for Wishart distribution
scale = 1          # scale for the covariance matrix


# generate a random sample x in R^D (from any distribution)
mu0    = c(rnorm(D, 3, 4))
Sigma0 = 3 * I_D 
x      = mvrnorm(N, mu0, Sigma0)

# sample mu_k ~ N(0, I_D)ye
mu = t(mvrnorm(K, rep(0, D), I_D)) # D x K : mu_k stored column wise 

# sample Sigma_k ~ W(df, 10 * I_D)
Sigma = rWishart(K, df, scale * I_D) # array of K matrices; Sigma[,,k]

# compute the mc-approximate value -- this will be the baseline, what we compare
# the variational bound to
R = 1e5
mc_approx = mcBound(x, K, D, R , mu, Sigma)
mc_approx # 41.65944

# --------------------- begin variational approximation ---------------------- #

maxIter = 200
bd = rep(Inf, maxIter) # store the upper bound for each iter, should decrease

## initialize values for alpha, xi_1,...,xi_K
alpha = 0         # alpha can be any real number
xi    = rep(0, K) # xi_k in [0, inf)

## coordinate ascent to minimize the upper bound
## start at 2 so the first upper bound is guaranteed to be less than Inf
for (i in 2:maxIter) {
    
    # update xi_1,...,xi_K
    mu_x      = t(x) %*% mu             # K x 1 vector, k-th element is mu_k'x
    for (k in 1:K) {
        xSigmax[k] = t(x) %*% Sigma[,,k] %*% x
    }
    xi = xSigmax + mu_x^2 + alpha^2 - 2 * alpha * mu_x
    
    # update alpha
    lambda_xi = 1 / (4 * xi) * tanh(0.5 * xi)            # K x 1 vector 
    alpha = (0.5 * (0.5 * K - 1) + sum(lambda_xi * mu_x)) / sum(lambda_xi)
    
    # compute bound, store in bd[i]
    bd[i] = computeBound(mu_x, xSigmax, lambda_xi, xi, alpha, K)
    # check convergence (difference in bound, difference in parameters)
    
} # end of coordinate ascent

# ------------------------- helper functions --------------------------------- #

# mu_x      : D X K        matrix of mean vectors stored columnwise
# xSigmax   : K x (D x D)  array of (K x K) covariance matrices
# xi        : K x 1        variational parameter used the bound
# lambda_xi : K x 1        1 / (4 * xi) * tanh(0.5 * xi)
# alpha     : 1 x 1        variational parameter used in the bound
# K         : 1 x 1        number of clusters
computeBound = function(mu_x, xSigmax, xi, lambda_xi, alpha, K) {
    const     = - alpha * (0.5 * K - 1) # 1 x 1 constant term in the bound  
    # lambda_xi = tanh(xi)                # K x 1 vector 
    # mu_x      = t(x) %*% mu             # K x 1 vector, k-th element is mu_k'x
    # xSigmax   = rep(0, K)
    # for (k in 1:K) {
    #     xSigmax[k] = t(x) %*% Sigma[,,k] %*% x
    # }
    
    bd = sum(lambda_xi * (xSigmax + mu_x^2) + 
                 (0.5 - 2 * alpha * lambda_xi) * mu_x -
                 0.5 * xi + lambda_xi * (alpha^2 - xi^2) + log(1 + exp(xi)))
    bound = bound + const
    return(bound)
} # end of computeBound() function


# mcBound() ---- computes a monte carlo estimate of E[log(sum(exp(x'beta)))] 
# input:
    # x      : D x 1        randomly generated vector
        # should x be the same for all R iterations? I would think so, since we
        # want the MC estimate to be a function of the betas
    # K      : 1 x 1        number of clusters
    # D      : 1 x 1        dimension of covariate
    # R      : 1 x 1        number of MC samples
    # mu     : D x K        mean vector for each of the K clusters
    # Sigma  : K x (D x D)  K-dim array of (D x D) cov matrices for each cluster
mcBound = function(x, K, D, R = 1e5, mu, Sigma, seed = 1) {
    set.seed(seed)
    mcSum = 0;
    b_k   =  matrix(0, D, K) # D x K matrix, coefficient vectors stored col-wise
    # outer MC loop
    for (r in 1:R) {
        
        for (k in 1:K) {
            # sample beta_k ~ N (mu_k, Sigma_k)
            b_k[,k] = mvrnorm(1, mu[,k], Sigma[,,k]) 
        }
        
        # mcSum = mcSum + log(sum(exp(t(x) %*% b_k)))     
        mcSum = mcSum + log_sum_exp(t(x) %*% b_k)    # more stable calculation
    } # end of outer MC loop
    
    return(1 / R * mcSum)
    
} # end of mcBound() function


# log_sum_exp():
# calculates expressions of the form log(sum(exp(x)))
log_sum_exp = function(x) { 
    offset = max(x)                         
    s = log(sum(exp(x - offset))) + offset  # scale by max to prevent overflow
    i = which(!is.finite(s))                # check for any overflow
    if (length(i) > 0) {                    # replace inf values with max
        s[i] = offset 
    }
    return(s)
} # end of log_sum_exp()



# end of bouchard.R
