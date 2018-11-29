
## bouchard.R

library(MASS)

## testing the bound for E[logSumExp] proposed in Bouchard (2008)


D   = 2          # dimension of covariate, 2 <= d <= 100
K   = 2          # number of clusters
N   = 1          # number of samples
I_D = diag(1, D) # D X D  identity matrix
df  = 100        # degrees of freedom for Wishart distribution

# generate a random sample x in R^D (from any distribution)
mu0    = c(rnorm(D, 3, 4))
Sigma0 = 3 * I_D 
x = mvrnorm(N, mu0, Sigma0)

# sample mu_k ~ N(0, I_D)ye
mu = t(mvrnorm(K, rep(0, D), I_D)) # D x K : mu_k stored column wise 

# sample Sigma_k ~ W(df, 10 * I_D)
Sigma = rWishart(K, df, 10 * I_D) # array of K matrices; Sigma[,,k]

maxIter = 200
bd = rep(Inf, maxIter) # store the upper bound for each iter, should decrease

## coordinate ascent to minimize the upper bound
## initialize values for alpha, xi_1,...,xi_K
## start at 2 so the first upper bound is guaranteed to be less than Inf
for (i in 2:maxIter) {
    
    # update xi_1,...,xi_K
    
    # update alpha
    
    # compute bound, store in bd[i]
    
    # check convergence (difference in bound, difference in parameters)
    
} # end of coordinate ascent



# mu    : D X K        matrix of mean vectors stored columnwise
# Sigma : K x (D x D)  array of (K x K) covariance matrices
# alpha : 1 x 1        variational parameter used in the bound
# xi    : K x 1        variational parameter used the bound
# K     : 1 x 1        number of clusters
computeBound = function(mu, Sigma, alpha, xi, K) {
    const     = - alpha * (0.5 * K - 1) # 1 x 1 constant term in the bound  
    lambda_xi = tanh(xi)                # K x 1 vector 
    mu_x      = t(x) %*% mu             # K x 1 vector, k-th element is mu_k'x
    for (k in 1:K) {
        xSigmax[k] = t(x) %*% Sigma[,,k] %*% x
    }
    bound = sum(lambda_xi * (xSigmax + mu_x^2) + 
                    (0.5 - 2 * alpha * lambda_xi) * mu_x -
                    0.5 * xi + lambda_xi * (alpha^2 - xi^2) + log(1 + exp(xi)))
    bound = bound + const
    return(bound)
}

xSigmax = numeric(K)








