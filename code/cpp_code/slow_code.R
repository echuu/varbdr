
# slow_code.R -- code to implement in C++ so that we can pass in params from R
#                and calculate results in C++ for perf gain

setwd("C:/Users/chuu/varbdr/code/cpp_code")

source("C:/Users/chuu/varbdr/code/cov_dep/misc.R")
source("debug_funcs.R")

library("Rcpp")
library("microbenchmark")
library("RcppEigen")
library("RcppArmadillo")
library("RcppParallel")

N = 100
K = 3
D = 2

set.seed(1)
X = matrix(rnorm(N * D), N, D)                   # (N x D)
mu = matrix(rnorm(D * K), D, K)                  # (D x K)

X_mu = X %*% mu                                  # (N x K)

xi = matrix(1, N, K)                             # (N x K)
lambda = lambda_xi(xi)                           # (N x K)

Qk = diag(rbeta(D, 1, 1))                        # (D x D)

alpha = numeric(N)                               # (N x 1)
xQx = numeric(K)                                 # (K x 1)


slow_func = function(lambda, X, X_mu, Qk, xi) {
    
    N = nrow(lambda)
    K = ncol(lambda)
    
    phi = numeric(N)
    
    for (n in 1:N) {
        alpha[n] = 1 / sum(lambda[n,]) * 
            (0.5 * (0.5 * K - 1) + crossprod(X_mu[n,], lambda[n,]))
        
        for (k in 1:K) {
            xQx[k] = sum(lambda[,k]) * crossprod(X[n,], Qk %*% X[n,])
        }
        
        
        xi[n,] = sqrt((X_mu[n,] - alpha[n])^2 + xQx)
        
    }
    
    lambda = lambda_xi(xi)
    
    
    # (0.4) compute phi (function of alpha, xi) --------------------------------
    for (n in 1:N) {
        phi[n] = sum((X_mu[n,] - alpha[n] - xi[n,]) / 2 + 
                               log(1 + exp(xi[n,])))
    } # end for() update for phi
    
    
    return(list(alpha = alpha, xQx = xQx, xi = xi, phi = phi, lambda = lambda))
}


sourceCpp("fast_functions.cpp")

tmp = testFuncs(X_mu, alpha, xQx) # used to test smaller functions


res = slow_func(lambda, X, X_mu, Qk, xi)        # slower version of the m-step
res_fast = mainFunc(lambda, X, X_mu, Qk, xi)    # build this to be the m-step

head(res$phi)
head(res_fast$phi)

checkEqual(res, res_fast)

microbenchmark(slow_func(lambda, X, X_mu, Qk, xi), 
               mainFunc(lambda, X, X_mu, Qk, xi))



## compare results to Rcpp -----------------------------------------------------

# install.packages("Rtools)



# (1) find skeleton for Rcpp Armadillo for general purpose vector/matrix mult
# (2) adapt code for the lines above
#     note: multiple returns
#           eventually will need list of matricies
#           scalar * matrix operations are ABUNDANT in R code



# m-step code ------------------------------------------------------------------

# for (n in 1:N) {
#    theta$alpha[n] = 1 / sum(theta$lambda[n,]) * 
#        (0.5 * (0.5 * K - 1) + crossprod(X_mu[n,], theta$lambda[n,]))

# (0.2) update xi (computed row-wise) ----------------------------------
#    xQx = numeric(K)
#    for (k in 1:K) {
#        xQx[k] = quadMult(X[n,], theta$Q_k_inv[,,k]) # function in misc.R
#    }

#    theta$xi[n,] = sqrt((X_mu[n,] - theta$alpha[n])^2 + xQx)
# }




