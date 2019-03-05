
# slow_code.R -- code to implement in C++ so that we can pass in params from R
#                and calculate results in C++ for perf gain

setwd("C:/Users/chuu/varbdr/code/cpp_code")

library("Rcpp")
library("microbenchmark")
library("RcppEigen")
library("RcppArmadillo")
library("RcppParallel")

N = 10
K = 3
D = 2

X = matrix(rnorm(N * D), N, D)
lambda = matrix(rnorm(N * K), N, K)

Qk = diag(rnorm(D))

alpha = numeric(N)
xQx = numeric(K)

for (n in 1:N) {
    
    alpha[n] = sum(lambda[n,]) * crossprod(X[n,])
    
    for (k in 1:K) {
        # matrix product that involves alpha_n
        xQx[k] = sum(lambda[,k]) * crossprod(X[n,], Qk %*% X[n,])
    }
    
}
   
slow_func = function(lambda, X, Qk) {
    
    N = nrow(lambda)
    K = ncol(lambda)
    
    for (n in 1:N) {
        alpha[n] = sum(lambda[n,]) * crossprod(X[n,])
        
        for (k in 1:K) {
            xQx[k] = sum(lambda[,k]) * crossprod(X[n,], Qk %*% X[n,])
        }
        
    }
    
    return(list(alpha = alpha, xQx = xQx))
}

res = slow_func(lambda, X, Qk)

res$alpha
res$xQx


sourceCpp("fast_functions.cpp")

Qk[1,] %*% Qk %*% Qk[1,]
multC(Qk)

res_fast = mainFunc(lambda, X, Qk)
res_fast$l_sum
res_fast$xQx

slow_func(lambda, X) == mainFunc(lambda, X, Qk)


microbenchmark(slow_func(lambda, X, Qk), mainFunc(lambda, X, Qk))


## compare results to Rcpp -----------------------------------------------------

# install.packages("Rtools)



# (1) find skeleton for Rcpp Armadillo for general purpose vector/matrix mult
# (2) adapt code for the lines above
#     note: multiple returns
#           eventually will need list of matricies
#           scalar * matrix operations are ABUNDANT in R code


matmultC(t(X), X)

set.seed(1)
A = matrix(rnorm(10000), 100, 100)
y = rnorm(100)

t(y) %*% A %*% y

crossprod(y, A %*% y)

microbenchmark(t(y) %*% A %*% y, 
               crossprod(y, A %*% y),
               quadMult(A, y))


microbenchmark(schur(y, y), y*y)


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




