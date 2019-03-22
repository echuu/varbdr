
# slow_code.R -- code to implement in C++ so that we can pass in params from R
#                and calculate results in C++ for perf gain

source("C:/Users/chuu/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)


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
D = 1

synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = synth_data_1d$X
intercept = FALSE
max_iter = 5000

mu = matrix(rnorm(D * K), D, K)                  # (D x K)

X_mu = X %*% mu                                  # (N x K)

xi = matrix(1, N, K)                             # (N x K)
lambda = lambda_xi(xi)                           # (N x K)

Qk = diag(rbeta(D, 1, 1))                        # (D x D)
Qk = rbeta(D, 1, 1)

alpha = numeric(N)                               # (N x 1)
xQx = numeric(K)                                 # (K x 1)

r_nk = matrix(runif(N * K), N, K)                # (N x K)
r_nk = r_nk / rowSums(r_nk)

N_k = colSums(r_nk);                             # (K x 1)

Lambda_0 = diag(rep(1, ncol(X)))
m_0 = c(colMeans(X))
Lambda0_m0 = Lambda_0 %*% m_0
m0_Lambda0_m0 = c(quadMult(m_0,  Lambda_0))

a_0 = 1
b_0 = 1

slow_func = function(lambda, X, X_mu, Qk, xi) {
    
    N = nrow(lambda)
    K = ncol(lambda)
    
    alpha = numeric(N)
    phi = numeric(N)
    xQx = numeric(K)
    
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

slow_update = function(X, y, r_nk, N_k, lambda, alpha, 
                       Lambda_0, Lambda0_m0, m0_Lambda0_m0, 
                       a_0, b_0) {
    
    N = nrow(X)
    D = ncol(X)
    K = ncol(lambda)
    
    I_D = diag(1, D)
    
    a_k = numeric(K)
    b_k = numeric(K)
    
    
    for (k in 1:K) {
        
        rl_nk_xx = matrix(0, D, D)               # used to calculate q(gamma)
        r_nk_xx =  matrix(0, D, D)               # used to calculate q(beta|tau)
        
        for (n in 1:N) {
            
            r_x = r_nk[n,k] * crossprod(t(X[n,]))
            
            rl_nk_xx = rl_nk_xx + r_x * lambda[n,k]
            
            r_nk_xx = r_nk_xx + r_x
        }
        
        
        # (1.1) update q(gamma) : Q_k, Q_k^{-1}, eta_k, mu_k ------------------
        Q_k     = I_D + 2 * rl_nk_xx
        Q_k_inv = solve(Q_k)
        eta_k    = t(X) %*% (r_nk[,k] * (0.5 + 2 * lambda[,k] * alpha))
        mu_k     =  Q_k_inv %*% eta_k
        
        # (1.2) update q(beta | tau) : V_k, V_k^{-1}, zeta_k, m_k --------------
        V_k     = as.matrix(Lambda_0 + r_nk_xx)
        V_k_inv = solve(V_k)
        zeta_k   = Lambda0_m0 + crossprod(X, (r_nk[,k] * y))
        m_k      = V_k_inv %*% zeta_k
        
    } # end of q(gamma), q(beta|tau) updates
    
    
    # (1.3) update q(tau): a_k, b_k --------------------------------------------
    a_k = a_0 + 0.5 * N_k # (K x 1) update
    for (k in 1:K) {
        b_k[k] = - quadMult(zeta_k, V_k_inv) + 
            sum(r_nk[,k] * (y^2))
    }
    b_k = b_0 + 0.5 * (b_k + m0_Lambda0_m0)
    
    
    # update the random variables using posterior means ------------------------
    
    # beta_k  : coefficient vector in the normal density
    beta_k = m_k                # (D x K)
    
    # tau_k   : precision parameter in the normal density
    tau_k = a_k / b_k     # (K x 1)
    
    # gamma_k : coefficient vector in the mixture weights
    gamma_k = mu_k              # (D x K)
    
    
    
    return(list(r_x = r_x, rl_nk_xx = rl_nk_xx, r_nk_xx = r_nk_xx,
                Q_k = Q_k, Q_k_inv = Q_k_inv, eta_k = eta_k, mu_k = mu_k,
                V_k = V_k, V_k_inv = V_k_inv, zeta_k = zeta_k, m_k = m_k,
                a_k = a_k, b_k = b_k,
                beta_k = beta_k, tau_k = tau_k, gamma_k = gamma_k))
    
} #  end slow_update() function


# ------------------------------------------------------------------------------


source("C:/Users/chuu/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)


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
D = 1

synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = synth_data_1d$X
intercept = FALSE
max_iter = 5000

sourceCpp("getVarParams.cpp")
theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)
str(theta_cpp)

theta_cpp$m_0           # should match colMeans(X)

colMeans(X)
head(theta_cpp$m_k)
head(theta_cpp$mu_k)
head(theta_cpp$lambda)

mat_list = mat_list_ops(K, 2, rep(1, 2))
str(mat_list)

# ------------------------------------------------------------------------------

sourceCpp("test_funcs.cpp")


theta_cpp = testConstructor(y, X, N, D, K, intercept, max_iter)
str(theta_cpp)

theta_cpp$L

testUpdate = slow_update(X, y, r_nk, N_k, lambda, alpha, 
                         Lambda_0, Lambda0_m0, m0_Lambda0_m0, 
                         a_0, b_0)

testUpdate$b_k
testUpdate$a_k
testUpdate$b_k
testUpdate$V_k_inv
testUpdate$tau_k
testUpdate$zeta_k

fastUpdate = updateVariational(X, y, r_nk, N_k, lambda, alpha, 
                               Lambda_0, Lambda0_m0, m0_Lambda0_m0, 
                               a_0, b_0)

fastUpdate$b_k


fastUpdate$a_k
fastUpdate$zeta_k
fastUpdate$V_k_inv
fastUpdate$gamma_k




checkEqual(testUpdate, fastUpdate)

microbenchmark(slow_update(X, y, r_nk, N_k, lambda, alpha, 
                           Lambda_0, Lambda0_m0, m0_Lambda0_m0, 
                           a_0, b_0), 
               updateVariational(X, y, r_nk, N_k, lambda, alpha, 
                                 Lambda_0, Lambda0_m0, m0_Lambda0_m0, 
                                 a_0, b_0))

# ------------------------------------------------------------------------------


sourceCpp("matrix_ops.cpp")


mat_list_ops(1, 3, c(1, 2, 3))





