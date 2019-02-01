
## mStep.R
## perform the variational e-step

source("misc.R")

library(matrixcalc)

# input: 
#          theta : list of variational parameters
# output:
#          theta : list of variational parameters with 
#                  variational parameters updated

mStep = function(theta, prior) {
    
    
    X = prior$X
    y = prior$y
    N = prior$N
    K = prior$K
    D = prior$D

    I_D        = diag(1, D)                       # (D X D) : identity matrix
    X_mu       = X %*% theta$mu_k                 # (N x K) : (N x D) * (D x K)
    # M          = theta$mu_k %*% t(theta$lambda)   # (D x N) : (D x K) * (K x N)
    Lambda0_m0 = prior$Lambda_0 %*% prior$m_0     # (D x 1) : Lambda_0 * m_0
    
    ## update alpha, xi, lambda, phi --------------------------------------------
    #     (0.1) alpha                 # (N x 1)
    #     (0.2) xi                    # (N x K)
    #     (0.3) lambda                # (N x K)
    #     (0.4) phi                   # (N x 1)
    
    # (0.1) update alpha -------------------------------------------------------------
    for (n in 1:N) {
        theta$alpha[n] = 1 / sum(theta$lambda[n,]) * 
            (0.5 * (0.5 * K - 1) + t(X_mu[n,]) %*% theta$lambda[n,])
    }
    
    
    # (0.2) update xi (computed row-wise) --------------------------------------------
    for (n in 1:N) {
        xQx = numeric(K)
        for (k in 1:K) {
            xQx[k] = quadMult(X[n,], theta$Q_k_inv[,,k]) # function in misc.R
        }
        
        theta$xi[n,] = sqrt((X_mu[n,] - theta$alpha[n])^2 + xQx)
    } # end of outer for() update for xi
    
    
    # (0.3) compute lambda(xi) using updated value of xi  ----------------------
    theta$lambda = lambda_xi(theta$xi)              # (N x K), misc.R
    
    
    # (0.4) compute phi (function of alpha, xi) --------------------------------
    for (n in 1:N) {
        theta$phi[n] = sum((X_mu[n,] - theta$alpha[n] - theta$xi[n,]) / 2 + 
                               log(1 + exp(theta$xi[n,])))
    } # end for() update for phi
    
    
    ## update variational distributions ----------------------------------------
    #     (1.1) q(gamma) 
    #     (1.2) q(beta|tau)
    #     (1.3) q(tau)
    
    # (1.1) update q(gamma) : Q_k, Q_k^{-1}, eta_k, mu_k -----------------------
    for (k in 1:K) {
        
        rl_nk_xx = matrix(0, D, D)
        for (n in 1:N) {
            rl_nk_xx = rl_nk_xx + 
                (theta$r_nk[n,k] * theta$lambda[n,k]) * crossprod(t(X[n,]))
        }
        
        theta$Q_k[,,k] = I_D + rl_nk_xx
        
        theta$Q_k_inv[,,k] = solve(theta$Q_k[,,k])
        
        theta$eta_k[,k] = t(X) %*% 
            (theta$r_nk[,k] * (0.5 + 2 * theta$lambda[,k] * theta$alpha))
        
        theta$mu_k[,k] =  theta$Q_k_inv[,,k] %*% theta$eta_k[,k]
    }
    
    # (1.2) update q(beta | tau) : V_k, V_k^{-1}, zeta_k, m_k ------------------
    for (k in 1:K) {
        
        r_nk_xx = matrix(0, nrow = D, ncol = D)
        # r_nk_yx = 0
        for (n in 1:N) {
            r_nk_xx = r_nk_xx + theta$r_nk[n,k] * crossprod(t(X[n,]))
            # r_nk_yx = r_nk_yx + theta$r_nk[n,k] * y[n] * X[n,]
        }
        
        # theta$V_k[,,k] = prior$Lambda_0 + t(X) %*% diag(theta$r_nk[,k]) %*% X
        
        theta$V_k[,,k]     = prior$Lambda_0 + r_nk_xx
        theta$V_k_inv[,,k] = solve(theta$V_k[,,k])
        theta$zeta_k[,k] = Lambda0_m0 + t(X) %*% (theta$r_nk[,k] * y)
        theta$m_k[,k] = theta$V_k_inv[,,k] %*% theta$zeta_k[,k]
    }
    
    # (1.3) update q(tau): a_k, b_k --------------------------------------------
    theta$a_k = prior$a_0 + 0.5 * theta$N_k
    for (k in 1:K) {
        # theta$b_k[k] = - t(theta$zeta_k[,k]) %*% theta$V_k_inv[,,k] %*% 
        #    theta$zeta_k[,k] + sum(theta$r_nk[,k] * y^2)
        theta$b_k[k] = - quadMult(theta$zeta_k[,k], theta$V_k_inv[,,k]) + 
            sum(theta$r_nk[,k] * y^2)
    }
    theta$b_k = prior$b_0 + 
        0.5 * (theta$b_k + c(quadMult(prior$m_0,  prior$Lambda_0)))
    
    
    # update the random variables using posterior means ------------------------
    
    # beta_k  : coefficient vector in the normal density
    theta$beta_k = theta$m_k                # (D x K)
    
    # tau_k   : precision parameter in the normal density
    theta$tau_k = theta$a_k / theta$b_k     # (K x 1)
    
    # gamma_k : coefficient vector in the mixture weights
    theta$gamma_k = theta$mu_k              # (D x K)
    
    
    # pi_k    : mixing weights (function of the gamma_k's)
    # maybe don't compute them here because we only need them when we compute
    # the density after CAVI has converged. If we compute here, then that's
    # wasted computation
    
    # update the current iteration
    theta$curr = theta$curr + 1
    
    
    return(theta)
    
} # end mStep()


# end of mStep.R


