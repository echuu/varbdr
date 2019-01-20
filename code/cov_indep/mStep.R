
## mStep.R
## perform the variational m-step -- covariate-INDEPENDENT case

library(matrixcalc)

# input: 
#          theta : list of variational parameters
# output:
#          theta : list of variational parameters with 
#                  variational parameters updated

mStep() = function(theta, prior) {
    
    I_D        = diag(1, D)                      # (D X D) : identity matrix
    X_mu       = X %*% mu                        # (N x K) : (N x D) * (D x K)
    M          = theta$mu_k %*% t(theta$lambda)  # (D x N) : (D x K) * (K x N)
    Lambda0_m0 = prior$Lambda_0 %*% prior$m_0    # (D x 1) : Lambda_0 * m_0
    
    X = prior$X
    y = prior$y
    N = nrow(X)
    K = ncol(X)
    
    
    # update q(pi) = Dir( pi | alpha_k )
    #    parameter to update: alpha_k
    theta$alpha_k = prior$alpha_0 + theta$N_k
    
    
    # update q(beta_k | tau_k) = N ( beta_k | m_k, (tau_k V_k)^{-1} ) 
    #     parameters to update: V_k, V_k^{-1}, zeta_k, m_k
    for (k in 1:K) {
        theta$V_k[,,k] = prior$Lambda_0 + t(X) %*% diag(theta$r_nk[,k]) %*% X
        theta$V_k_inv[,,k] = solve(theta$V_k[,,k])
        theta$zeta_k[,k] = Lambda0_m0 + t(X) %*% (theta$r_nk[,k] * y)
        theta$m_k[,k] = theta$V_k_inv[,,k] %*% theta$zeta_k[,k]
    }
    
    
    # update q(tau): a_k, b_k
    theta$a_k = prior$a_0 + theta$N_k
    for (k in 1:K) {
        theta$b_k[k] = - t(theta$zeta_k[,k]) %*% theta$V_k_inv[,,k] %*% 
            theta$zeta_k[,k] + sum(theta$r_nk[,k] * y^2)
    }
    
    theta$b_k = (b_k + t(prior$m_0) %*% prior$Lambda_0 %*% prior$m_0) / 2 + b_0
    
    # update the current iteration
    theta$curr = theta$curr + 1
    
    return(theta)
    
} # end mStep()


# end of mStep.R



