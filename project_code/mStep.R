
## mStep.R
## perform the variational e-step

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
    
    # alpha, xi, phi, related quantities ---------------------------------------
    # theta$alpha                 # (N x 1)
    # theta$ xi                   # (N x K)
    # theta$phi                   # (N x 1)
    
    # update alpha
    for (n in 1:N) {
        theta$alpha[n] = (0.5 * (0.5 * K - 1) + 
                              X[n,] %*% M[,n]) / sum(theta$lambda[n,])
    }
    
    # update xi (computed column-wise)
    for (k in 1:K) {
        
        x_Qk_inv_x = numeric(N)
        for (n in 1:N) {
            x_Qk_inv_x[n] = X[n,] %*% theta$Q_k_inv[,,k] %*% t(X[n,])
        }
        
        # populate the k-th column of xi
        xi[,k] = (X_mu[,k] - theta$alpha)^2  + x_Qk_inv_x  
    }
    
    # compute lambda(xi) using updated value of xi
    theta$lambda = 1 / (4 * theta$xi) * tanh(0.5 * theta$xi)    # (N x K)
    
    # compute phi (function of alpha, xi)
    for (n in 1:N) {
        theta$phi[n] = sum((X_mu[n,] - theta$alpha[n] - theta$xi[n,]) / 2 + 
                               log(1 + exp(theta$xi[n,])))
    }
    
    
    # update variational distributions -----------------------------------------
    
    # update q(gamma) : Q_k, Q_k^{-1}, eta_k, mu_l
    for (k in 1:K) {
        theta$Q_k[,,k] = I_D + 
            2 * t(X) %*% diag(theta$r_nk[,k] * lambda[,k]) %*% X
        theta$Q_k_inv[,,k] = solve(theta$Q_k[,,k])
        theta$eta_k[,k] = t(X) %*% 
            c(theta$r_nk[,k] * (0.5 + 2 * theta$lambda[,k] * theta$alpha))
        theta$mu_k[,k] =  theta$Q_k_inv[,,k] %*% theta$eta_k[,k]
    }
    
    # update q(beta | tau) : V_k, V_k^{-1}, zeta_k, m_k
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
            theta$zeta_k[,k] + sum(r_nk[,k]* y^2)
    }
    
    theta$b_k = (b_k + t(prior$m_0) %*% prior$Lambda_0 %*% prior$m_0) / 2 + b_0
    
    return(theta)
    
} # end mStep()


# end of mStep.R


