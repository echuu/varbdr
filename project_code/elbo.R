
## elbo.R
## calculate the variational lower bound using the current variational params

library(matrixcalc)

vb_elbo = function(theta, prior) {
    
    # initialize the 7 expectations to be computed
    e_ln_p_y = e_ln_p_z = e_ln_p_gamma = e_ln_p_beta_tau = 0
    e_ln_q_z = e_ln_q_beta_tau = e_ln_q_gamma = 0
    
    N = nrow(X)
    K = ncol(X)
    X = theta$X
    y = theta$y
    
    # commonly computed terms
    psi_a = digamma(theta$a_k)    # (K x 1)
    psi_b = digamma(theta$b_k)    # (K x 1)
    
    
    # INCOMPLETE: alpha, xi, phi, related quantities ---------------------------
    alpha = numeric(N)                # (N x 1)
    xi    = matrix(0, N, K)           # (N x K)
    M     = theta$mu_k %*% lambda     # (D x N) : (D x K) * (K x N)
    
    
    # recompute alpha
    for (n in 1:N) {
        alpha[n] = X[n,] %*% M[,n]
    }
    
    # recompute xi (computed column-wise)
    X_mu = X %*% mu    # (N x K) : (N x D) * (D x K)
    for (k in 1:K) {
        
        x_Qk_x = numeric(N)
        for (n in 1:N) {
            x_Qk_x[n] = X[n,] %*% theta$Q_k_inv[,,k] %*% t(X[n,])
        }
        
        xi[,k] = (X_mu[,k] - alpha)^2  + x_Qk_x  # n-dim col vector in k-th col
        
    }
    
    # compute phi
    
    
    # compute 7 expectations ---------------------------------------------------
    
    
    # E [ ln p(y | X, beta, tau, Z) ]
    e1 = numeric(N)
    for (n in 1:N) {
        x_Vinv_x = numeric(K)
        for (k in 1:K) {
            x_Vinv_x[k] = X[n,] %*% theta$V_inv[,,k] %*% t(X[n,]) # scalar
        } # end of inner for()
        
        # perform the inner k-summation
        e1[n] = sum(r[n,] * log(2 * pi) - psi_a - psi_b + x_Vinv_x + 
                        a_k / b_k * (y[n] - t(X[n,] %*% theta$m_k)))
        
    }
    e_ln_p_y = sum(e1)
    
    
    # E [ ln p(Z | X, gamma) ]
    e2 = numeric(K)
    for (k in 1:K) {
        # exchange the summations; perform the inner n-summation
        e2[k] = sum(r[,k] * (X_mu[,k] - alpha - phi))
    }
    e_ln_p_z = sum(e2)

    
    # E [ ln p(gamma) ]
    e_ln_p_gamma = - 0.5 * K * D * log(2 * pi) - 0.5 * sum(theta$mu_k^2)
    
    # E [ ln p(beta, tau) ]
    e3 = sum((prior$a_0 + 0.5 * D - 1) * (psi_a - psi_b)) -
        K * (0.5 * D * log(2 * pi) - log(det(prior$Lambda_0)) - 
                 prior$a_0 * log(prior$b_0) + lgamma(prior$a_0))
    
    diff = sweep(m_k, MARGIN = 1, STATS = prior$m_0, FUN = '-')    # (D x K)
    e3_k = numeric(K)
    for (k in 1:K) {
        e3_k[k] = t(diff[,k]) %*% prior$Lambda_0 %*% diff[,k] +
            matrix.trace(prior$Lambda_0 %*% theta$V_k_inv[,,k])
    }
    e3_k = theta$a_k / theta$b_k * e3_k
    
    e_ln_p_beta_tau = e3 + sum(e3_k)
    
    # E [ ln q(Z) ]
    e_ln_q_z = sum(theta$r * theta$log_r) # sum all elements of (N x K) matrix
    
    # E [ ln q(beta, tau) ]
    e6_k = (0.5 * D + a_k - 1) * (psi_a - psi_b) + 
        theta$a_k * log(theta$b_k) - theta$a_k - lgamma(theta$a_k)  # (K x 1)
    for (k in 1:K) {
        e6_k[k] = e6_k[k] + log(det(theta$V_k))
    }
    e_ln_q_beta_tau = sum(e6_k) - 0.5 * K * D * (log(2 * pi) + 1)
    
    
    # E [ ln q(gamma) ]
    e7 = - 0.5 * K * D * (log(2 * pi) + 1)
    for (k in 1:K) {
        e7 = e7 + log(det(theta$Q_k))
    }
    
    e_ln_q_gamma = e7
    
    # compute the elbo ---------------------------------------------------------
    vlb = e_ln_p_y + e_ln_p_z + e_ln_p_gamma + e_ln_p_beta_tau +
        e_ln_q_z + e_ln_q_gamma + e_ln_q_beta_tau
    
    return(vlb)
    
} # end vb_elbo() function


