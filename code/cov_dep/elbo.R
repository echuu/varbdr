
## elbo.R -- covariate-DEPENDENT case
## calculate the variational lower bound using the current variational params

source(paste(COV_DEP, MISC_FUNCS, sep = '/'))

elbo = function(theta, prior) {
    
    # initialize the 7 expectations to be computed
    e_ln_p_y = e_ln_p_z = e_ln_p_gamma = e_ln_p_beta_tau = 0
    e_ln_q_z = e_ln_q_beta_tau = e_ln_q_gamma = 0
    
    X = prior$X
    y = prior$y
    N = prior$N
    D = prior$D
    K = prior$K
    
    
    # commonly computed terms
    psi_a = digamma(theta$a_k)    # (K x 1)
    psi_b = digamma(theta$b_k)    # (K x 1)
    ak_bk = theta$a_k / theta$b_k # (K x 1) : a_k / b_k, elementwise division
    X_mu  = X %*% theta$mu_k      # (N x K) : [Xmu_1, Xmu_2, ... , Xmu_K]
    
    
    # compute 7 expectations ---------------------------------------------------
    
    # (1) E [ ln p(y | X, beta, tau, Z) ] --------------------------------------
    e1 = numeric(N)
    for (n in 1:N) {
        x_Vinv_x = numeric(K)
        for (k in 1:K) {
            x_Vinv_x[k] = quadMult(X[n,], theta$V_k_inv[,,k]) # x_n'V_k^{-1}x_n
        } # end of inner for()
        
        # perform the inner k-summation
        e1[n] = sum(theta$r_nk[n,] * 
                        (log(2 * pi) - psi_a + psi_b + x_Vinv_x + 
                             ak_bk * (y[n] - t(X[n,] %*% theta$m_k))^2))
        
    }
    
    
    e_ln_p_y = -0.5 * sum(e1)
    theta$e_ln_p_y = -0.5 * sum(e1)  # comment out later
    
    # (2) E [ ln p(Z | X, gamma) ] ---------------------------------------------
    e2 = numeric(K)
    for (k in 1:K) {
        # exchange the summations; perform the inner n-summation
        e2[k] = sum(theta$r_nk[,k] * (X_mu[,k] - theta$alpha - theta$phi))
    }
    e_ln_p_z = sum(e2) # outer k-summation
    theta$e_ln_p_z = sum(e2) # comment out later
    
    # (3) E [ ln p(gamma) ] ----------------------------------------------------
    trace_qk = 0
    for (k in 1:K) {
        trace_qk = trace_qk + matrix.trace(as.matrix(theta$Q_k_inv[,,k]))
    }
    e_ln_p_gamma = - 0.5 * K * D * log(2 * pi) - 0.5 * sum(theta$mu_k^2) - 
        0.5 * trace_qk
    theta$e_ln_p_gamma = - 0.5 * K * D * log(2 * pi) - 0.5 * sum(theta$mu_k^2)

    # (4) E [ ln p(beta, tau) ] ------------------------------------------------
    e3 = sum((prior$a_0 + 0.5 * D - 1) * (psi_a - psi_b)) -
        K * (0.5 * D * log(2 * pi) - 0.5 * log(det(prior$Lambda_0)) - 
                 prior$a_0 * log(prior$b_0) + lgamma(prior$a_0))
    
    # print(paste("e3 constant =", e3))
    
    # subtract each column of m_k by m_0 (centering each mean component)
    diff = sweep(theta$m_k, MARGIN = 1, STATS = prior$m_0, FUN = '-')  # (D x K)
    e3_diff_k  = numeric(K)  # store (m_k - m_0)' Lambda_0 (m_k - m_0)
    e3_trace_k = numeric(K)  # store trace(Lambda_0 V_K^{-1})
    for (k in 1:K) {
        e3_diff_k[k]  = quadMult(diff[,k], prior$Lambda_0)
        e3_trace_k[k] = matrix.trace(prior$Lambda_0 %*% theta$V_k_inv[,,k])
    }
    e3_diff_k = ak_bk * (e3_diff_k + prior$b_0) # element-wise multiplication
    e_ln_p_beta_tau = e3 - 0.5 * sum(e3_diff_k + e3_trace_k)
    theta$e_ln_p_beta_tau = e3 - 0.5 * sum(e3_diff_k + e3_trace_k) # comment out later
    
    
    # (5) E [ ln q(Z) ] -------------------------------------------------------- 
    e_ln_q_z = sum(theta$r_nk * theta$log_r_nk) # sum over (N x K) matrix
    theta$e_ln_q_z = sum(theta$r_nk * theta$log_r_nk) # comment out later
    
    
    # (6) E [ ln q(beta, tau) ] ------------------------------------------------
    e6_k = (0.5 * D + theta$a_k - 1) * (psi_a - psi_b) + 
        theta$a_k * log(theta$b_k) - theta$a_k - lgamma(theta$a_k)  # (K x 1)
    for (k in 1:K) {
        # print(class(theta$V_k[,,k]))
        e6_k[k] = e6_k[k] + log(det(as.matrix(theta$V_k[,,k])))
    }
    e_ln_q_beta_tau = sum(e6_k) - 0.5 * K * D * (log(2 * pi) + 1)
    theta$e_ln_q_beta_tau = sum(e6_k) - 0.5 * K * D * (log(2 * pi) + 1) # comment out later
    
    # (7) E [ ln q(gamma) ] ----------------------------------------------------
    e7 = - 0.5 * K * D * (log(2 * pi) + 1)
    for (k in 1:K) {
        e7 = e7 + log(det(as.matrix(theta$Q_k[,,k])))
    }
    e_ln_q_gamma = 0.5 * e7
    theta$e_ln_q_gamma = e7 # comment out later
    # compute the elbo ---------------------------------------------------------
    vlb = e_ln_p_y + e_ln_p_z + e_ln_p_gamma + e_ln_p_beta_tau -
        e_ln_q_z - e_ln_q_gamma - e_ln_q_beta_tau
    
    theta$L[theta$curr] = vlb
    
    return(theta)
    
} # end vb_elbo() function
