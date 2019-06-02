

# initVarParams.R -- approximation + sparsity assumption



# initVarParams() -- initialize the variational parameters

initVarParams = function(y, X, N, D, K, intercept = FALSE, max_iter, 
                         m_k = NULL, mu_k = NULL) {

    I_D   = diag(1, D)                 # (D X D)  identity matrix
    I_K   = diag(1, K)                 # (K x K)  identity matrix
    
    L     = rep(-Inf, max_iter)        # Store the variational lower bounds
    
    
    # variables for inverence
    beta_k    = matrix(0, D, K)        # (D x K) : each beta_k stored col-wise
    tau_k     = rep(1, K)              # (1 x K) : scale for precision, V_k
    gamma_k   = matrix(0, D, K)        # (D x K) : each gamma_k stored col-wise
    
    
    # variables included for sparsity assumption
        # TODO: the initial value for probability of inclusion might have
        # to be tuned later
    lambda_d = rep(0.5, D)            # (D x 1) : probability of inclusion
    
    
    # e-step parameters: 
    # (0) q(z_nk) = r_nk^{z_nk}; log_r_nk, N_k
    
    r_nk       = matrix(0, N, K)       # (N x K) : normalized responsibilities
    log_r_nk   = matrix(0, N, K)       # (N x K) : log(r_nk)
    N_k        = colSums(r_nk) + 1e-10 # (K x 1) : sum of wts for each cluster
    
    # m-step parameters:
    
    ## START variable selection parameters -------------------------------------
    
    # (1a) q(beta_d | w_d) = w_d N(beta_d | m_d, Q_d^{-1}) + 
    #                                    (1 - w_d) * delta_0 (beta_d)
    Q_d = array(I_K, c(K, K, D))       # D x (K x K) : array of cov matricies
                                       #               for each beta_d (K x 1)
    Q_d = array(I_K, c(K, K, D))       # D x (K x K)
    
    # quantities used to compute Q_d, m_d
    U_d    = array(I_K, c(K, K, D))    # D x (K x K)
    R_dj   = array(I_K, c(K, K, D))    # D x (K x K)
    zeta_d = matrix(0, D, K)           # (D x K)
    eta_d  = matrix(0, D, K)           # (D x K)
    
    
    # (1b) q(w_d) = Ber(w_d | lambda_d)
    omega_d  = rep(1, D)              # (D x 1) : indicator for inclusion
    
    ## END variable selection parameters ---------------------------------------
    
    
    
    # (2) q(gamma_k) ~ N(gamma_k | mu_k, Q_k^{-1})
    
    # (2a) : variational parameters for approx (Bouchard)
    #     xi     : initialization cannot be 0 because lambda involves 1/xi
    #     lambda : compute here since alpha is a function of lambda
    alpha   = rep(1, N)                # N x 1 : used to compute xi_{n,1:K}
    xi      = matrix(1, N, K)          # N x K : nth row used to compute alpha_n
    lambda  = lambda_xi(xi)            # N x K : matrix of lambda(xi)
    phi     = numeric(N)               # N x 1 : function of alpha, xi
    
    
    # (2b) variational parameters for q(gamma_1:K) ~ N(gamma_k | mu_k, Q_k^{-1})
    Q_k     = array(I_D, c(D, D, K))   # K x (D x D)
    Q_k_inv = array(I_D, c(D, D, K))   # K x (D x D) precisions for gamma_k
    eta_k   = matrix(0, D, K)          # D x K       Q_k_inv * eta_k = mu_k
    
    
    
    # (3) variational parameters for tau_k
    a_k = rep(1, K)                    # K x 1 : shape param for tau_k
    b_k = rep(1, K)                    # K x 1 : rate param for tau_k  
    
    
    # initialize m_d, mu_k
    if (is.null(m_d) && is.null(mu_k)) {
        
        m_d     = matrix(0, D, K)          # D x K       : mean of beta_d
        mu_k    = matrix(0, D, K)          # D x K       : mean of gamma_k
        
        y_kmeans  = kmeans(y, K, nstart = 25)
        y_k_index = y_kmeans$cluster       # cluster index
        
        for (k in 1:K) {                   # mle of beta_k --> m_k
            y_k       = y[y_k_index == k]  # response vec for the k-th cluster
            X_k       = X[y_k_index == k,] # design mat for the k-th cluster                   
            beta_mle  = solve(t(X_k) %*% X_k, t(X_k) %*% y_k)
            m_d[,k]   = beta_mle
            mu_k[,k]  = beta_mle
        }
        
        
    } else if (is.null(m_d) || is.null(mu_k)) {
        warning("Both m_d, mu_k values must be provided, or both left NULL.")
    } else {
        
        # TODO: implement random initialization (can't just be around 0 right?)
        
        print("Using random initialization of m_k, mu_k")
    }
        
    
    # current iteration of CAVI: start at 2, because we look at *previous* elbo
    #                            elements when checking for convergence
    curr = 2 
    
    # list containing all variational parameters
    theta = list(beta_k = beta_k, tau_k = tau_k, gamma_k = gamma_k, 
                 omega_d, lambda_d,
                 r_nk = r_nk, log_r_nk = log_r_nk, N_k = N_k, 
                 Q_d = Q_d, Q_d_inv = Q_d_inv, m_d = m_d,
                 U_d = U_d, R_dj = R_dj, zeta_d = zeta_d, eta_d = eta_d,
                 alpha = alpha, xi = xi, lambda = lambda, phi = phi,
                 Q_k = Q_k, Q_k_inv = Q_k_inv, mu_k = mu_k, 
                 eta_k = eta_k, 
                 a_k = a_k, b_k = b_k, 
                 L = L, curr = curr, intercept = intercept)
    
} # end initVarParams() function


# end of initiVarParams.R file
