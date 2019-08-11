

# initVP.R -- initilaize variatonal parameters for covariate-independent weights
#             w/ variable selection implemented for gaussian componenets


initVarParams_indep = function(y, X, N, D, K, m_d = NULL,
                               intercept = FALSE, max_iter = 1e4, 
                               tol = 1e-3, VERBOSE = TRUE) {
    
    I_D   = diag(1, D)                 # (D X D)  identity matrix
    I_K   = diag(1, K)                 # (K x K)  identity matrix
    
    
    pip_D = rep(0.3, D)                # (D x 1)  prob of inclusion 
    
    
    # variational lower bound
    L = rep(-Inf, max_iter)
    
    # random variables variables for inference
    pi_k      = numeric(K)             # (1 x K) : conceptration parameters
    beta_k    = matrix(0, D, K)        # (D x K) : each beta_k stored col-wise
    tau_k     = rep(1, K)              # (1 x K) : scale for precision, V_k
    omega_d   = rep(1, D)              # (1 x D) : indicator for inclusion
    
    
    #  ----------------------------------------------------------------------  #
    ## e-step update parameters 
    
    # (0) q(z_nk)
    r_nk      =  matrix(0, N, K)       # (N x K) : normalized responsibilities
    log_r_nk  = matrix(0, N, K)
    N_k       = colSums(r_nk) + 1e-10
    
    
    #  ----------------------------------------------------------------------  #
    ## m-step update parameters 
    
    # (1) q(pi) = \prod_k q(pi_k) = Dir(pi | alpha_1,...,alpha_K)
    # pi = (pi_1, pi_2, ... , pi_K)
    alpha_k = rep(1 / K, K)            # (K x 1)
        
    # (2a) q(beta_d | omega_d) = omega_d N( beta_d | m_d, Q_d^{-1}) + 
    #                              (1 - omega_d) delta_0 (beta_d)
    
    Q_d = array(I_K, c(K, K, D))       # D x (K x K) : array of cov matricies
                                       #               for each beta_d (K x 1)
    Q_d_inv = array(I_K, c(K, K, D))   # D x (K x K)
    
    # quantities used to compute Q_d, m_d 
    U_d    = matrix(0, K, D)           # (K x D) : d-th col reps diag K x K mat
    R_dj   = matrix(0, K, D)           # (K x D) : d-th col reps diag K x K mat
    zeta_d = matrix(0, K, D)           # (K x D) : d-th col reps diag K x K mat
    eta_d  = matrix(0, K, D)           # (K x D) : d-th col reps diag K x K mat
    
    
    # (2b) q(omega_d) = Ber(omega_d | lambda_d)
    lambda_d = pip_D                   # (D x 1) : inclusion probabilities
    omega_d  = rbinom(D, 1, lambda_d)  # (D x 1) : generate random indicators
    
    # (3) q(tau_k) = Gamma(tau_k | a_k, b_k)
    a_k = rep(1, K)                    # K x 1 : shape param for tau_k
    b_k = rep(1, K)                    # K x 1 : rate param for tau_k 
    
    # initialize m_d
    if (is.null(m_d)) {
        
        m_d       = matrix(0, K, D)            # (K x D) : mean of beta_d
        y_kmeans  = kmeans(y, K, nstart = 25)
        y_k_index = y_kmeans$cluster           # cluster index
        
        for (k in 1:K) {                       # mle of beta_k =: m_k
            y_k       = y[y_k_index == k]      # response vec for k-th cluster
            X_k       = X[y_k_index == k,]     # design mat for k-th cluster                   
            beta_mle  = solve(t(X_k) %*% X_k, t(X_k) %*% y_k)
            m_d[k,]   = beta_mle
        }
        
        
    } else {
        # TODO: implement random initialization (can't just be around 0 right?)
        print("Using random initialization of m_k, mu_k")
    }
    
    # m_k := E(beta_k) \in R^D
    # Sigma_k := Cov(beta_k) \in R^{D x D}
    Sigma_k = array(I_D, c(D, D, K)) # K x (D x D) : cov matrix for each beta_k
    
    # TODO: this should not just be the transpose of m_d because m_d is the
    # conditional mean of beta_d, whereas m_k is the unconditional mean 
    # of beta_k
    m_k     = t(m_d)                 # (D x K)     : mean components for beta_k, 
                                     #               stored col-wise
    
    # current iteration of CAVI: start at 2, because we look at *previous* elbo
    #                            elements when checking for convergence
    curr = 2
    
    # convergence status
    converge = FALSE
    
    # list containing all variational parameters
    theta = list(beta_k = beta_k, tau_k = tau_k, pi_k = pi_k, 
                 omega_d = omega_d, lambda_d = lambda_d,
                 r_nk = r_nk, log_r_nk = log_r_nk, N_k = N_k, 
                 Q_d = Q_d, Q_d_inv = Q_d_inv, m_d = m_d,
                 U_d = U_d, R_dj = R_dj, zeta_d = zeta_d, eta_d = eta_d,
                 alpha_k = alpha_k, 
                 a_k = a_k, b_k = b_k, 
                 Sigma_k = Sigma_k, m_k = m_k,
                 L = L, curr = curr, intercept = intercept, max_iter = max_iter,
                 tol = tol, converge = converge, VERBOSE = VERBOSE)
}



