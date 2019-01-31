
# initVarParams.R -- covariate-DEPENDENT case

   

# initVarParams() -- initialize the variational parameters
    # beta_k  : matrix of all 1's (not 0 since we start with X * beta != 0)
    # tau_k   : vector of 1's
    # gamma_k : matrix of all 0's
    # V_k_inv : (K x K) identity matrix
    # m_k     : matrix of 0-mean vectors for beta_k
    # a_k     : vector of all 1's
    # b_k     : vector of all 1's
    # alpha   : vector of all 1's
    # xi      : matrix of all 0's (these will be re-estimated using alphas)
    # Q_k_inv : (K x K) identity matrix
    # mu_k    : matrix of 0-mean vectors for gamma_k
initVarParams = function(y, X, N, D, K, max_iter) {
    
    I_D   = diag(1, D)                 # (D X D)  identity matrix
    L     = rep(-Inf, max_iter)        # Store the variational lower bounds
 
    
    
    # explicit random variables -- 
    #    don't thnk these are used in CAVI, these are used later to 
    #    generate the model parameters
    # pi_k      = rep(1 / K, K)        # (1 x K) : mixing weights
    beta_k    = matrix(0, D, K)        # (D x K) : each beta_k stored col-wise
    tau_k     = rep(1, K)              # (1 x K) : scale for precision, V_k
    gamma_k   = matrix(0, D, K)        # (D x K) : each gamma_k stored col-wise
    
    # E[z_nk] = r_nk; these quantities are updated during the var e-step
        # Note: rho_nk / sum_j rho_nj
        # can    calculate this by exponentiating log_rho_nk, 
        # dividing each element by the sum of the ROW it belongs in
    
    r_nk       = matrix(0, N, K)       # (N x K) : normalized responsibilities
    log_r_nk   = matrix(0, N, K)       # (N x K) : log(r_nk)
    N_k        = colSums(r_nk) + 1e-10 # (K x 1) : sum of wts for each cluster
    
    # model parameters for each of the random variables above
    
    # (1) variational parameters for beta_k | tau_k ~ N(m_k, (tau_k V_k)^{-1})
    V_k     = array(I_D, c(D, D, K))   # K x (D x D)
    V_k_inv = array(I_D, c(D, D, K))   # K x (D x D) : scaled precision
    zeta_k  = matrix(0, D, K)          # D x K       : V_k_inv * zeta_k = m_k
    m_k     = matrix(0, D, K)          # D x K       : mean of gamma_k
    
    # (2) variational parameters for tau_k
    a_k = rep(1, K)                    # K x 1 : shape param for tau_k
    b_k = rep(1, K)                    # K x 1 : rate param for tau_k
    
    # (3) variational parameters for gamma_k (via Bouchard)
    alpha   = rep(1, N)                # N x 1 : used to compute xi_{n,1:K}
    xi      = matrix(0, N, K)          # N x K : nth row used to compute alpha_n
    lambda  = matrix(0, N, K)          # N x K : matrix of lambda(xi)
    phi     = numeric(N)               # N x 1 : function of alpha, xi
        
    # variational parameters for gamma_1:K ~ N(gamma_k | mu_k, Q_k^{-1})
    Q_k     = array(I_D, c(D, D, K))   # K x (D x D)
    Q_k_inv = array(I_D, c(D, D, K))   # K x (D x D) precisions for gamma_k
    eta_k   = matrix(0, D, K)          # D x K       Q_k_inv * eta_k = mu_k
    mu_k    = matrix(0, D, K)          # D x K       mean of gamma_k
    
    
    # provide better starting values for m_k, mu_k
    # for now, we use same starting values for both m_k, mu_k
    # same as in cov-indpt case: replace m_k with mle estimate for beta_k
    # k is determined by first doing k-means on the y-values
    
    y_kmeans  = kmeans(y, K, nstart = 25)
    y_k_index = y_kmeans$cluster       # cluster index
    for (k in 1:K) {                   # mle of beta_k --> m_k
        y_k       = y[y_k_index == k]  # response vector for the k-th cluster
        X_k       = X[y_k_index == k,] # design matrix for the k-th cluster                   
        beta_mle  = solve(t(X_k) %*% X_k, t(X_k) %*% y_k)
        m_k[,k]   = beta_mle
        mu_k[,k]  = beta_mle
    }
    
    # current iteration of CAVI
    curr = 0
    
    # list containing all variational parameters
    theta = list(beta_k = beta_k, tau_k = tau_k, gamma_k = gamma_k,
                 r_nk = r_nk, log_r_nk = log_r_nk, N_k = N_k, 
                 V_k = V_k, V_k_inv = V_k_inv, zeta_k = zeta_k, m_k = m_k, 
                 a_k = a_k, b_k = b_k, 
                 alpha = alpha, xi = xi, lambda = lambda, phi = phi,
                 Q_k = Q_k, Q_k_inv = Q_k_inv, eta_k, mu_k = mu_k, 
                 L = L, curr = curr)
    
    return(theta)
    
} # end of initVarParams() function



# end of initVarParams.R
