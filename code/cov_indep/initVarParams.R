
## initVarParams.R -- covariate-INDEPENDENT case

# initialize the variational parameters

# input:
#         N        : (1 x 1) -- # of response values
#         D        : (1 x 1) -- dimension of covariates
#         K        : (1 x 1) -- # of mixture componenets
# output: 
#         theta    : list containing the variational parameters

initVarParams = function(y, X, N, D, K, max_iter) {

    I_D   = diag(1, D)                 # (D X D)  identity matrix
    L     = rep(-Inf, max_iter)        # store the variational lower bonuds
    dc    = vector("list", max_iter)   # store density curves for each iter
    
    # explicit random variables -- ---------------------------------------------
    #    don't thnk these are used in CAVI, these are used later to 
    #    generate the model parameters
    pi_k      = numeric(K)             # (1 x K) : concentration parameters
    beta_k    = matrix(1, D, K)        # (D x K) : each beta_k stored col-wise
    tau_k     = rep(1, K)              # (1 x K) : scale for precision, V_k


    # E[z_nk] = r_nk; these quantities are updated during the var e-step
    #     Note: rho_nk / sum_j rho_nj
    #     can calculate this by exponentiating log_rho_nk, 
    #     dividing each element by the sum of the ROW it belongs in
    r_nk       = matrix(0, N, K)       # (N x K) : normalized responsibilities
    log_r_nk   = matrix(0, N, K)       # (N x K) : log(r_nk)
    N_k        = colSums(r_nk) + 1e-10 # (K x 1) : sum of wts for each cluster
    
    # model parameters for the random variables involved in the model
    
    # (1) variational parameter for pi -----------------------------------------
    # q(pi) = Dir( alpha_1, ... , alpha_K )
    alpha_k = rep(1 / K, K)            # K x 1
    
    # (2) variational parameters for beta_k | tau_k ----------------------------
    # q ( beta_k | tau_k ) ~ N ( m_k, (tau_k V_k)^{-1} )
    V_k     = array(I_D, c(D, D, K))   # K x (D x D)
    V_k_inv = array(I_D, c(D, D, K))   # K x (D x D) : scaled precision
    zeta_k  = matrix(0, D, K)          # (D x K)     : V_k_inv * zeta_k = m_k
    m_k     = matrix(0, D, K)        # (D x K)     : mean of gamma_k

    # (3) variational parameters for tau_k -------------------------------------
    # q(tau_k) = Ga ( a_k, b_k )
    a_k = rep(1, K)                    # (K x 1) : shape param for tau_k
    b_k = rep(1, K)                    # (K x 1) : rate param for tau_k
    

    # previous initialization doesn't work (as seen when we try to do something 
    # similar in the GMM model with mean params initialized to all 0)
    # instead: we replace m_k with the mle estimate for beta_k
    # k is determined by first doing k-means on the y-values
    
    y_kmeans  = kmeans(y, K, nstart = 25)
    y_k_index = y_kmeans$cluster                      # cluster index
    for (k in 1:K) {                                  # mle of beta_k --> m_k
        y_k       = y[y_k_index == k]
        X_k       = X[y_k_index == k,]                   
        beta_mle  = solve(t(X_k) %*% X_k, t(X_k) %*% y_k)
        m_k[,k]   = beta_mle
    }

    # current iteration of CAVI
    curr = 1
    
    # create object with all variational parameters ----------------------------

    theta = list(pi_k = pi_k, beta_k = beta_k, tau_k = tau_k,    # RVs
                 log_r_nk = log_r_nk, r_nk = r_nk, N_k = N_k,    # r_nk's
                 alpha_k = alpha_k,                              # dir params
                 V_k_inv = V_k_inv, zeta_k = zeta_k, m_k = m_k,  # gaus. params
                 a_k = a_k, b_k = b_k,                           # gam. params
                 L = L, curr = curr)                             # elbo, curr
    
    
} # end initVarParam() function
    