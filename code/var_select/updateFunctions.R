



# spike and slab update for D covariates
#     input            : 
#                prior : list of prior (and misc.) parameters
#                theta : list of current values of variational parameters 
#     output           : theta, list of variational params w/ following updated
#              Q_d     : D x (K x K)  inverse of precision matrix
#              Q_d_inv : D x (K x K)  precision matrix for beta_d
#              m_d     : (K x D)      mean vector for beta_d stored col-wise
#              R_dj    : D x (K x K)  intermediate matrix for eta_d
#              U_d     : D x (K x K)  intermediate matrix for Q_d
#              zeta_d  : (K x D)      intermediate vector for eta_d
#              eta_d   : (K x D)      intermediate vector for m_d
spikeSlabUpdate = function(prior, theta) {

    # dimensions
    N = prior$N
    K = prior$K
    D = prior$D
    
    # design matrix, response vector
    X = prior$X
    y = prior$y
    
    
    # commonly used matrices, vectors (more to be added as necessary)
    I_K       = diag(1, K)                       # (K x K) : identity matrix
    r_nk      = theta$r_nk
    
    # define final matricies, vectors
    Q_d       = array(I_K, c(K, K, D))           # D x (K x K): inv precision
    Q_d_inv   = array(I_K, c(K, K, D))           # D x (K x K): precision 
    m_d       = matrix(0, K, D)                  # (K x D) : mean comps 
    lambda_d  = numeric(D)                       # (D x 1) : inclusion probs.
    
    # allocate space for intermediate matrics
    
    # R_dj is a (K x D) matrix unique to the d-th covariate; for each d in 1:D, 
    # we need to compute (D-1) diagonal matrices, R_dj
    # We store the K diagonal entries in each of the j columns; (d-th col empty)
    R_dj   = matrix(0, K, D)
    
    # U_d is a (K x D) matrix, where the d-th col stores the K diagonal entries 
    # for the d-th covariate
    U_d    = matrix(0, K, D)
    
    # zeta_d is a (K x D) matrix, where the d-th column corresponds to the 
    # d-th covariate -- used for intermediate calculation
    zeta_d = matrix(0, K, D)
    
    # zeta_d is a (K x D) matrix, where the d-th column corresponds to the 
    # d-th covariate -- used for intermediate calculation
    eta_d  = matrix(0, K, D)
    
    
    for (d in 1:D) {
        
        j_index = c(1:D)[-d] # 1 <= j <= D, j != d
        
        xd_sq = X[,d]^2
        xd_y  = X[,d] * y 
        
        ## compute intermediate quantities:
        
        # (1) compute R_dj, U_d, zeta_d, eta_d
        
        ## (1.1) compute R_dj
        for (j in j_index) {
            
            x_dj = X[,d] * X[,j] # (N x 1)        
            
            for (k in 1:K) { # compute elements of the R_dj matrix
                R_dj[k, j] = theta$tau_k[k] * crossprod(r_nk[,k], x_dj)
            } 
            print(paste("d = ", d, "; ", "update R_dj", sep = ''))
        } # end of R_dj computation
        
        # (1.2) compute U_d
        for (k in 1:K) {
            U_d[k, d] = theta$tau_k[k] * crossprod(r_nk[,k], xd_sq)
            print(paste("d = ", d, "; ", "update U_d", sep = ''))
        } # end of U_d computation
        
        # (1.3) compute zeta_d
        for (k in 1:K) {
            zeta_d[k, d] = theta$tau_k[k] * crossprod(r_nk[,k], xd_y)
            print(paste("d = ", d, "; ", "update zeta_d", sep = ''))
        } # end of zeta_d computation
        
        # (1.4) compute eta_d
        eta_tmp = matrix(0, K, 1)
        for (j in j_index) {
            eta_tmp = eta_tmp + diag(R_dj[,j]) %*% theta$m_d[,j]
        }
        
        eta_d[,d] = zeta_d[,d] - 0.5 * eta_tmp
        
        # end of intermediate calculations -------------------------------------
        
        ## (2) update variational parameters: Q_d, m_d, lambda_d
        
        # (2.1) update Q_d
        Q_d[,,d]     = I_K + prior$xi_0 * diag(U_d[,d])
        Q_d_inv[,,d] = solve(Q_d[,,d])
        
        print(paste("d = ", d, "; ", "update Q_d", sep = ''))
        
        # (2.2) update m_d
        m_d = Q_d_inv[,,d] %*% eta_d[,d]
        print(paste("d = ", d, "; ", "update m_d", sep = ''))
        
        # (2.3) update lambda_d : sigmoid of log(lambda_d / (1 - lambda_d))
        lambda_d[d] = sigmoid(prior$log_odds[d] + K / 2 * log(prior$xi_0) + 
            0.5 * crossprod(m_d[,d], eta_d[,d]))
        
        # end of variational updates -------------------------------------------
        
    } # end outer loop
    
    # primary quantities
    theta$Q_d      = Q_d
    theta$Q_d_inv  = Q_d_inv      # V [ beta_d | omega_d = 1 ] = Q_d^(-1)
    theta$m_d      = m_d          # E [ beta_d | omega_d = 1 ] = m_d
    theta$lambda_d = lambda_d     # posterior inclusion probabilities
    
    # secondary quantities (may come up in ELBO later)
    theta$R_dj   = R_dj
    theta$U_d    = U_d
    theta$zeta_d = zeta_d
    theta$eta_d  = eta_d
    
    return(theta)    
} # end spikeSlabUpdate() function


# precisionUpdate() : update parmaeters for q(tau) = Ga(tau | a_k, b_k)
#     input     : 
#         prior : list of prior (and misc.) parameters
#         theta : list of current values of variational parameters 
#     output    : theta, list of variational params w/ following updated
#         a_k   : (K x 1) vector of shape parameters for the K precision comps
#         b_k   : (K x 1) vector of rate parameters for the K precision comps
precisionUpdate = function(prior, theta) {
    
    K = prior$K
    X = prior$X
    y = prior$y
    
    a_k = numeric(K)
    b_k = numeric(K)
    
    a_k = prior$a_0 + 0.5 * theta$N_k - 1
    
    X_mk    = X %*% theta$m_k  # (N x K) : (N x D) * (D x K)
    b_k_tmp = numeric(K)       # (K x 1) : 1st term in summation
    xSigmax = numeric(K)       # (K x 1) : 2nd term in summation
    
    for (k in 1:K) {
         b_k_tmp[k] = b_k_tmp[k] + sum((y - X_mk[,k])^2)
    }
    
    for (k in 1:K) {
        for (n in 1:N) {
            xSigmax[k] = xSigmax[k] + t(X[n,]) %*% theta$Sigma_k[,,k] %*% X[n,]
        }
    }
    
    b_k = b_0 + 0.5 * (b_k_tmp + xSigmax)
    
    theta$a_k = a_k
    theta$b_k = b_k
    
    return(theta)
} # end of precisionUpdate() functio


# precisionUpdate() : update parmaeters for q(gamma) = N(gamma | mu_k, Q_k_inv)
#     input       : 
#         prior   : list of prior (and misc) parameters
#         theta   : list of current values of variational parameters 
#     output      : theta, list of variational params w/ following updated
#         alpha   : (N x 1)      additional var. param for upper bound
#         xi      : (N x K)      additional var. param for upper bound
#         lambda  : (N x K)      function of xi
#         phi     : (N x 1)      function of alpha, xi
#         V_k     : K x (D x D)  inverse of precision matrix for gamma_k
#         V_k_inv : K x (D x D)  precision matrix for gamma_k
#         mu_k    : (D x K)      mean vectors for gamma_k's, stored col-wise
#         eta_k   : (D x K)      intermediate vector for computing mu_k
weightUpdate = function(prior, theta) {
    
    X = prior$X
    y = prior$y
    N = prior$N
    K = prior$K
    D = prior$D
    
    I_D        = diag(1, D)                       # (D X D) : identity matrix
    X_mu       = X %*% theta$mu_k                 # (N x K) : (N x D) * (D x K)
    
    ## TODO: code below needs to be adapted to the variables used during v.s.
    
    ## update additional variational params used in Bouchard bound -------------
    #     (0.1) alpha                 # (N x 1)
    #     (0.2) xi                    # (N x K)
    #     (0.3) lambda                # (N x K)
    #     (0.4) phi                   # (N x 1)
    for (n in 1:N) {
        
        # (0.1) update alpha ---------------------------------------------------
        theta$alpha[n] = 1 / sum(theta$lambda[n,]) * 
            (0.5 * (0.5 * K - 1) + crossprod(X_mu[n,], theta$lambda[n,]))
        
        # (0.2) update xi (computed row-wise) ----------------------------------
        xVx = numeric(K)
        for (k in 1:K) {
            xVx[k] = quadMult(X[n,], theta$V_k_inv[,,k]) # function in misc.R
        }
        
        theta$xi[n,] = sqrt((X_mu[n,] - theta$alpha[n])^2 + xVx)
    }
    
    # (0.3) compute lambda(xi) using updated value of xi  ----------------------
    theta$lambda = lambda_xi(theta$xi)              # (N x K), misc.R
    
    # (0.4) compute phi (function of alpha, xi) --------------------------------
    for (n in 1:N) {
        theta$phi[n] = sum(0.5 * (X_mu[n,] - theta$alpha[n] - theta$xi[n,]) + 
                               log(1 + exp(theta$xi[n,])))
    } # end update for phi
    
    
    # end of Bouchard variational parameter updates ----------------------------
    
    
    ## update q(gamma) = N(gamma | mu_k, V_k^{-1} ) ----------------------------

    # (1.1) update q(gamma) : V_k, V_k^{-1}, eta_k, mu_k -----------------------
    for (k in 1:K) {
        
        r_x = 0
        rl_nk_xx = matrix(0, D, D)               # used to calculate q(gamma)

        for (n in 1:N) {
            ## TODO: save x_n * x_n' calculation since this is done every 
            # iteration of CAVI -> store N x (N x N) matrices (worth it?)
            # may be too expensive for large N
            rl_nk_xx = rl_nk_xx + 
                theta$r_nk[n,k] * theta$lambda[n,k] * crossprod(t(X[n,]))
        }
        
        
        # (1.1) update q(gamma) : V_k, V_k^{-1}, eta_k, mu_k ------------------
        theta$V_k[,,k] = I_D + 2 * rl_nk_xx
        theta$V_k_inv[,,k] = solve(theta$V_k[,,k])
        
        theta$eta_k[,k] = t(X) %*% 
            (theta$r_nk[,k] * (0.5 + 2 * theta$lambda[,k] * theta$alpha))
        
        theta$mu_k[,k] =  theta$V_k_inv[,,k] %*% theta$eta_k[,k]
        
        
    } # end of q(gamma)
    
    # --------------------------------------------------------------------------
    
    return(theta)
} # end weightUpdate() function










