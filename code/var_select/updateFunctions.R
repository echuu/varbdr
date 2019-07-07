



# spike and slab update for D covariates

# output:
#        spikeSlabObj  :              obj containing parameters for s&s update
#              Q_d     : D x (K x K)  inverse of precision matrix
#              Q_d_inv : D x (K x K)  precision matrix for beta_d
#              m_d     : (K x D)      mean vector for beta_d stored col-wise
#              R_dj    : D x (K x K)  intermediate matrix for eta_d
#              U_d     : D x (K x K)  intermediate matrix for Q_d
#              zeta_d  : (K x D)      intermediate vector for eta_d
#              eta_d   : (K x D)      intermediate vector for m_d

spikeSlabUpdate = function(prior, theta) {
    
    # TODO: figure out if it's better to modify in place for final/intermediate
    # vectors/matricies; if pre-define everything, that's a lot of space needed.
    # For now, define/allocate as needed

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
    m_d       = matrix(0, K, D)                  # (K x D)
    
    
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
        
        ## (2) update variational parameters: Q_d, m_d
        
        # (2.1) update Q_d
        Q_d[,,d]     = I_K + prior$xi_0 * diag(U_d[,d])
        Q_d_inv[,,d] = solve(Q_d[,,d])
        
        print(paste("d = ", d, "; ", "update Q_d", sep = ''))
        
        # (2.2) update m_d
        m_d = Q_d_inv[,,d] %*% eta_d[,d]
        print(paste("d = ", d, "; ", "update m_d", sep = ''))
        
        # end of variatioal updates --------------------------------------------
        
    } # end outer loop
    
    # prepare output
    spikeSlabObj = list(Q_d = Q_d, Q_d_inv = Q_d_inv, R_dj = R_dj, 
                        U_d = U_d, zeta_d = zeta_d, eta_d = eta_d)
    
    return(spikeSlabObj)    
}


# precisionUpdate() : update parmaeters for q(tau) = Ga(tau | a_k, b_k)
#     input     : 
#         prior : list of prior (and misc.) parameters
#         theta : list of updated variational parameters 
#     output    : list containing shape, rate parameters for q(tau)
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
    
    return(list(a_k = a_k, b_k = b_k))
} # end of precisionUpdate() functio



weightUpdate = function(prior, theta) {
    
    alpha   = 0
    xi      = 0
    lambda  = 0
    phi     = 0
    V_k     = 0
    V_k_inv = 0
    mu_k    = 0
    eta_k   = 0  # consider renaming so avoid confusion - nu?
    
    
    gammaUpdate = list(alpha = alpha, xi = xi, lambda = lambda, phi = phi,
                       Q_k = Q_k, Q_k_inv = Q_k_inv, mu_k = mu_k, 
                       eta_k = eta_k)
    
}





