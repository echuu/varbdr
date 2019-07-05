





# spike and slab update for D covariates

# output:
#        spikeSlabObj  : object containing parameters for s&s update
#              Q_d     : D x (K x K) inverse of precision matrix
#              Q_d_inv : D x (K x K) precision matrix for beta_d
#              m_d     : (K x D) mean vector for beta_d stored col-wise
#              R_dj    : D x (K x K) intermediate matrix for eta_d
#              U_d     : D x (K x K) intermediate matrix for Q_d
#              zeta_d  : (K x D) intermediate vector for eta_d
#              eta_d   : (K x D) intermediate vector for m_d

spikeSlabUpdate = function(prior, theta) {
    
    # TODO: figure out if it's better to modify in place for final/intermediate
    # vectors/matricies; if pre-define everything, that's a lot of space needed.
    # For now, define/allocate as needed

    # dimensions
    N = prior$N
    K = prior$K
    D = prior$D
    
    # commonly used matrices, vectors (more to be added as necessary)
    
    I_K       = diag(1, K)                       # (K x K) : identity matrix
    
    
    # define final matricies, vectors
    Q_d       = array(I_K, c(K, K, D))           # D x (K x K): inv precision
    Q_d_inv   = array(I_K, c(K, K, D))           # D x (K x K): precision 
    m_d       = matrix(0, K, D)                  # (K x D)
    
    
    # allocate space for intermediate matrics
    R_dj   = array (I_K, c(K, K, D))             # D x (K x K)
    U_d    = array(I_K, c(K, K, D))              # D x (K x K)
    zeta_d = matrix(0, K, D)                     # (K x D)
    eta_d  = matrix(0, K, D)                     # (K x D)
    
    
    for (d in 1:D) {
        ## compute intermediate quantities:
        
        # (1) compute U_d, zeta_d, R_dj
        
        # (2) compute eta_d
        
        ## update variational parameters
        
        # (3) update Q_d
        Q_d[,,d] = I_K + prior$xi_0 * U_d[,,d]
        
        # (4) update m_d
        m_d = Q_d[,,d] %*% eta_d[,d]
    }
    
    # prepare output
    spikeSlabObj = list(Q_d, Q_d_inv, R_dj, U_d, zeta_d, eta_d)
    
    return(spikeSlabObj)    
}






