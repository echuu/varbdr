

## updateFunctions.R --- implementation of the the (closed-form) variational
##                       updates for covariate-independent weights with
##                       sparsity assumption in the gaussian components




## update q(z)
rnkUpdate = function(prior, theta) {
    
    # extract design, response, and dimensions for local use
    X = prior$X                                 # (N x D) design matrix
    y = prior$y
    N = prior$N                                 # num of observations
    K = prior$K                                 # num of clusters
    D = prior$D                                 # dimension of covariates
    
    
    # initialize the 3 matrices used to compute r_nk
    r_nk       = matrix(0, N, K)
    log_r_nk   = matrix(0, N, K)
    log_rho_nk = matrix(0, N, K)
    
    # commonly computed quantities
    ak_bk = theta$a_k / theta$b_k
    X_mk  = X %*% theta$m_k        # (N x K)
    
    # precompute digamma functions
    dig_a = digamma(theta$a_k)
    dig_b = digamma(theta$b_k)
    
    # precompute E[ln pi_k] \in R^K
    eln_pi_k = digamma(theta$alpha_k) - digamma(sum(theta$alpha_k))
    
    for (n in 1:N) { # fill in ln_rho_nk row-wise
        
        xSigmax = numeric(K)
        
        for (k in 1:K) {
            xSigmax[k] = quadMult(X[n,], theta$Sigma_k[,,k])
        }
        
        log_rho_nk[n,] = - 0.5 * log(2 * pi) + 0.5 * (dig_a - dig_b) - 
            0.5 * ak_bk * ((y[n] - X_mk[n,])^2 + xSigmax) + 
            eln_pi_k
        
    } # end outer for()
    
    logZ     = apply(log_rho_nk, 1, log_sum_exp)  # log of normalizing constant
    log_r_nk = log_rho_nk - logZ                  # log of r_nk
    r_nk     = apply(log_r_nk, 2, exp)            # exponentiate to recover r_nk
    
    theta$log_rho_nk  = log_rho_nk      # (N x K)
    theta$log_r_nk    = log_r_nk        # (N x K)
    theta$r_nk        = r_nk            # (N x K)
    theta$N_k         = colSums(r_nk)   # (K x 1)
    
    return(theta)
    
} # end rnkUpdate() function




## q(pi)
wtUpdate = function(prior, theta) {
    
    theta$alpha_k = prior$alpha_0 + theta$N_k    # (K x 1)
    return(theta)    
    
} # end wtUpdate() function




## update q(beta_d, omega_d)
ssUpdate = function(prior, theta) {
    
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
    
    ak_bk = theta$a_k / theta$b_k
    
    
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
                # R_dj[k, j] = theta$tau_k[k] * crossprod(r_nk[,k], x_dj)
                R_dj[k, j] = ak_bk[k] * crossprod(r_nk[,k], x_dj)
            } 
            # print(paste("d = ", d, "; ", "update R_dj", sep = ''))
        } # end of R_dj computation
        
        # (1.2) compute U_d
        for (k in 1:K) {
            U_d[k, d] = ak_bk[k] * crossprod(r_nk[,k], xd_sq)
            # print(paste("d = ", d, "; ", "update U_d", sep = ''))
        } # end of U_d computation
        
        # (1.3) compute zeta_d
        for (k in 1:K) {
            zeta_d[k, d] = ak_bk[k] * crossprod(r_nk[,k], xd_y)
            # print(paste("d = ", d, "; ", "update zeta_d", sep = ''))
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
        
        # print(paste("d = ", d, "; ", "update Q_d", sep = ''))
        
        # (2.2) update m_d
        m_d[,d] = Q_d_inv[,,d] %*% eta_d[,d]
        # print(paste("d = ", d, "; ", "update m_d", sep = ''))
        
        # (2.3) update lambda_d : sigmoid of log(lambda_d / (1 - lambda_d))
        lambda_d[d] = sigmoid(prior$log_odds[d] + 0.5 * K * log(prior$xi_0) -
                                  0.5 * log(det(as.matrix(Q_d[,,d]))) + 
                                  0.5 * crossprod(m_d[,d], eta_d[,d]))
        
        print(paste("lambda_d = ", lambda_d[d], sep = ''))
        
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
    
} # end ssUpdate() function




## update q(tau)
precUpdate = function(prior, theta) {
    
    K = prior$K
    X = prior$X
    y = prior$y
    
    a_k = numeric(K)
    b_k = numeric(K)
    
    a_k = prior$a_0 + 0.5 * theta$N_k # (K x 1)
    
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
    
    b_k = b_0 + 0.5 * (b_k_tmp + xSigmax) # (K x 1)
    
    theta$a_k = a_k
    theta$b_k = b_k
    
    return(theta)
    
} # end precUpdate() function

