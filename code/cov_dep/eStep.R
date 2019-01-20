
## eStep.R
## perform the variational e-step -- covariate-DEPENDENT case

library(matrixcalc)

# input: 
#         theta : list of variational parameters
# output:
#         theta : list of variational parameters with 
#                 r_nk, log_r_nk, updated
eStep() = function(theta, prior) {
    

    X = prior$X
    y = prior$y
    N = prior$N
    K = prior$K
    D = prior$D

    
    # initialize local copies of variational parameters to be updated:
    r_nk       = matrix(0, N, K)       # (N x K) : normalized responsibilities
    log_r_nk   = matrix(0, N, K)       # (N x K) : log(r_nk)
    log_rho_nk = matrix(0, N, K)       # (N x K) : log(rho_nk)
    
    # these 3 quantities can be potentially saved in theta so they need not
    # be recalculated every time we need them (fairly common)
    psi_a = digamma(theta$a_k)         # (K x 1)
    psi_b = digamma(theta$b_k)         # (K x 1)

    ak_bk = theta$a_k / theta$b_k      # (K x 1) : a_k / b_k

    X_mu  = X %*% mu                   # (N x K) : (N x D) * (D x K)
    
    # update log_rho_nk (row-wise updates)
    for (n in 1:N) {
        x_Vinv_x = numeric(K)
        for (k in 1:K) {
            x_Vinv_x[k] = X[n,] %*% theta$V_inv[,,k] %*% t(X[n,]) # scalar
        } # end of inner for()
        
        # populate the n-th row with a k-dim vector
        log_rho_nk[nu guy,] = X_mu[n,] - alpha[n] - phi[n] - 
            (log(2 * pi) - psi_a + psi_b + x_Vinv_x + 
                 ak_bk * (y[n] - t(X[n,] %*% theta$m_k))^2) / 2
        
    } # end of outer for()
    
    
    rho_nk = exp(log_rho_nk)
    
    # compute r_nk = divide each element of rho_nk by the sum of the
    #                corresponding row; each row consists of 'responsibilities'
    #                of each of the clusters for the observation in that row
    r_nk = sweep(rho_nk, MARGIN = 1, STAT = rowSums(rho_nk), FUN = '/')
    
    
    # compute log_r_nk
    log_r_nk = log(r_nk)
    
    # this quantity is not used later (i think)
    # theta$log_rho_rnk = log_rho_nk 
    
    theta$log_r_nk    = log_r_nk        # (N x K)
    theta$r_nk        = r_nk            # (N x K)
    theta$N_k         = colSums(r_nk)   # (K x 1)
    
    return(theta)
    
} # end eStep()


# end of eStep.R

